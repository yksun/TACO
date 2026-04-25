"""Pipeline runner for TACO genome assembly."""
import os
import sys
import subprocess
import shutil
import gzip
import socket
import csv
import glob
import hashlib
import json
import shlex
import time
import platform as plat
from datetime import datetime

from taco.cli import STEP_NAMES


class TeeWriter:
    """Write to both a file and the original stream (tee-like behaviour)."""

    def __init__(self, log_fh, original_stream):
        self.log_fh = log_fh
        self.original = original_stream

    def write(self, data):
        if data:
            self.original.write(data)
            self.log_fh.write(data)

    def flush(self):
        self.original.flush()
        self.log_fh.flush()

    # Forward any attribute lookups to the original stream so that
    # code that checks e.g. sys.stdout.isatty() still works.
    def __getattr__(self, name):
        return getattr(self.original, name)


class PipelineRunner:
    """Main pipeline execution engine for TACO."""

    PIPELINE_NAME = "TACO-1.3.1"

    def __init__(self, args):
        # Core parameters
        self.genomesize = args.genomesize
        self.threads = args.threads
        self.fastq = os.path.realpath(args.fastq)
        self.motif = args.motif
        self.taxon = getattr(args, 'taxon', 'other')
        self.platform = args.platform
        self.reference_fasta = args.reference
        self.steps = args.steps
        self.assembly_only = args.assembly_only
        self.benchmark = getattr(args, 'benchmark', False)

        # Telomere parameters — taxon-aware defaults
        # Fungi: smaller windows (telomere arrays are short, 50-300 bp)
        # Plants/vertebrates: larger windows (telomere arrays can span 5-20 kb)
        self.telomere_mode = getattr(args, 'telomere_mode', 'hybrid')
        self.telo_end_window = getattr(args, 'telo_end_window', 5000)
        self.telo_kmer_min = getattr(args, 'telo_kmer_min', 4)
        self.telo_kmer_max = getattr(args, 'telo_kmer_max', 30)

        # Score window: user override takes precedence; otherwise taxon-aware
        user_score_window = getattr(args, 'telo_score_window', None)
        if user_score_window is not None and user_score_window != 500:
            self.telo_score_window = user_score_window
        elif self.taxon == "fungal":
            self.telo_score_window = 300   # fungal telomere arrays are short
        elif self.taxon in ("plant", "vertebrate", "animal"):
            self.telo_score_window = 1000  # longer arrays in these taxa
        else:
            self.telo_score_window = 500   # balanced default

        # Post-refinement options
        self.no_purge_dups = getattr(args, 'no_purge_dups', False)
        self.no_polish = getattr(args, 'no_polish', False)
        self.no_coverage_qc = getattr(args, 'no_coverage_qc', False)
        self.allow_t2t_replace = getattr(args, 'allow_t2t_replace', False)

        # Backbone selection
        self.auto_mode = args.auto_mode
        self.assembler = None
        self.choose_flag = False
        if args.choose is not None:
            self.choose_flag = True
            if args.choose != "__prompt__":
                self.assembler = args.choose

        # BUSCO
        self.busco_lineage = args.busco  # may be None if --taxon other and no --busco
        self.run_busco = args.busco is not None
        if self.busco_lineage is None and self.taxon != "other":
            # Fallback: use taxon default (should be set in cli.py)
            from taco.cli import TAXON_BUSCO_LINEAGE
            self.busco_lineage = TAXON_BUSCO_LINEAGE.get(self.taxon, "ascomycota_odb10")
            self.run_busco = True

        # Merqury: explicitly disabled with --no-merqury; explicitly enabled
        # with --merqury or --merqury-db.  Otherwise, auto-enable when an
        # existing .meryl database is discoverable or when merqury.sh + meryl
        # are installed so the read database can be built automatically.
        self.merqury_db = getattr(args, 'merqury_db', None)
        mk = getattr(args, 'merqury_k', "auto")
        if mk == "auto":
            self.merqury_k = "auto"
        else:
            try:
                self.merqury_k = int(mk)
            except (ValueError, TypeError):
                self.merqury_k = 21
        self.merqury_build_db = False
        if getattr(args, 'no_merqury', False):
            self.merqury_enable = False
        elif args.merqury or self.merqury_db:
            self.merqury_enable = True
        else:
            # Auto-detect: check if merqury.sh is installed and a .meryl db exists
            import shutil as _shutil
            merqury_bin = _shutil.which("merqury.sh")
            auto_db = None
            if merqury_bin:
                for cand in ["reads.meryl", "meryl/reads.meryl",
                             "merqury/reads.meryl"]:
                    if os.path.isdir(cand):
                        auto_db = cand
                        break
                if auto_db is None:
                    import glob as _glob
                    found = _glob.glob("*.meryl")
                    if found:
                        auto_db = found[0]
            if merqury_bin and auto_db:
                self.merqury_enable = True
                self.merqury_db = auto_db
            else:
                self.merqury_enable = False

        # Auto-enable Merqury for ALL platforms when merqury.sh + meryl are
        # installed.  Merqury is most accurate with high-accuracy reads, but
        # ONT/CLR can still provide a relative comparison signal.
        if not self.merqury_enable:
            merqury_bin = shutil.which("merqury.sh")
            meryl_bin = shutil.which("meryl")
            if merqury_bin and meryl_bin:
                self.merqury_enable = True
                self.merqury_build_db = True
                self.log_info(f"Merqury auto-enabled for {self.platform} "
                              f"(will build a reads .meryl database from input reads)")
            elif merqury_bin:
                self.log_warn("merqury.sh found but meryl not installed — "
                              "Merqury disabled")
        # Warn for non-high-accuracy primary read platforms.
        if self.merqury_enable and self.platform != "pacbio-hifi":
            platform_label = {
                "nanopore": "Oxford Nanopore",
                "pacbio": "PacBio CLR",
            }.get(self.platform, self.platform)
            self.log_warn(
                f"Merqury QV estimates are most reliable with high-accuracy "
                f"reads (PacBio HiFi or Illumina). With {platform_label} reads, "
                f"QV values may underestimate true assembly quality. "
                f"Merqury completeness and relative QV ranking across "
                f"assemblers should be interpreted cautiously. "
                f"Disable with --no-merqury.")

        # Derived paths
        fastq_name = os.path.basename(self.fastq)
        self.project = fastq_name.replace('.fastq.gz', '').replace('.fq.gz', '').replace('.fastq', '').replace('.fq', '')
        self.taco_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.working_dir = os.getcwd()

        # Logging and optional benchmarking
        self.logs_dir = os.path.join(self.working_dir, "logs")
        self.benchmark_dir = os.path.join(self.working_dir, "benchmark_logs")
        os.makedirs(self.logs_dir, exist_ok=True)
        self.run_start_ts = datetime.now()
        self.run_end_ts = None
        self.run_start_monotonic = None
        self.run_id = self.run_start_ts.strftime("%Y%m%d_%H%M%S")
        self.bench_step_tsv = os.path.join(self.benchmark_dir, "step_benchmark.tsv")
        self.bench_run_tsv = os.path.join(self.benchmark_dir, "run_metadata.tsv")
        self.bench_summary = os.path.join(self.benchmark_dir, "run_summary.txt")
        self.bench_manifest_json = os.path.join(self.benchmark_dir, "run_manifest.json")
        self.bench_tools_tsv = os.path.join(self.benchmark_dir, "software_versions.tsv")
        self.bench_outputs_tsv = os.path.join(self.benchmark_dir, "output_manifest.tsv")
        self.bench_methods_txt = os.path.join(self.benchmark_dir, "methods_note.txt")
        self._sha256_cache = {}

    # ------------------------------------------------------------------ #
    # Logging
    # ------------------------------------------------------------------ #
    @staticmethod
    def _ts():
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def log(self, msg):
        print(f"[{self._ts()}] {msg}", flush=True)

    def log_info(self, msg):
        print(f"[info] {msg}", flush=True)

    def log_warn(self, msg):
        print(f"[warn] {msg}", file=sys.stderr, flush=True)

    def log_error(self, msg):
        print(f"[error] {msg}", file=sys.stderr, flush=True)

    # ------------------------------------------------------------------ #
    # Command execution
    # ------------------------------------------------------------------ #
    def run_cmd(self, cmd, desc=None, check=True):
        """Run an external command. cmd may be a string (shell) or list."""
        if desc:
            self.log(desc)
        if isinstance(cmd, list):
            display = " ".join(str(c) for c in cmd)
        else:
            display = cmd
        self.log(f"$ {display}")
        result = subprocess.run(
            cmd,
            shell=isinstance(cmd, str),
            capture_output=False,
        )
        if check and result.returncode != 0:
            self.log_error(f"Command failed (exit {result.returncode}): {display}")
            sys.exit(result.returncode)
        return result

    def require_cmd(self, name):
        """Return True if *name* is on PATH, else print warning and return False."""
        if shutil.which(name):
            return True
        self.log_warn(f"Required command not found: {name}")
        return False

    def require_cmd_or_exit(self, name):
        """Exit with error if *name* is not on PATH."""
        if not shutil.which(name):
            self.log_error(f"Required command not found: {name}")
            sys.exit(127)

    # ------------------------------------------------------------------ #
    # Version logging
    # ------------------------------------------------------------------ #
    # Per-tool version extraction: (tool_name, binary, flags_to_try)
    # Some tools don't support --version; each entry specifies the best strategy.
    VERSION_COMMANDS = [
        # Assemblers
        ("canu",           "canu",             ["--version"]),
        ("nextDenovo",     "nextDenovo",       ["--version"]),
        ("peregrine",      "pg_asm",           ["--version"]),
        ("ipa",            "ipa",              ["--version"]),
        ("flye",           "flye",             ["--version"]),
        ("hifiasm",        "hifiasm",          ["--version"]),
        ("lja",            "lja",              ["--version"]),
        ("mbg",            "MBG",              ["--version"]),
        ("raven",          "raven",            ["--version"]),
        # Analysis
        ("busco",          "busco",            ["--version"]),
        ("quast",          "quast",            ["--version"]),
        ("minimap2",       "minimap2",         ["--version"]),
        ("samtools",       "samtools",         ["--version"]),
        ("seqtk",          "seqtk",            []),  # seqtk prints version on bare call
        ("bwa",            "bwa",              []),   # bwa prints version on bare call
        # Merging / dedup / polishing
        ("merge_wrapper",  "merge_wrapper.py", ["--version"]),
        ("purge_dups",     "purge_dups",       []),   # prints version on bare call
        ("nextPolish2",    "nextPolish2",      ["--version"]),
        ("yak",            "yak",              ["version"]),
        ("racon",          "racon",            ["--version"]),
        ("medaka",         "medaka",           ["--version"]),
        # Merqury / Meryl
        ("merqury",        "merqury.sh",       []),   # no version flag
        ("meryl",          "meryl",            ["--version"]),
        # Runtime
        ("python",         "python3",          ["--version"]),
    ]

    CONDA_PACKAGE_ALIASES = {
        "merge_wrapper": ["quickmerge"],
        "merqury": ["merqury"],
        "nextPolish2": ["nextpolish2"],
        "python": ["python"],
    }

    @staticmethod
    def _version_output_is_error(text):
        """Return True for command output that is clearly not a version."""
        bad_phrases = (
            "failed to parse command line parameters",
            "unrecognized option",
            "unrecognized arguments",
            "invalid option",
            "unknown command",
            "no such option",
            "missing argument",
            "required argument",
            "traceback",
        )
        lowered = text.lower()
        if any(p in lowered for p in bad_phrases):
            return True

        first = text.splitlines()[0].strip().lower() if text.splitlines() else ""
        return first.startswith(("usage:", "error:", "exception:"))

    @staticmethod
    def _version_flags(flags_to_try):
        """Combine preferred tool flags with conservative generic fallbacks."""
        generic = ["--version", "-V", "-v", "version", ""]
        preferred = list(flags_to_try or [])
        out = []
        for flag in preferred + generic:
            if flag not in out:
                out.append(flag)
        return out

    def _conda_package_version(self, label, package_names=None):
        """Return package version from the active conda/mamba environment."""
        names = package_names or self.CONDA_PACKAGE_ALIASES.get(label, [label])
        names = {n.lower().replace("_", "-") for n in names}

        prefixes = []
        for prefix in [os.environ.get("CONDA_PREFIX"), sys.prefix]:
            if prefix and prefix not in prefixes:
                prefixes.append(prefix)

        for prefix in prefixes:
            meta_dir = os.path.join(prefix, "conda-meta")
            if not os.path.isdir(meta_dir):
                continue
            for meta_json in glob.glob(os.path.join(meta_dir, "*.json")):
                try:
                    with open(meta_json) as f:
                        meta = json.load(f)
                except Exception:
                    continue
                pkg_name = (meta.get("name") or "").lower().replace("_", "-")
                if pkg_name in names:
                    version = meta.get("version") or "unknown"
                    return f"{version} (conda package: {meta.get('name')})"

        for manager in ["conda", "mamba", "micromamba"]:
            exe = shutil.which(manager)
            if not exe:
                continue
            try:
                result = subprocess.run(
                    [exe, "list", "--json"],
                    capture_output=True, text=True, timeout=20,
                )
            except Exception:
                continue
            if result.returncode != 0 or not result.stdout:
                continue
            try:
                packages = json.loads(result.stdout)
            except json.JSONDecodeError:
                continue
            for pkg in packages:
                pkg_name = (pkg.get("name") or "").lower().replace("_", "-")
                if pkg_name in names:
                    return f"{pkg.get('version', 'unknown')} (conda package: {pkg.get('name')})"

        return None

    def log_version(self, label, cmd, flags_to_try=None):
        """Extract version string from a tool, handling various output formats."""
        import re as _re

        # Try preferred flags first; also try bare command for tools that print
        # version/help on stderr when invoked without arguments.
        for flag in self._version_flags(flags_to_try):
            try:
                if isinstance(cmd, str) and any(c.isspace() for c in cmd):
                    base_cmd = shlex.split(cmd)
                else:
                    base_cmd = [cmd]
                run_args = base_cmd + ([flag] if flag else [])
                r = subprocess.run(
                    run_args,
                    capture_output=True, text=True, timeout=10,
                )
                text = ((r.stdout or "") + "\n" + (r.stderr or "")).strip()
                if not text:
                    continue
                if self._version_output_is_error(text):
                    continue

                # Search for version-like patterns in output
                # Match: tool_name X.Y.Z, vX.Y.Z, Version X.Y, etc.
                m = _re.search(
                    r'(?:version[:\s]*|v)(\d+\.\d+(?:\.\d+)?(?:[-.]\S*)?)',
                    text, _re.IGNORECASE)
                if m:
                    return m.group(0).strip()

                # Match: bare version number at start of line
                for line in text.split("\n"):
                    line = line.strip()
                    if not line:
                        continue
                    m2 = _re.match(r'^(\d+\.\d+(?:\.\d+)?(?:[-.]\S*)?)', line)
                    if m2:
                        return m2.group(1)

                # If we got output but no version pattern, check for useful first line
                first = text.split("\n")[0].strip()
                if first and "unrecognized" not in first.lower() \
                   and "invalid option" not in first.lower() \
                   and "unknown command" not in first.lower() \
                   and "failed to parse" not in first.lower() \
                   and not first.lower().startswith("usage:"):
                    return first

            except Exception:
                pass

        conda_version = self._conda_package_version(label)
        if conda_version:
            return conda_version
        return "unknown"

    def _tool_version_records(self):
        """Collect machine-readable software version records."""
        records = []
        for label, binary, flags in self.VERSION_COMMANDS:
            bin_path = shutil.which(binary)
            if not bin_path:
                # Try lowercase variant
                bin_path = shutil.which(binary.lower())
            if bin_path:
                ver = self.log_version(label, bin_path, flags)
                records.append({
                    "tool": label,
                    "binary": binary,
                    "path": bin_path,
                    "version": ver,
                    "status": "installed",
                })
            else:
                ver = self._conda_package_version(label)
                records.append({
                    "tool": label,
                    "binary": binary,
                    "path": "",
                    "version": ver or "",
                    "status": "conda_package_only" if ver else "not_installed",
                })
        return records

    def write_versions(self):
        """Write software version summary to version.txt."""
        vf = os.path.join(self.working_dir, "version.txt")
        lines = [
            f"Pipeline: {self.PIPELINE_NAME}",
            f"Run ID:   {self.run_id}",
            f"Date:     {self._ts()}",
            f"Host:     {socket.gethostname()}",
            "",
            "Software versions:",
        ]
        self.version_records = self._tool_version_records()
        for rec in self.version_records:
            if rec["status"] == "installed" or rec["version"]:
                lines.append(f"  {rec['tool']}: {rec['version']}")
            else:
                lines.append(f"  {rec['tool']}: not installed")
        with open(vf, "w") as f:
            f.write("\n".join(lines) + "\n")

    # ------------------------------------------------------------------ #
    # Requirement checking
    # ------------------------------------------------------------------ #
    def check_requirements(self):
        """Check that all required tools are available."""
        from taco.utils import is_assembler_compatible

        core_required = ["python3", "minimap2", "samtools"]
        if self.run_busco:
            core_required.append("busco")
        if not shutil.which("quast.py") and not shutil.which("quast"):
            core_required.append("quast.py/quast")

        assembler_bins = {
            "canu": "canu",
            "nextDenovo": "nextDenovo",
            "peregrine": "pg_asm",
            "ipa": "ipa",
            "flye": "flye",
            "hifiasm": "hifiasm",
            "lja": "lja",
            "mbg": "MBG",
            "raven": "raven",
        }
        platform_bins = []
        for asm, binary in assembler_bins.items():
            if is_assembler_compatible(asm, self.platform):
                if asm == "mbg":
                    if not shutil.which("MBG") and not shutil.which("mbg"):
                        platform_bins.append("MBG/mbg")
                else:
                    platform_bins.append(binary)

        missing = [c for c in core_required if "/" in c or not shutil.which(c)]
        missing.extend(c for c in platform_bins if "/" in c or not shutil.which(c))
        # quast can be quast.py or quast
        missing = sorted(set(missing))
        if missing:
            self.log_warn(
                f"Missing tools for selected platform ({self.platform}): "
                f"{', '.join(missing)}"
            )
            self.log_warn(
                "Some steps may fail. Create and activate the TACO conda "
                "environment first, or install tools manually."
            )

    # ------------------------------------------------------------------ #
    # Reference FASTA resolution
    # ------------------------------------------------------------------ #
    def resolve_reference_fasta(self):
        """Download / decompress --reference if needed; update self.reference_fasta."""
        if not self.reference_fasta:
            return
        src = self.reference_fasta
        # URL download
        if src.startswith("http://") or src.startswith("https://") or src.startswith("ftp://"):
            self.log_info(f"Downloading reference FASTA: {src}")
            local = os.path.join(self.working_dir, "reference_input.fasta")
            if shutil.which("curl"):
                self.run_cmd(f'curl -L "{src}" -o "{local}"')
            elif shutil.which("wget"):
                self.run_cmd(f'wget -O "{local}" "{src}"')
            else:
                self.log_error("Neither curl nor wget found for downloading")
                sys.exit(1)
            src = local
        # Decompress gzip
        if src.endswith(".gz"):
            self.log_info(f"Decompressing: {src}")
            out = src[:-3]
            with gzip.open(src, "rb") as fi, open(out, "wb") as fo:
                shutil.copyfileobj(fi, fo)
            src = out
        if not os.path.isfile(src) or os.path.getsize(src) == 0:
            self.log_error(f"Reference FASTA not found or empty: {src}")
            sys.exit(1)
        self.reference_fasta = src

    # ------------------------------------------------------------------ #
    # NextDenovo configuration
    # ------------------------------------------------------------------ #
    def write_nextdenovo_config(self):
        """Write NextDenovo run config and input fofn."""
        nd_read_map = {"pacbio-hifi": "hifi", "nanopore": "ont", "pacbio": "clr"}
        nd_read_type = nd_read_map.get(self.platform, "hifi")

        fofn = os.path.join(self.working_dir, f"input_{self.project}.fofn")
        with open(fofn, "w") as f:
            f.write(f"{self.fastq}\n")

        lst = os.path.join(self.working_dir, f"reads_{self.project}.lst")
        with open(lst, "w") as f:
            f.write(f"{self.fastq}\n")

        cfg = os.path.join(self.working_dir, f"run_{self.project}.cfg")
        with open(cfg, "w") as f:
            f.write(f"""[General]
job_type = local
job_prefix = nextDenovo
task = assemble
rewrite = yes
deltmp = yes
parallel_jobs = {self.threads}
input_type = raw
read_type = {nd_read_type}
input_fofn = {fofn}
workdir = NextDenovo

[correct_option]
read_cutoff = 1k
genome_size = {self.genomesize}
sort_options = -m 20g -t {self.threads}
minimap2_options_raw = -t {self.threads}
pa_correction = 3
correction_options = -p {self.threads}

[assemble_option]
minimap2_options_cns = -t {self.threads}
nextgraph_options = -a 1
""")

    # ------------------------------------------------------------------ #
    # Benchmarking
    # ------------------------------------------------------------------ #
    def _ensure_benchmark_dir(self):
        """Create benchmark output directory only when benchmarking is enabled."""
        if not self.benchmark:
            return False
        os.makedirs(self.benchmark_dir, exist_ok=True)
        return True

    def _run_elapsed_seconds(self):
        """Return current whole-run wall-clock duration if available."""
        if self.run_start_monotonic is not None:
            return time.monotonic() - self.run_start_monotonic
        if self.run_start_ts:
            return (datetime.now() - self.run_start_ts).total_seconds()
        return None

    @staticmethod
    def _normalize_exit_code(code):
        """Convert SystemExit.code into a numeric process exit code."""
        if code is None:
            return 0
        if isinstance(code, int):
            return code
        try:
            return int(code)
        except (TypeError, ValueError):
            return 1

    @staticmethod
    def _format_dt(dt):
        if not dt:
            return ""
        return dt.isoformat(timespec="seconds")

    def _file_sha256_if_requested(self, path):
        """Compute SHA-256 only when explicitly requested for publication records."""
        if os.environ.get("TACO_BENCHMARK_SHA256", "0") != "1":
            return ""
        if not path or not os.path.isfile(path):
            return ""
        abs_path = os.path.abspath(path)
        cache_key = (abs_path, os.path.getsize(abs_path), os.path.getmtime(abs_path))
        if cache_key in self._sha256_cache:
            return self._sha256_cache[cache_key]
        h = hashlib.sha256()
        with open(abs_path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1024 * 1024), b""):
                h.update(chunk)
        digest = h.hexdigest()
        self._sha256_cache[cache_key] = digest
        return digest

    def _path_metadata(self, path):
        """Return stable metadata for an input/output path."""
        if not path:
            return {
                "path": "",
                "basename": "",
                "exists": False,
                "is_file": False,
                "size_bytes": "",
                "mtime": "",
                "sha256": "",
            }
        abs_path = os.path.abspath(path)
        exists = os.path.exists(abs_path)
        is_file = os.path.isfile(abs_path)
        size = os.path.getsize(abs_path) if is_file else ""
        mtime = ""
        if exists:
            mtime = datetime.fromtimestamp(os.path.getmtime(abs_path)).isoformat(timespec="seconds")
        return {
            "path": abs_path,
            "basename": os.path.basename(abs_path),
            "exists": exists,
            "is_file": is_file,
            "size_bytes": size,
            "mtime": mtime,
            "sha256": self._file_sha256_if_requested(abs_path),
        }

    def _git_metadata(self):
        """Return repository provenance for the TACO checkout."""
        meta = {
            "available": False,
            "repo": self.taco_dir,
            "commit": "",
            "branch": "",
            "dirty": "",
            "status_line_count": "",
        }
        if not shutil.which("git") or not os.path.isdir(os.path.join(self.taco_dir, ".git")):
            return meta

        def git_out(args):
            result = subprocess.run(
                ["git"] + args,
                cwd=self.taco_dir,
                capture_output=True,
                text=True,
                timeout=10,
            )
            if result.returncode != 0:
                return ""
            return (result.stdout or "").strip()

        status = git_out(["status", "--porcelain"])
        meta.update({
            "available": True,
            "commit": git_out(["rev-parse", "HEAD"]),
            "branch": git_out(["rev-parse", "--abbrev-ref", "HEAD"]),
            "dirty": "yes" if status else "no",
            "status_line_count": str(len(status.splitlines())) if status else "0",
        })
        return meta

    def _benchmark_rows(self):
        """Build publication-oriented run metadata rows."""
        fastq_meta = self._path_metadata(self.fastq)
        ref_meta = self._path_metadata(self.reference_fasta)
        git_meta = self._git_metadata()
        elapsed = self._run_elapsed_seconds()
        return [
            ("field", "value"),
            ("pipeline", self.PIPELINE_NAME),
            ("run_id", self.run_id),
            ("run_start", self._format_dt(self.run_start_ts)),
            ("run_end", self._format_dt(self.run_end_ts)),
            ("elapsed_wall_sec", f"{elapsed:.3f}" if elapsed is not None else ""),
            ("date", self._ts()),
            ("command_line", " ".join(shlex.quote(a) for a in sys.argv)),
            ("python_executable", sys.executable),
            ("host", socket.gethostname()),
            ("working_directory", self.working_dir),
            ("taco_package_dir", self.taco_dir),
            ("taco_git_commit", git_meta["commit"]),
            ("taco_git_branch", git_meta["branch"]),
            ("taco_git_dirty", git_meta["dirty"]),
            ("taco_git_status_line_count", git_meta["status_line_count"]),
            ("conda_env", os.environ.get("CONDA_DEFAULT_ENV", "N/A")),
            ("threads", str(self.threads)),
            ("genome_size", self.genomesize),
            ("fastq", fastq_meta["path"]),
            ("fastq_basename", fastq_meta["basename"]),
            ("fastq_exists", str(fastq_meta["exists"])),
            ("fastq_size_bytes", str(fastq_meta["size_bytes"])),
            ("fastq_mtime", fastq_meta["mtime"]),
            ("fastq_sha256", fastq_meta["sha256"] or "not_computed"),
            ("reference", ref_meta["path"]),
            ("reference_basename", ref_meta["basename"]),
            ("reference_exists", str(ref_meta["exists"])),
            ("reference_size_bytes", str(ref_meta["size_bytes"])),
            ("reference_mtime", ref_meta["mtime"]),
            ("reference_sha256", ref_meta["sha256"] or "not_computed"),
            ("platform", self.platform),
            ("taxon", self.taxon),
            ("busco_lineage", self.busco_lineage or ""),
            ("motif", self.motif or "auto"),
            ("telomere_mode", self.telomere_mode),
            ("telo_end_window", str(self.telo_end_window)),
            ("telo_score_window", str(self.telo_score_window)),
            ("telo_kmer_min", str(self.telo_kmer_min)),
            ("telo_kmer_max", str(self.telo_kmer_max)),
            ("assembly_only", str(self.assembly_only)),
            ("benchmark", str(self.benchmark)),
            ("merqury_enabled", str(self.merqury_enable)),
            ("merqury_db", self.merqury_db or "auto"),
            ("merqury_k", str(self.merqury_k)),
            ("no_purge_dups", str(self.no_purge_dups)),
            ("no_polish", str(self.no_polish)),
            ("no_coverage_qc", str(self.no_coverage_qc)),
            ("allow_t2t_replace", str(self.allow_t2t_replace)),
            ("steps", ",".join(str(s) for s in self.steps)),
            ("os", plat.platform()),
            ("cpu_count", str(os.cpu_count() or "N/A")),
            ("benchmark_sha256_requested", os.environ.get("TACO_BENCHMARK_SHA256", "0")),
        ]

    def write_run_metadata(self):
        """Write run metadata TSV."""
        if not self._ensure_benchmark_dir():
            return
        rows = self._benchmark_rows()
        with open(self.bench_run_tsv, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(rows)
        self.write_benchmark_manifest()

    def write_benchmark_manifest(self):
        """Write JSON run manifest for reproducibility and publication audits."""
        if not self._ensure_benchmark_dir():
            return
        row_map = {k: v for k, v in self._benchmark_rows()[1:]}
        manifest = {
            "pipeline": self.PIPELINE_NAME,
            "run_id": self.run_id,
            "run": {
                "start": row_map.get("run_start", ""),
                "end": row_map.get("run_end", ""),
                "elapsed_wall_sec": row_map.get("elapsed_wall_sec", ""),
                "command_line": row_map.get("command_line", ""),
                "working_directory": row_map.get("working_directory", ""),
            },
            "inputs": {
                "fastq": self._path_metadata(self.fastq),
                "reference": self._path_metadata(self.reference_fasta),
                "genome_size": self.genomesize,
                "platform": self.platform,
                "taxon": self.taxon,
            },
            "parameters": {
                "threads": self.threads,
                "steps": list(self.steps),
                "assembly_only": self.assembly_only,
                "busco_lineage": self.busco_lineage,
                "telomere_mode": self.telomere_mode,
                "motif": self.motif or "auto",
                "telo_end_window": self.telo_end_window,
                "telo_score_window": self.telo_score_window,
                "telo_kmer_min": self.telo_kmer_min,
                "telo_kmer_max": self.telo_kmer_max,
                "merqury_enabled": self.merqury_enable,
                "merqury_db": self.merqury_db or "auto",
                "merqury_k": self.merqury_k,
                "no_purge_dups": self.no_purge_dups,
                "no_polish": self.no_polish,
                "no_coverage_qc": self.no_coverage_qc,
                "allow_t2t_replace": self.allow_t2t_replace,
            },
            "environment": {
                "host": row_map.get("host", ""),
                "os": row_map.get("os", ""),
                "cpu_count": row_map.get("cpu_count", ""),
                "python_executable": row_map.get("python_executable", ""),
                "conda_env": row_map.get("conda_env", ""),
            },
            "code": self._git_metadata(),
        }
        with open(self.bench_manifest_json, "w") as f:
            json.dump(manifest, f, indent=2, sort_keys=True)

    def write_benchmark_tool_versions(self):
        """Write machine-readable software versions for methods tables."""
        if not self._ensure_benchmark_dir():
            return
        records = getattr(self, "version_records", None)
        if records is None:
            records = self._tool_version_records()
            self.version_records = records
        with open(self.bench_tools_tsv, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                delimiter="\t",
                fieldnames=["tool", "binary", "path", "version", "status"],
            )
            writer.writeheader()
            writer.writerows(records)

    def write_benchmark_methods_note(self):
        """Write a short methods note that points to the auditable outputs."""
        if not self._ensure_benchmark_dir():
            return
        try:
            from taco.utils import ALL_ASSEMBLERS, is_assembler_compatible
            compatible = [
                asm for asm in ALL_ASSEMBLERS
                if asm != "reference" and is_assembler_compatible(asm, self.platform)
            ]
        except Exception:
            compatible = []

        lines = [
            f"TACO run {self.run_id} ({self.PIPELINE_NAME})",
            "",
            "Recommended citation/reporting note:",
            (
                f"TACO {self.PIPELINE_NAME.replace('TACO-', '')} was run on "
                f"{self.platform} reads with genome size {self.genomesize}, "
                f"{self.threads} threads, taxon preset '{self.taxon}', BUSCO lineage "
                f"'{self.busco_lineage or 'not set'}', and telomere mode "
                f"'{self.telomere_mode}'."
            ),
            (
                "Scientific assembly comparison metrics are reported in "
                "assemblies/assembly_info.csv and, for full refinement runs, "
                "final_results/final_result.csv."
            ),
            (
                "The benchmark_logs directory records run provenance and timing only; "
                "it should be interpreted alongside the biological QC tables rather "
                "than as a replacement for them."
            ),
            "",
            f"Exact command: {' '.join(shlex.quote(a) for a in sys.argv)}",
            f"Compatible assemblers for this platform: {', '.join(compatible) if compatible else 'see logs'}",
            f"Merqury enabled: {self.merqury_enable}; k={self.merqury_k}; db={self.merqury_db or 'auto'}",
        ]
        if self.merqury_enable and self.platform != "pacbio-hifi":
            lines.append(
                "Merqury caution: QV estimates from Nanopore or PacBio CLR reads "
                "can underestimate true assembly quality; compare relative rankings "
                "within the same read set cautiously."
            )
        lines.extend([
            "",
            "Publication-friendly companion files:",
            f"- {os.path.relpath(self.bench_run_tsv, self.working_dir)}",
            f"- {os.path.relpath(self.bench_manifest_json, self.working_dir)}",
            f"- {os.path.relpath(self.bench_tools_tsv, self.working_dir)}",
            f"- {os.path.relpath(self.bench_outputs_tsv, self.working_dir)}",
            f"- {os.path.relpath(self.bench_step_tsv, self.working_dir)}",
        ])
        with open(self.bench_methods_txt, "w") as f:
            f.write("\n".join(lines) + "\n")

    def write_benchmark_output_manifest(self):
        """Write existence/size/mtime records for key result files."""
        if not self._ensure_benchmark_dir():
            return
        paths = [
            "version.txt",
            "assemblies/assembly_info.csv",
            "assemblies/assembly.busco.csv",
            "assemblies/assembly.telo.csv",
            "assemblies/assembly.quast.csv",
            "assemblies/assembly.merqury.csv",
            "assemblies/merged.busco.csv",
            "assemblies/merged.telo.csv",
            "assemblies/merged.quast.csv",
            "assemblies/merged.merqury.csv",
            "assemblies/final.merged.fasta",
            "final_results/final_result.csv",
            "final_results/assembly_info.csv",
            "final_results/assembly_only_result.csv",
            "final_results/final.merged.fasta",
            "final_results/final_assembly.fasta",
            "final_results/final.merged.provenance.gff3",
            "final_results/pool_contig_provenance.tsv",
            "final_results/quickmerge_validation.tsv",
            "final_results/telomere_pool_decisions.tsv",
            "final_results/rescue_trial_summary.tsv",
            "final_results/rescue_rejection_summary.txt",
            "final_results/single_tel.candidates.tsv",
            "final_results/single_tel.replaced.debug.tsv",
            "final_results/selection_decision.txt",
            "final_results/selection_debug.tsv",
            "final_results/coverage_summary.tsv",
            "final_results/weak_regions.tsv",
            "final_results/weak_regions.gff3",
            "final_results/refinement_warning.txt",
            "telomere_pool/pool_contig_provenance.tsv",
            "telomere_pool/telomere_pool_decisions.tsv",
            "telomere_pool/quickmerge_validation.tsv",
            "telomere_pool/protected_telomere_contigs.fasta",
            os.path.relpath(self.bench_run_tsv, self.working_dir),
            os.path.relpath(self.bench_manifest_json, self.working_dir),
            os.path.relpath(self.bench_tools_tsv, self.working_dir),
            os.path.relpath(self.bench_step_tsv, self.working_dir),
            os.path.relpath(self.bench_summary, self.working_dir),
            os.path.relpath(self.bench_methods_txt, self.working_dir),
        ]
        with open(self.bench_outputs_tsv, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["path", "exists", "is_file", "size_bytes", "mtime"])
            for rel_path in paths:
                abs_path = rel_path if os.path.isabs(rel_path) else os.path.join(self.working_dir, rel_path)
                meta = self._path_metadata(abs_path)
                writer.writerow([
                    rel_path,
                    meta["exists"],
                    meta["is_file"],
                    meta["size_bytes"],
                    meta["mtime"],
                ])

    def init_benchmark_step_table(self):
        """Create header for step benchmark TSV."""
        if not self._ensure_benchmark_dir():
            return
        with open(self.bench_step_tsv, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow([
                "run_id", "step", "step_name", "start_time", "end_time",
                "runtime_sec", "status", "exit_code", "log_file",
            ])

    def append_step_benchmark(self, step, start_ts, end_ts, status, log_file, exit_code=0):
        """Append a row to the step benchmark table."""
        if not self.benchmark:
            return
        runtime = (end_ts - start_ts).total_seconds()
        sname = STEP_NAMES.get(step, "Unknown")
        with open(self.bench_step_tsv, "a", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow([
                self.run_id,
                step,
                sname,
                start_ts.strftime('%Y-%m-%d %H:%M:%S'),
                end_ts.strftime('%Y-%m-%d %H:%M:%S'),
                f"{runtime:.3f}",
                status,
                exit_code,
                log_file,
            ])

    def write_benchmark_summary(self):
        """Write a human-readable benchmark summary."""
        if not self.benchmark:
            return
        lines = [
            f"Pipeline: {self.PIPELINE_NAME}",
            f"Run ID: {self.run_id}",
            f"Date: {self._ts()}",
            f"Benchmark table: {self.bench_step_tsv}",
            f"Run metadata: {self.bench_run_tsv}",
            "",
        ]
        total_runtime = 0.0
        step_count = 0
        counts = {}
        if os.path.isfile(self.bench_step_tsv):
            with open(self.bench_step_tsv, newline="") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    if row.get("run_id") != self.run_id:
                        continue
                    step_count += 1
                    status = row.get("status", "unknown")
                    counts[status] = counts.get(status, 0) + 1
                    try:
                        total_runtime += float(row.get("runtime_sec", "0") or 0)
                    except ValueError:
                        pass
        lines.append(f"Steps recorded: {step_count}")
        lines.append(f"Total step runtime: {total_runtime:.3f} sec")
        elapsed = self._run_elapsed_seconds()
        if elapsed is not None:
            lines.append(f"Run wall time: {elapsed:.3f} sec")
        for status in sorted(counts):
            lines.append(f"{status}: {counts[status]} step(s)")
        lines.extend([
            "",
            f"Run manifest: {self.bench_manifest_json}",
            f"Software versions: {self.bench_tools_tsv}",
            f"Output manifest: {self.bench_outputs_tsv}",
            f"Methods note: {self.bench_methods_txt}",
        ])
        with open(self.bench_summary, "w") as f:
            f.write("\n".join(lines) + "\n")
        self.write_run_metadata()
        self.write_benchmark_methods_note()
        self.write_benchmark_output_manifest()

    # ------------------------------------------------------------------ #
    # Main pipeline dispatch
    # ------------------------------------------------------------------ #
    def run(self):
        """Execute the pipeline."""
        self.run_start_ts = datetime.now()
        self.run_start_monotonic = time.monotonic()
        self.run_end_ts = None
        self.run_id = self.run_start_ts.strftime("%Y%m%d_%H%M%S")
        self.log(f"Starting {self.PIPELINE_NAME}")
        self.log_info(f"Sequencing platform: {self.platform}")
        self.log_info(f"Telomere detection mode: {self.telomere_mode}")
        if self.motif:
            self.log_info(f"User telomere motif: {self.motif}")
        else:
            self.log_info("No user motif supplied; using auto-discovery + built-in families")
        if self.assembly_only:
            self.log_info("Assembly-only mode enabled")
        if self.benchmark:
            self.log_info("Benchmark logging enabled: benchmark_logs/")

        self.check_requirements()
        self.resolve_reference_fasta()
        self.write_versions()
        self.write_run_metadata()
        self.write_benchmark_tool_versions()
        self.write_benchmark_methods_note()
        self.init_benchmark_step_table()
        self.write_nextdenovo_config()

        from taco.steps import STEP_FUNCTIONS

        for step_num in self.steps:
            func = STEP_FUNCTIONS.get(step_num)
            if func is None:
                self.log_error(f"Unknown step: {step_num}")
                sys.exit(1)

            sname = STEP_NAMES.get(step_num, f"Step {step_num}")
            log_file = os.path.join(self.logs_dir, f"step_{step_num}.log")

            self.log(f"===== STEP {step_num} START: {sname} =====")
            start = datetime.now()

            # Tee stdout/stderr into the per-step log file
            log_fh = open(log_file, "w")
            old_stdout, old_stderr = sys.stdout, sys.stderr
            sys.stdout = TeeWriter(log_fh, old_stdout)
            sys.stderr = TeeWriter(log_fh, old_stderr)

            try:
                func(self)
                status = "success"
                exit_code = 0
            except SystemExit as e:
                exit_code = self._normalize_exit_code(e.code)
                status = "success" if exit_code == 0 else "failed"
                end = datetime.now()
                self.run_end_ts = end
                sys.stdout, sys.stderr = old_stdout, old_stderr
                log_fh.close()
                self.append_step_benchmark(step_num, start, end, status, log_file, exit_code)
                if exit_code == 0:
                    self.log_info(f"Step {step_num} exited early. See {log_file}")
                else:
                    self.log_error(f"Step {step_num} failed. See {log_file}")
                self.write_benchmark_summary()
                sys.exit(exit_code)
            except Exception as e:
                status = "failed"
                exit_code = 1
                end = datetime.now()
                self.run_end_ts = end
                sys.stdout, sys.stderr = old_stdout, old_stderr
                log_fh.close()
                self.append_step_benchmark(step_num, start, end, status, log_file, exit_code)
                self.log_error(f"Step {step_num} failed with exception: {e}")
                self.write_benchmark_summary()
                raise
            finally:
                # Restore streams in the normal (non-exception) path too
                if sys.stdout is not old_stdout:
                    sys.stdout = old_stdout
                if sys.stderr is not old_stderr:
                    sys.stderr = old_stderr
                if not log_fh.closed:
                    log_fh.close()

            end = datetime.now()
            self.append_step_benchmark(step_num, start, end, status, log_file, exit_code)
            self.log(f"===== STEP {step_num} END: {sname} ({(end - start).total_seconds():.0f}s) =====")

        self.run_end_ts = datetime.now()
        self.write_benchmark_summary()
        self.log(f"{self.PIPELINE_NAME} completed successfully")
