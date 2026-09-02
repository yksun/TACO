"""Pipeline runner for TACO genome assembly."""
import os
import sys
import subprocess
import shutil
import signal
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
from taco.policy import (AssemblyPolicy, DELIVERABLE_MODES,
                         DEFAULT_MODE, MODE_BOTH, SCORING_PROFILES)


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

    PIPELINE_NAME = "TACO-1.5.2"

    def __init__(self, args):
        # Core parameters
        self.genomesize = args.genomesize
        self.threads = args.threads
        self.fastq = os.path.realpath(args.fastq)
        self.motif = args.motif
        self.taxon = getattr(args, 'taxon', 'other')
        self.platform = args.platform
        self.reference_fasta = args.reference
        # Compare-only FASTA — passive QC + contig-vs-final report.  Never
        # used for backbone selection, quickmerge, polish, or refinement.
        self.compare_fasta = getattr(args, 'compare', None)
        # User-supplied final assembly to use in steps 13/14 in place of
        # assemblies/final.merged.fasta.  Used by --final-fa.
        self.final_fa_override = getattr(args, 'final_fa', None)
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

        # Score window: an explicit user value (including 500) takes
        # precedence; otherwise fall back to the taxon-aware default.  The CLI
        # default is None precisely so a deliberate --telo-score-window 500 is
        # not mistaken for "unset".
        user_score_window = getattr(args, 'telo_score_window', None)
        if user_score_window is not None:
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
        self.concordance_mode = getattr(args, 'concordance_mode', 'exclude')
        self.no_contam_screen = getattr(args, 'no_contam_screen', False)
        # --no-contam-screen stays an alias for --purify-mode off, so a v1.3.7
        # command line behaves the way its author intended instead of silently
        # gaining a stage that now modifies the assembly.
        self.purify_mode = ('off' if self.no_contam_screen
                            else getattr(args, 'purify_mode', 'on'))
        self.metagenome = getattr(args, 'metagenome', False)

        # Assembly representation mode (v1.5.0).  Everything that scores an
        # assembly or removes sequence for looking redundant consults this one
        # object, so the two objectives cannot drift apart across call sites.
        # Built after no_purge_dups because the policy honours it.
        self.assembly_mode = getattr(args, 'assembly_mode', DEFAULT_MODE)
        self.policy = AssemblyPolicy(mode=self.assembly_mode, taxon=self.taxon,
                                     user_no_purge_dups=self.no_purge_dups)

        self.chimera_action = getattr(args, 'chimera_action', 'split')
        self.spanning_anchor = getattr(args, 'spanning_anchor', 1000)
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
            # Fallback: use taxon default (should be set in cli.py).
            # Note: "fungal" defaults to fungi_odb10 (the broad fungal lineage);
            # use --busco ascomycota_odb10 (or similar) for sub-lineage scoring.
            from taco.cli import TAXON_BUSCO_LINEAGE
            self.busco_lineage = TAXON_BUSCO_LINEAGE.get(self.taxon, "fungi_odb10")
            self.run_busco = True

        # Where BUSCO should look for / cache lineage datasets.  If neither
        # --busco-download-path nor BUSCO_DOWNLOAD_PATH is set, BUSCO uses
        # ./busco_downloads relative to the working directory (its default).
        cli_dl = getattr(args, 'busco_download_path', None)
        env_dl = os.environ.get("BUSCO_DOWNLOAD_PATH") or None
        self.busco_download_path = cli_dl or env_dl
        # Whether BUSCO is allowed to download missing lineages from the net.
        # Default: yes (this is the behavior most users expect; offline-only
        # users must opt in).  This replaces the earlier
        # STEP12_BUSCO_ALLOW_DOWNLOAD=1 opt-in.
        self.busco_offline_only = bool(getattr(args, 'busco_offline_only', False))

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
    def run_cmd(self, cmd, desc=None, check=True, timeout=None):
        """Run an external command. cmd may be a string (shell) or list.

        *timeout* bounds the wall-clock time in seconds. On expiry the command's
        whole PROCESS GROUP is killed, not just the immediate child: a shell
        pipeline leaves its real worker as a grandchild, so killing only the
        child abandons that worker to run forever. This was observed on a real
        machine, where MUMmer `delta-filter` processes from abandoned `dnadiff`
        runs had each been burning a core for 8 and 11 days.

        A timed-out command is reported as a failure (returncode 124, as
        coreutils `timeout` does) and obeys *check* like any other failure.
        """
        if desc:
            self.log(desc)
        if isinstance(cmd, list):
            display = " ".join(shlex.quote(str(c)) for c in cmd)
        else:
            display = cmd
        self.log(f"$ {display}")
        if timeout:
            self.log_info(f"Time limit: {timeout}s")
        start = time.time()
        # start_new_session puts the child in its own process group so the whole
        # tree can be signalled on timeout.
        proc = subprocess.Popen(
            cmd,
            shell=isinstance(cmd, str),
            start_new_session=True,
        )
        timed_out = False
        try:
            proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            self._kill_process_group(proc, display, timeout)
        elapsed = time.time() - start
        returncode = 124 if timed_out else proc.returncode
        self.log(f"Command exit={returncode} elapsed={elapsed:.1f}s")
        if check and returncode != 0:
            self.log_error(f"Command failed (exit {returncode}): {display}")
            sys.exit(returncode)
        return subprocess.CompletedProcess(cmd, returncode)

    def _kill_process_group(self, proc, display, timeout):
        """SIGTERM then SIGKILL the timed-out command's whole process group."""
        self.log_warn(
            f"Command exceeded its {timeout}s time limit and was abandoned: "
            f"{display}")
        # Read the group id ONCE. Re-reading it after signalling races with the
        # process exiting, and a ProcessLookupError on this path would turn a
        # handled timeout into a crash.
        try:
            pgid = os.getpgid(proc.pid)
        except ProcessLookupError:
            return                      # already gone; nothing to signal
        for sig, label, grace in ((signal.SIGTERM, "SIGTERM", 10),
                                  (signal.SIGKILL, "SIGKILL", 5)):
            try:
                os.killpg(pgid, sig)
            except (ProcessLookupError, PermissionError):
                return
            self.log_info(f"Sent {label} to process group {pgid}")
            try:
                proc.wait(timeout=grace)
                return
            except subprocess.TimeoutExpired:
                continue

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
        ("seqkit",         "seqkit",           ["version"]),  # Step 0 read/base counting
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
        "purge_dups": ["purge_dups", "purge-dups"],
        "raven": ["raven-assembler"],
        "python": ["python"],
    }

    @staticmethod
    def _version_output_is_error(text):
        """Return True for command output that is clearly not a version."""
        bad_phrases = (
            "failed to parse command line parameters",
            "[e::",
            "can not open",
            "cannot open",
            "failed to open",
            "no such file",
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
    @staticmethod
    def _tool_available(tool_spec):
        """Return True if a tool or one of several alternative binaries exists."""
        if isinstance(tool_spec, (tuple, list)):
            return any(shutil.which(t) for t in tool_spec)
        return shutil.which(tool_spec) is not None

    def check_requirements(self):
        """Check tools needed by the requested steps only."""
        from taco.utils import is_assembler_compatible

        step_set = set(self.steps)
        required = {}
        optional = {}

        def add_required(tool_spec, reason):
            required.setdefault(tuple(tool_spec) if isinstance(tool_spec, list) else tool_spec, set()).add(reason)

        def add_optional(tool_spec, reason):
            optional.setdefault(tuple(tool_spec) if isinstance(tool_spec, list) else tool_spec, set()).add(reason)

        assembler_steps = {
            1: ("canu", "canu", True),
            2: ("nextDenovo", "nextDenovo", True),
            3: ("peregrine", "pg_asm", True),
            4: ("ipa", "ipa", True),
            5: ("flye", "flye", True),
            6: ("hifiasm", "hifiasm", True),
            7: ("lja", "lja", True),
            8: ("mbg", ["MBG", "mbg"], False),
            9: ("raven", "raven", True),
        }

        for step, (asm, binary, required_tool) in assembler_steps.items():
            if step not in step_set:
                continue
            if not is_assembler_compatible(asm, self.platform):
                self.log_info(
                    f"Step {step} ({asm}) is incompatible with platform "
                    f"{self.platform}; it will be skipped.")
                continue
            if required_tool:
                add_required(binary, f"Step {step} {asm} assembly")
            else:
                add_optional(binary, f"Step {step} {asm} assembly")

        if 10 in step_set:
            if self.run_busco:
                add_required("busco", "Step 10 BUSCO comparison")
            add_required(["quast.py", "quast"], "Step 10 QUAST comparison")
            if self.merqury_enable:
                add_optional("merqury.sh", "Step 10 Merqury comparison")
                if getattr(self, 'merqury_build_db', False) or not self.merqury_db:
                    add_optional("meryl", "Step 10 Merqury read database build")

        if 11 in step_set:
            add_optional("minimap2", "Step 11 telomere-pool clustering/quickmerge validation")
            add_optional("merge_wrapper.py", "Step 11 optional quickmerge T2T recovery")

        if 12 in step_set:
            add_optional("minimap2", "Step 12 rescue alignment, polishing, coverage QC, purge_dups")
            if self.run_busco:
                add_optional("busco", "Step 12 BUSCO trial validation")
            if self.merqury_enable:
                add_optional("merqury.sh", "Step 12 backbone scoring Merqury refresh")
                if getattr(self, 'merqury_build_db', False) or not self.merqury_db:
                    add_optional("meryl", "Step 12 Merqury read database build")
            if not getattr(self, 'no_coverage_qc', False):
                add_optional("samtools", "Step 12 final assembly coverage QC")
            # Under 'both' the primary deliverable runs purge_dups and the full
            # one does not, so the tool is needed if ANY deliverable wants it.
            if any(p.purge_dups_enabled for p in self.policy.sub_policies()):
                for tool in ["purge_dups", "pbcstat", "calcuts", "split_fa", "get_seqs"]:
                    add_optional(tool, "Step 12 purge_dups cleanup")
            if not getattr(self, 'no_polish', False):
                if self.platform == "pacbio-hifi":
                    add_optional("nextPolish2", "Step 12 HiFi polishing")
                    add_optional("yak", "Step 12 HiFi polishing")
                    add_optional("samtools", "Step 12 HiFi polishing read alignment")
                elif self.platform == "nanopore":
                    add_optional("medaka_consensus", "Step 12 Nanopore polishing")
                    add_optional("racon", "Step 12 Nanopore polishing fallback")
                else:
                    add_optional("racon", "Step 12 PacBio CLR polishing")

        if 13 in step_set:
            if self.run_busco:
                add_required("busco", "Step 13 final BUSCO QC")
            add_required(["quast.py", "quast"], "Step 13 final QUAST QC")
            if self.merqury_enable:
                add_optional("merqury.sh", "Step 13 final Merqury QC")
                if getattr(self, 'merqury_build_db', False) or not self.merqury_db:
                    add_optional("meryl", "Step 13 Merqury read database build")

        missing_required = {
            self._format_tool_spec(tool): sorted(reasons)
            for tool, reasons in required.items()
            if not self._tool_available(tool)
        }
        missing_optional = {
            self._format_tool_spec(tool): sorted(reasons)
            for tool, reasons in optional.items()
            if not self._tool_available(tool)
        }

        if missing_required:
            self.log_warn(
                "Missing required tools for requested steps: " +
                "; ".join(
                    f"{tool} ({', '.join(reasons)})"
                    for tool, reasons in sorted(missing_required.items())
                )
            )
            self.log_warn(
                "Only requested steps were checked. Install these tools or "
                "remove the affected steps; resume runs do not require "
                "assembler binaries unless their assembly steps are selected."
            )

        if missing_optional:
            self.log_warn(
                "Optional tools unavailable for requested steps: " +
                "; ".join(
                    f"{tool} ({', '.join(reasons)})"
                    for tool, reasons in sorted(missing_optional.items())
                )
            )
            self.log_warn(
                "TACO will skip or fall back for optional analyses where "
                "possible; install these tools for the most complete QC."
            )

    @staticmethod
    def _format_tool_spec(tool_spec):
        if isinstance(tool_spec, tuple):
            return "/".join(tool_spec)
        return str(tool_spec)

    @staticmethod
    def _any_path_exists(patterns):
        for pattern in patterns:
            if glob.glob(pattern):
                return True
        return False

    def _copy_resume_input(self, dst, candidates):
        """Restore an active working file from cleanup output if needed."""
        if os.path.exists(dst):
            return False
        for src in candidates:
            if os.path.isfile(src):
                dst_dir = os.path.dirname(dst)
                if dst_dir:
                    os.makedirs(dst_dir, exist_ok=True)
                shutil.copy2(src, dst)
                self.log_warn(
                    f"Resume input restored: copied {src} -> {dst}")
                return True
        return False

    def restore_resume_inputs_for_step(self, step):
        """Restore inputs moved by older cleanup runs before a resumed step.

        After step 14 cleanup, raw assembler dirs (e.g. ``hicanu``) move
        into ``temp/assemblers/``, normalized FASTAs persist in
        ``assemblies/``, stable CSVs are copied into ``final_results/``,
        and the telomere pool products move into ``telomere_pool/``.
        This helper restores whichever inputs the resumed step needs from
        those structured locations so steps 8–14 can be re-run without
        rebuilding the upstream stages.
        """
        # ---- Universal: assemblies/*.result.fasta should be available
        # for any QC / refinement step.  Restore from temp/assemblers/ if
        # cleanup has already moved the raw assembler directories.
        if step in (8, 9, 10, 11, 12, 13, 14):
            from taco.utils import ALL_ASSEMBLERS
            assembler_paths = {
                "canu": ["temp/assemblers/hicanu/canu.contigs.fasta",
                         "hicanu/canu.contigs.fasta"],
                "nextDenovo": ["temp/assemblers/NextDenovo/03.ctg_graph/nd.asm.fasta",
                               "NextDenovo/03.ctg_graph/nd.asm.fasta"],
                "peregrine": ["temp/assemblers/peregrine-2021/asm_ctgs_m_p.fa",
                              "peregrine-2021/asm_ctgs_m_p.fa"],
                "ipa": ["temp/assemblers/ipa/assembly-results/final.p_ctg.fasta",
                        "ipa/assembly-results/final.p_ctg.fasta"],
                "flye": ["temp/assemblers/flye/assembly.fasta",
                         "flye/assembly.fasta"],
                "hifiasm": ["temp/assemblers/hifiasm/hifiasm.fasta",
                            "hifiasm/hifiasm.fasta"],
                "hifiasm_hap": ["temp/assemblers/hifiasm/hifiasm.hap.fasta",
                                "hifiasm/hifiasm.hap.fasta"],
                "ipa_alt": ["temp/assemblers/ipa/assembly-results/ipa.alt.fasta",
                            "ipa/assembly-results/ipa.alt.fasta"],
                "lja": ["temp/assemblers/lja_out/assembly.fasta",
                        "lja_out/assembly.fasta"],
                "mbg": ["temp/assemblers/mbg_out/mbg.fasta",
                        "mbg_out/mbg.fasta"],
                "raven": ["temp/assemblers/raven_out/raven.fasta",
                          "raven_out/raven.fasta"],
            }
            for name, candidates in assembler_paths.items():
                dst = os.path.join("assemblies", f"{name}.result.fasta")
                # Don't overwrite a normalized result file if it's already
                # there; only restore when missing.
                self._copy_resume_input(dst, candidates)

        # ---- assembly_info.csv (used in steps 11/12/13/14)
        if step in (11, 12, 13, 14):
            self._copy_resume_input(
                "assemblies/assembly_info.csv",
                ["final_results/assembly_info.csv",
                 "final_results/assembly_only_result.csv"],
            )

        # ---- per-assembly metric CSVs needed for assembly-only summary
        if step == 14 and self.assembly_only:
            for metric_csv in [
                "assembly.busco.csv",
                "assembly.quast.csv",
                "assembly.telo.csv",
                "assembly.merqury.csv",
            ]:
                self._copy_resume_input(
                    os.path.join("assemblies", metric_csv),
                    [os.path.join("final_results", metric_csv)],
                )

        # ---- merged metric CSVs needed for the full final report
        if step == 14 and not self.assembly_only:
            for metric_csv in [
                "merged.busco.csv",
                "merged.quast.csv",
                "merged.telo.csv",
                "merged.merqury.csv",
            ]:
                self._copy_resume_input(
                    os.path.join("assemblies", metric_csv),
                    [os.path.join("final_results", metric_csv)],
                )

        # ---- selection / decision artefacts (step 12 reads these)
        if step in (11, 12, 13, 14):
            self._copy_resume_input(
                "assemblies/selection_decision.txt",
                ["final_results/selection_decision.txt"],
            )
            self._copy_resume_input(
                "assemblies/selection_debug.tsv",
                ["final_results/selection_debug.tsv"],
            )

        # ---- telomere-pool products (step 12 input)
        if step == 12:
            self._copy_resume_input(
                "pool_contig_provenance.tsv",
                ["telomere_pool/pool_contig_provenance.tsv",
                 "final_results/pool_contig_provenance.tsv"],
            )
            self._copy_resume_input(
                "protected_telomere_contigs.fasta",
                ["telomere_pool/protected_telomere_contigs.fasta"],
            )
            self._copy_resume_input(
                "protected_telomere_mode.txt",
                ["telomere_pool/protected_telomere_mode.txt",
                 "final_results/protected_telomere_mode.txt"],
            )
            for fname in ["t2t.fasta", "single_tel.fasta",
                          "telomere_supported.fasta",
                          "allmerged_telo_sort.fasta",
                          "t2t.list", "single_tel.list",
                          "telomere_supported.list",
                          "telomere_cluster_summary.tsv",
                          "t2t_cluster_summary.tsv",
                          "telomere_support_summary.csv"]:
                self._copy_resume_input(
                    fname,
                    [os.path.join("telomere_pool", fname),
                     os.path.join("final_results", fname)],
                )

        # ---- final.merged.fasta (steps 13/14)
        if step == 13 or (step == 14 and not self.assembly_only):
            override = getattr(self, "final_fa_override", None)
            dst = "assemblies/final.merged.fasta"
            if override:
                # --final-fa is authoritative: copy it in even if a stale
                # final.merged.fasta exists in this working directory.
                if os.path.isfile(override):
                    os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
                    if os.path.realpath(override) != os.path.realpath(dst) \
                            if os.path.exists(dst) else True:
                        shutil.copy2(override, dst)
                        self.log_warn(
                            f"--final-fa: copied {override} -> {dst}")
                else:
                    self.log_error(
                        f"--final-fa points to a missing file: {override}")
                    sys.exit(1)
            else:
                self._copy_resume_input(
                    dst,
                    ["final_results/final.merged.fasta",
                     "final_results/final_assembly.fasta"],
                )

    def warn_missing_step_inputs(self, step):
        """Warn clearly when a selected/resumed step lacks upstream files."""
        fastq_patterns = [self.fastq] if self.fastq else []
        raw_assembler_patterns = [
            "hicanu/canu.contigs.fasta",
            "NextDenovo/03.ctg_graph/nd.asm.fasta",
            "peregrine-2021/asm_ctgs_m_p.fa",
            "ipa/assembly-results/final.p_ctg.fasta",
            "flye/assembly.fasta",
            "hifiasm/hifiasm.fasta",
            "lja_out/assembly.fasta",
            "mbg_out/mbg.fasta",
            "raven_out/raven.fasta",
            "temp/assemblers/hicanu/canu.contigs.fasta",
            "temp/assemblers/NextDenovo/03.ctg_graph/nd.asm.fasta",
            "temp/assemblers/peregrine-2021/asm_ctgs_m_p.fa",
            "temp/assemblers/ipa/assembly-results/final.p_ctg.fasta",
            "temp/assemblers/flye/assembly.fasta",
            "temp/assemblers/hifiasm/hifiasm.fasta",
            "temp/assemblers/lja_out/assembly.fasta",
            "temp/assemblers/mbg_out/mbg.fasta",
            "temp/assemblers/raven_out/raven.fasta",
        ]
        assembler_or_normalized_patterns = (
            ["assemblies/*.result.fasta"] + raw_assembler_patterns
        )
        component_metric_patterns = [
            "assemblies/assembly.busco.csv",
            "assemblies/assembly.quast.csv",
            "assemblies/assembly.telo.csv",
            "assemblies/assembly.merqury.csv",
            "final_results/assembly.busco.csv",
            "final_results/assembly.quast.csv",
            "final_results/assembly.telo.csv",
            "final_results/assembly.merqury.csv",
        ]

        checks = {
            1: [
                ("input FASTQ",
                 fastq_patterns,
                 "the --fastq input / Step 0"),
            ],
            2: [
                ("input FASTQ",
                 fastq_patterns,
                 "the --fastq input / Step 0"),
            ],
            3: [
                ("input FASTQ",
                 fastq_patterns,
                 "the --fastq input / Step 0"),
            ],
            4: [
                ("input FASTQ",
                 fastq_patterns,
                 "the --fastq input / Step 0"),
            ],
            5: [
                ("input FASTQ",
                 fastq_patterns,
                 "the --fastq input / Step 0"),
            ],
            6: [
                ("input FASTQ",
                 fastq_patterns,
                 "the --fastq input / Step 0"),
            ],
            7: [
                ("input FASTQ",
                 fastq_patterns,
                 "the --fastq input / Step 0"),
            ],
            8: [
                ("input FASTQ",
                 fastq_patterns,
                 "the --fastq input / Step 0"),
            ],
            9: [
                ("input FASTQ",
                 fastq_patterns,
                 "the --fastq input / Step 0"),
            ],
            10: [
                ("assembler FASTA files to normalize",
                 assembler_or_normalized_patterns,
                 "Steps 1-9 assembler outputs"),
            ],
            11: [
                ("assembly comparison table and telomere FASTAs",
                 ["assemblies/assembly_info.csv",
                  "assemblies/*.telo.fasta", "assemblies/*.result.fasta"],
                 "Step 10 (normalize + QC)"),
            ],
            12: [
                ("assembly comparison table",
                 ["assemblies/assembly_info.csv"],
                 "Step 10"),
                ("protected telomere contigs",
                 ["protected_telomere_contigs.fasta",
                  "telomere_pool/protected_telomere_contigs.fasta"],
                 "Step 11 (telomere pool)"),
                ("normalized assembler FASTA files",
                 ["assemblies/*.result.fasta"],
                 "Step 10"),
            ],
            13: [
                ("final merged assembly",
                 ["assemblies/final.merged.fasta",
                  "final_results/final.merged.fasta",
                  "final_results/final_assembly.fasta"],
                 "Step 12 (refinement)"),
            ],
            14: [
                ("assembly comparison table",
                 ["assemblies/assembly_info.csv",
                  "final_results/assembly_info.csv",
                  "final_results/assembly_only_result.csv"] + component_metric_patterns,
                 "Step 10 (normalize + QC)"),
            ],
        }
        if step == 14 and not self.assembly_only:
            checks[14].extend([
                ("final merged assembly",
                 ["assemblies/final.merged.fasta",
                  "final_results/final.merged.fasta",
                  "final_results/final_assembly.fasta"],
                 "Step 12 (refinement)"),
                ("final BUSCO metric table",
                 ["assemblies/merged.busco.csv",
                  "final_results/merged.busco.csv"],
                 "Step 13 (final QC)"),
                ("final telomere metric table",
                 ["assemblies/merged.telo.csv",
                  "final_results/merged.telo.csv"],
                 "Step 13 (final QC)"),
                ("final QUAST metric table",
                 ["assemblies/merged.quast.csv",
                  "final_results/merged.quast.csv"],
                 "Step 13 (final QC)"),
            ])
        if step == 13 and self.merqury_enable:
            checks[13].append(
                ("Merqury database for final QC",
                 ["merqury/reads.meryl", "merqury/reads.k*.meryl",
                  "*.meryl"],
                 "Step 10 or Step 13 (auto-built from reads)")
            )
        missing = []
        for desc, patterns, producer in checks.get(step, []):
            if not self._any_path_exists(patterns):
                missing.append((desc, patterns, producer))
        if not missing:
            return

        self.log_warn(
            f"Step {step} is being run without some expected upstream files. "
            "This is common after selecting individual steps or resuming in a "
            "fresh directory, but the step may fail if it cannot rebuild them.")
        for desc, patterns, producer in missing:
            self.log_warn(
                f"  Missing {desc} from {producer}; expected one of: "
                f"{', '.join(patterns)}")
            self.log_warn(
                f"  Finish {producer} first, or provide one of the expected "
                "files before running this step alone.")

    # ------------------------------------------------------------------ #
    # Reference FASTA resolution
    # ------------------------------------------------------------------ #
    @staticmethod
    def _download_local_name(url, stem):
        """Local filename for a downloaded FASTA, preserving a .gz suffix."""
        url_path = url.split("?", 1)[0].rstrip("/")
        return f"{stem}.fasta.gz" if url_path.endswith(".gz") else f"{stem}.fasta"

    def _maybe_gunzip(self, src):
        """Decompress *src* when it is gzip-compressed, detected by magic
        bytes rather than filename.  Remote servers routinely serve gzip
        regardless of the URL suffix, so a name-only check misses them and
        leaves raw gzip bytes to be mis-read as text downstream.  Returns the
        path to the plain-text FASTA (unchanged if not gzip)."""
        try:
            with open(src, "rb") as fh:
                is_gz = fh.read(2) == b"\x1f\x8b"
        except OSError:
            return src
        if not is_gz:
            return src
        out = src[:-3] if src.endswith(".gz") else src + ".decompressed.fasta"
        self.log_info(f"Decompressing gzip FASTA: {src} -> {out}")
        with gzip.open(src, "rb") as fi, open(out, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        return out

    def _download_file(self, url, local):
        """Download *url* to *local* using curl/wget via argv (no shell)."""
        if shutil.which("curl"):
            self.run_cmd(["curl", "-L", url, "-o", local])
        elif shutil.which("wget"):
            self.run_cmd(["wget", "-O", local, url])
        else:
            self.log_error("Neither curl nor wget found for downloading")
            sys.exit(1)

    def resolve_reference_fasta(self):
        """Download / decompress --reference if needed; update self.reference_fasta."""
        if not self.reference_fasta:
            return
        src = self.reference_fasta
        # URL download (argv form so shell metacharacters in the URL are safe)
        if src.startswith("http://") or src.startswith("https://") or src.startswith("ftp://"):
            self.log_info(f"Downloading reference FASTA: {src}")
            local = os.path.join(self.working_dir,
                                 self._download_local_name(src, "reference_input"))
            self._download_file(src, local)
            src = local
        # Decompress by content (gzip magic bytes), not filename
        src = self._maybe_gunzip(src)
        if not os.path.isfile(src) or os.path.getsize(src) == 0:
            self.log_error(f"Reference FASTA not found or empty: {src}")
            sys.exit(1)
        self.reference_fasta = src

    def resolve_compare_fasta(self):
        """Download / decompress --compare if needed; update self.compare_fasta.

        Mirrors resolve_reference_fasta but for --compare.  Compare-only
        assemblies receive QC and a contig-to-contig comparison report
        against TACO's final.merged.fasta but never participate in
        backbone selection or refinement.
        """
        if not self.compare_fasta:
            return
        src = self.compare_fasta
        if src.startswith("http://") or src.startswith("https://") or src.startswith("ftp://"):
            self.log_info(f"Downloading compare FASTA: {src}")
            local = os.path.join(self.working_dir,
                                 self._download_local_name(src, "compare_input"))
            self._download_file(src, local)
            src = local
        src = self._maybe_gunzip(src)
        if not os.path.isfile(src) or os.path.getsize(src) == 0:
            self.log_error(f"Compare FASTA not found or empty: {src}")
            sys.exit(1)
        self.compare_fasta = src

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
            ("assembly_mode", self.assembly_mode),
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
                "assembly_mode": self.assembly_mode,
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
                f"{self.threads} threads, taxon preset '{self.taxon}', assembly "
                f"representation mode '{self.assembly_mode}', BUSCO lineage "
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
        default_full_steps = list(range(0, 15))
        if self.assembly_only:
            self.log_info("Assembly-only mode enabled")
        elif self.steps != default_full_steps:
            self.log_info(
                "Selected-step/resume mode enabled: " +
                ",".join(str(s) for s in self.steps))
        if self.benchmark:
            self.log_info("Benchmark logging enabled: benchmark_logs/")

        self.check_requirements()
        self.resolve_reference_fasta()
        self.resolve_compare_fasta()
        self._absolutise_input_paths()
        self.write_versions()
        self.write_run_metadata()
        self.write_benchmark_tool_versions()
        self.write_benchmark_methods_note()
        self.init_benchmark_step_table()
        self.write_nextdenovo_config()

        from taco.steps import STEP_FUNCTIONS

        if self.policy.is_both:
            self._run_both_modes(STEP_FUNCTIONS)
        else:
            self._run_step_list(self.steps, STEP_FUNCTIONS)

        self.run_end_ts = datetime.now()
        self.write_benchmark_summary()
        self.log(f"{self.PIPELINE_NAME} completed successfully")

    # ── step execution ──────────────────────────────────────────────────────

    #: Steps whose products describe the READS and the assemblers, so they are
    #: identical whichever representation is being delivered and are computed
    #: once under ``--assembly-mode both``.
    SHARED_STEPS = frozenset(range(0, 11))
    #: Steps that build a deliverable from a chosen backbone.  These depend on
    #: the objective and are run once per representation.
    PER_MODE_STEPS = frozenset(range(11, 15))

    #: Shared step 0-10 products each per-mode workspace needs to see.  Large
    #: FASTAs are symlinked; small tables are copied because the per-mode steps
    #: rewrite some of them.
    SHARED_LINK_GLOBS = ("assemblies/*.result.fasta",
                         "assemblies/*.telo.fasta")
    SHARED_COPY_GLOBS = ("assemblies/*.telo.list",
                         "assemblies/*.telomere_end_scores.tsv",
                         "assemblies/*.telo_metrics.tsv",
                         "assemblies/assembly_info.csv",
                         "assemblies/assembly.*.csv",
                         "assemblies/t2t_concordance.tsv")
    #: Directories of shared per-assembly results.  Their CONTENTS are linked
    #: rather than the directory itself, because the per-mode steps add their
    #: own entries (busco/final, merqury/final) alongside them.
    SHARED_LINK_DIRS = ("busco", "merqury", "assemblies/concordance")
    #: Label steps 13 and 14 use for the assembly a run DELIVERS, as opposed to
    #: the assemblies it evaluates.  Products under this label belong to one
    #: representation and must never be shared between them: linking them would
    #: give both deliverables the same final QC results, and the first attempt to
    #: refresh one would try to remove a symlink.
    DELIVERABLE_LABEL = "final"

    @classmethod
    def _is_deliverable_product(cls, name):
        base = os.path.basename(name)
        return (base == cls.DELIVERABLE_LABEL
                or base.startswith(cls.DELIVERABLE_LABEL + ".")
                or base.startswith(cls.DELIVERABLE_LABEL + "_"))

    def _prepare_mode_workspace(self, mode):
        """Create ``mode_<mode>/`` sharing every step 0-10 product.

        The per-mode steps address their inputs and outputs by paths relative to
        the working directory, so giving each representation its own working
        directory keeps two deliverables from overwriting one another without
        having to thread an output prefix through all of steps 11-14.
        """
        ws = os.path.abspath(f"mode_{mode}")
        os.makedirs(os.path.join(ws, "assemblies"), exist_ok=True)

        def _link(src_abs, dst):
            if os.path.lexists(dst):
                return
            os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
            try:
                os.symlink(src_abs, dst)
            except OSError as e:
                self.log_warn(f"could not link {src_abs} -> {dst}: {e}")

        n_link = n_copy = n_skip = 0
        for pattern in self.SHARED_LINK_GLOBS:
            for src in glob.glob(pattern):
                if self._is_deliverable_product(src):
                    n_skip += 1
                    continue
                _link(os.path.abspath(src), os.path.join(ws, src))
                n_link += 1
        for pattern in self.SHARED_COPY_GLOBS:
            for src in glob.glob(pattern):
                if self._is_deliverable_product(src):
                    n_skip += 1
                    continue
                dst = os.path.join(ws, src)
                os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
                shutil.copy2(src, dst)
                n_copy += 1
        for d in self.SHARED_LINK_DIRS:
            if not os.path.isdir(d):
                continue
            os.makedirs(os.path.join(ws, d), exist_ok=True)
            for entry in os.listdir(d):
                if self._is_deliverable_product(entry):
                    n_skip += 1
                    continue
                _link(os.path.abspath(os.path.join(d, entry)),
                      os.path.join(ws, d, entry))
                n_link += 1

        self.log_info(f"mode_{mode}/: {n_link} shared products linked, "
                      f"{n_copy} copied, {n_skip} deliverable product(s) left "
                      f"for this representation to build")
        return ws

    #: Runner attributes holding paths supplied from outside the working
    #: directory.  ``--assembly-mode both`` runs the per-mode steps with the CWD
    #: inside ``mode_<name>/``, so a relative value here would silently resolve
    #: against the wrong directory.  Absolutising unconditionally also makes a
    #: single-mode run independent of any later chdir.
    EXTERNAL_PATH_ATTRS = ("reference_fasta", "compare_fasta", "final_fa_override",
                           "busco_download_path", "merqury_db")

    def _absolutise_input_paths(self):
        for attr in self.EXTERNAL_PATH_ATTRS:
            val = getattr(self, attr, None)
            if val and not os.path.isabs(val):
                setattr(self, attr, os.path.abspath(val))
                self.log_info(f"{attr}: resolved to {getattr(self, attr)}")

    def _run_both_modes(self, step_functions):
        """Deliver a primary and a full genome from one set of assemblies.

        Steps 0-10 describe the reads and the assemblers, so they run once in the
        run root.  Steps 11-14 depend on the objective and run once per
        representation inside ``mode_primary/`` and ``mode_full/``.
        """
        shared = [s for s in self.steps if s in self.SHARED_STEPS]
        per_mode = [s for s in self.steps if s in self.PER_MODE_STEPS]

        self.log(f"Assembly representation mode: {MODE_BOTH} — delivering "
                 f"{' and '.join(DELIVERABLE_MODES)}")
        if shared:
            self.log_info(f"Shared steps (run once): "
                          f"{','.join(str(s) for s in shared)}")
        if per_mode:
            self.log_info(f"Per-representation steps (run once per mode): "
                          f"{','.join(str(s) for s in per_mode)}")

        if shared:
            self._run_step_list(shared, step_functions)

        if not per_mode:
            self.log_info("No per-representation steps requested; "
                          "nothing further to deliver.")
            return

        root = os.getcwd()
        root_logs, root_mode = self.logs_dir, self.assembly_mode
        root_policy = self.policy
        self._mark_root_results_in_flight(
            [p.mode for p in root_policy.sub_policies()])
        delivered, failed = [], []
        try:
            for sub in root_policy.sub_policies():
                mode = sub.mode
                ws = self._prepare_mode_workspace(mode)
                self.log(f"===== DELIVERABLE {mode}: "
                         f"{SCORING_PROFILES[mode]['description']} =====")
                os.chdir(ws)
                # Point the runner at this representation for the duration.
                self.assembly_mode, self.policy = mode, sub
                self.logs_dir = os.path.join(ws, os.path.basename(root_logs))
                os.makedirs(self.logs_dir, exist_ok=True)
                try:
                    self._run_step_list(per_mode, step_functions)
                    delivered.append((mode, ws))
                # The representations are independent products, so one failing
                # must not cancel the other -- these runs take days, and a
                # failure in the first previously meant the second never began.
                # SystemExit is caught explicitly because _run_step_list exits
                # rather than raising on a failed step, and SystemExit does not
                # derive from Exception.
                except (Exception, SystemExit) as e:
                    failed.append((mode, e))
                    self.log_error(
                        f"Deliverable {mode} failed ({e}); continuing with the "
                        f"remaining representation(s). Its partial output is in "
                        f"{ws}.")
                finally:
                    os.chdir(root)
                    self.assembly_mode, self.policy = root_mode, root_policy
                    self.logs_dir = root_logs
        finally:
            os.chdir(root)
            self.assembly_mode, self.policy = root_mode, root_policy
            self.logs_dir = root_logs

        self._write_both_modes_report(delivered)

        if failed:
            for mode, e in failed:
                self.log_error(f"Deliverable {mode} did not complete: {e}")
            if delivered:
                self.log_warn(
                    f"{len(delivered)} of {len(delivered) + len(failed)} "
                    f"representations completed: "
                    f"{', '.join(m for m, _ in delivered)}. Re-run the failed "
                    f"one with --assembly-mode <mode> once the cause is fixed; "
                    f"completed work is preserved in its mode_* directory.")
            self.run_end_ts = datetime.now()
            self.write_benchmark_summary()
            sys.exit(1)

    #: Written at the run root while the per-mode phase is in flight, so a
    #: pre-existing final_results/ from an earlier run cannot be mistaken for
    #: this one's output.  Removed once the combined report is written.
    IN_FLIGHT_MARKER = "RESULTS_NOT_READY.txt"

    def _mark_root_results_in_flight(self, modes):
        """Say plainly that the root results are not this run's output yet.

        Under ``both`` the per-deliverable reports are written inside
        ``mode_<name>/final_results/``, so anything already sitting in the root
        ``final_results/`` belongs to a previous run.  Without this marker that
        stale file reads as current for as long as the run takes.
        """
        os.makedirs("final_results", exist_ok=True)
        stale = sorted(f for f in os.listdir("final_results")
                       if f != self.IN_FLIGHT_MARKER)
        path = os.path.join("final_results", self.IN_FLIGHT_MARKER)
        with open(path, "w") as f:
            f.write(
                "THIS RUN IS STILL IN PROGRESS — the files in this directory are\n"
                "NOT its output.\n\n"
                f"Run          : {self.run_id} ({self.PIPELINE_NAME})\n"
                f"Mode         : {MODE_BOTH} ({' + '.join(modes)})\n\n"
                "Under --assembly-mode both, each representation is reported inside\n"
                "its own directory while the run proceeds:\n"
                + "".join(f"    mode_{m}/final_results/\n" for m in modes) +
                "\nWhen the run finishes, this directory receives a combined\n"
                "final_result.csv carrying every assembler alongside a\n"
                + "".join(f"    merged_{m}\n" for m in modes) +
                "column, plus assembly_modes_comparison.tsv, and this file is\n"
                "removed.\n")
            if stale:
                f.write(
                    "\nFiles present here NOW are left over from an earlier run and\n"
                    "are superseded, not current:\n"
                    + "".join(f"    {n}\n" for n in stale))
        self.log_warn(
            f"Root final_results/ holds {len(stale)} file(s) from an earlier run; "
            f"see final_results/{self.IN_FLIGHT_MARKER}. This run reports into "
            f"{', '.join('mode_' + m + '/' for m in modes)} until it completes.")

    @staticmethod
    def _read_merged_column(csv_path):
        """Metric -> value for the 'merged' column of a final_result.csv."""
        vals = {}
        if not os.path.isfile(csv_path):
            return vals
        with open(csv_path, newline="") as f:
            r = csv.reader(f)
            hdr = next(r, None) or []
            lower = [h.strip().lower() for h in hdr]
            try:
                mi = lower.index("merged")
            except ValueError:
                mi = len(hdr) - 1
            for row in r:
                if row and row[0].strip():
                    vals[row[0].strip()] = row[mi] if len(row) > mi else ""
        return vals

    def _write_both_modes_report(self, delivered):
        """Combine the deliverables into one table at the run root.

        Produces ``final_results/final_result.csv`` carrying every assembler
        column from step 10 beside one ``merged_<mode>`` column per
        representation, so the assemblers and both deliverables can be read off a
        single table, plus a compact side-by-side of the headline metrics.
        """
        if not delivered:
            return
        os.makedirs("final_results", exist_ok=True)
        per_mode = {m: self._read_merged_column(
                        os.path.join(ws, "final_results", "final_result.csv"))
                    for m, ws in delivered}
        modes = [m for m, _ in delivered]

        # Metric order: step 10's table first, then anything only a deliverable
        # reports (selection metadata, pool counts) in first-seen order.
        rows, assemblers = [], []
        info = os.path.join("assemblies", "assembly_info.csv")
        if os.path.isfile(info):
            with open(info, newline="") as f:
                r = csv.reader(f)
                hdr = next(r, None) or ["Metric"]
                assemblers = [h.strip() for h in hdr[1:]]
                rows = [row for row in r if row and row[0].strip()]
        seen = {row[0].strip() for row in rows}
        extra = []
        for m in modes:
            for metric in per_mode[m]:
                if metric not in seen:
                    seen.add(metric)
                    extra.append(metric)

        out = os.path.join("final_results", "final_result.csv")
        with open(out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Metric"] + assemblers
                       + [f"merged_{m}" for m in modes])
            for row in rows:
                metric = row[0].strip()
                body = list(row[1:]) + [""] * (len(assemblers) - len(row[1:]))
                w.writerow([metric] + body
                           + [per_mode[m].get(metric, "") for m in modes])
            for metric in extra:
                w.writerow([metric] + [""] * len(assemblers)
                           + [per_mode[m].get(metric, "") for m in modes])
        self.log(f"Wrote {out} "
                 f"({len(rows) + len(extra)} metrics x "
                 f"{len(assemblers)} assemblers + {len(modes)} deliverables)")

        # Compact side-by-side of the headline metrics.
        fields = ["Assembly representation mode", "Selected assembler",
                  "Selection score", "BUSCO C (%)", "BUSCO S (%)", "BUSCO D (%)",
                  "Merqury QV", "Merqury completeness (%)", "Total length",
                  "# contigs", "N50", "Telomere strict T2T contigs",
                  "Telomere single-end strong contigs"]
        cmp_path = os.path.join("final_results", "assembly_modes_comparison.tsv")
        with open(cmp_path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["Metric"] + modes)
            for field in fields:
                vals = [(m if field == "Assembly representation mode"
                         else per_mode[m].get(field, "")) for m in modes]
                w.writerow([field] + vals)
        self.log(f"Wrote {cmp_path}")

        for mode, ws in delivered:
            fa = os.path.join(ws, "final_results", "final.merged.fasta")
            alt = os.path.join(ws, "assemblies", "final.merged.fasta")
            self.log(f"  {mode:<8} genome: "
                     f"{fa if os.path.isfile(fa) else alt}")
        marker = os.path.join("final_results", self.IN_FLIGHT_MARKER)
        if os.path.isfile(marker):
            os.remove(marker)
        self.log_info(
            "The two genomes are different REPRESENTATIONS of one sample, not "
            "competing attempts at the same thing. Neither is the 'correct' one; "
            "pick by what the downstream analysis needs. The full "
            "representation retains alternate sequence and is NOT phased.")

    def _run_step_list(self, steps, step_functions):
        for step_num in steps:
            func = step_functions.get(step_num)
            if func is None:
                self.log_error(f"Unknown step: {step_num}")
                sys.exit(1)

            sname = STEP_NAMES.get(step_num, f"Step {step_num}")
            log_file = os.path.join(self.logs_dir, f"step_{step_num}.log")

            # Tee stdout/stderr into the per-step log file
            log_fh = open(log_file, "w")
            old_stdout, old_stderr = sys.stdout, sys.stderr
            sys.stdout = TeeWriter(log_fh, old_stdout)
            sys.stderr = TeeWriter(log_fh, old_stderr)
            start = datetime.now()

            try:
                self.log(f"===== STEP {step_num} START: {sname} =====")
                self.restore_resume_inputs_for_step(step_num)
                self.warn_missing_step_inputs(step_num)
                func(self)
                status = "success"
                exit_code = 0
            except SystemExit as e:
                exit_code = self._normalize_exit_code(e.code)
                status = "success" if exit_code == 0 else "failed"
                end = datetime.now()
                self.run_end_ts = end
                self.log(
                    f"===== STEP {step_num} END: {sname} "
                    f"({(end - start).total_seconds():.0f}s, "
                    f"status={status}, exit={exit_code}) =====")
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
                self.log_error(f"Step {step_num} exception: {e}")
                self.log(
                    f"===== STEP {step_num} END: {sname} "
                    f"({(end - start).total_seconds():.0f}s, "
                    f"status={status}, exit={exit_code}) =====")
                sys.stdout, sys.stderr = old_stdout, old_stderr
                log_fh.close()
                self.append_step_benchmark(step_num, start, end, status, log_file, exit_code)
                self.log_error(f"Step {step_num} failed with exception: {e}")
                self.write_benchmark_summary()
                raise
            else:
                end = datetime.now()
                self.log(
                    f"===== STEP {step_num} END: {sname} "
                    f"({(end - start).total_seconds():.0f}s, "
                    f"status={status}, exit={exit_code}) =====")
                sys.stdout, sys.stderr = old_stdout, old_stderr
                log_fh.close()
                self.append_step_benchmark(step_num, start, end, status, log_file, exit_code)
            finally:
                # Restore streams if an exception path left them active.
                if sys.stdout is not old_stdout:
                    sys.stdout = old_stdout
                if sys.stderr is not old_stderr:
                    sys.stderr = old_stderr
                if not log_fh.closed:
                    log_fh.close()
