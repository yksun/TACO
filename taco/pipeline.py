"""Pipeline runner for TACO genome assembly."""
import os
import sys
import subprocess
import shutil
import gzip
import socket
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

    PIPELINE_NAME = "TACO-1.2.0"

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

        # Backbone selection
        self.auto_mode = args.auto_mode
        self.assembler = None
        self.choose_flag = False
        if args.choose is not None:
            self.choose_flag = True
            if args.choose != "__prompt__":
                self.assembler = args.choose

        # BUSCO
        self.busco_lineage = args.busco if args.busco else "ascomycota_odb10"
        self.run_busco = args.busco is not None

        # Merqury — auto-detect if installed and a .meryl database is found.
        # Explicitly disabled with --no-merqury; explicitly enabled with
        # --merqury or --merqury-db.  Otherwise, auto-enable when merqury.sh
        # is on PATH and a .meryl directory is discoverable.
        self.merqury_db = getattr(args, 'merqury_db', None)
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

        # Derived paths
        fastq_name = os.path.basename(self.fastq)
        self.project = fastq_name.replace('.fastq.gz', '').replace('.fq.gz', '').replace('.fastq', '').replace('.fq', '')
        self.taco_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.working_dir = os.getcwd()

        # Logging and benchmarking
        self.logs_dir = os.path.join(self.working_dir, "logs")
        self.benchmark_dir = os.path.join(self.working_dir, "benchmark_logs")
        os.makedirs(self.logs_dir, exist_ok=True)
        os.makedirs(self.benchmark_dir, exist_ok=True)
        self.run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.bench_step_tsv = os.path.join(self.benchmark_dir, "step_benchmark.tsv")
        self.bench_run_tsv = os.path.join(self.benchmark_dir, "run_metadata.tsv")
        self.bench_summary = os.path.join(self.benchmark_dir, "run_summary.txt")

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
    def log_version(self, label, cmd):
        """Try to capture the version of *cmd* and return the first line."""
        for flag in ["--version", "-V", "-v", "version"]:
            try:
                r = subprocess.run(
                    f"{cmd} {flag}",
                    shell=True,
                    capture_output=True,
                    text=True,
                    timeout=10,
                )
                out = (r.stdout or r.stderr or "").strip()
                if out:
                    return out.split("\n")[0]
            except Exception:
                pass
        return "unknown"

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
        tools = [
            "canu", "nextDenovo", "pg_asm", "ipa", "flye", "hifiasm",
            "seqtk", "busco", "minimap2", "bwa", "samtools",
            "merge_wrapper.py", "python3",
        ]
        for t in tools:
            if shutil.which(t) or shutil.which(t.replace(".py", "")):
                ver = self.log_version(t, t)
                lines.append(f"{t}: {ver}")
            else:
                lines.append(f"{t}: NOT FOUND")
        # quast special case
        qcmd = shutil.which("quast.py") or shutil.which("quast")
        if qcmd:
            ver = self.log_version("quast", qcmd)
            lines.append(f"quast: {ver}")
        else:
            lines.append("quast: NOT FOUND")
        with open(vf, "w") as f:
            f.write("\n".join(lines) + "\n")

    # ------------------------------------------------------------------ #
    # Requirement checking
    # ------------------------------------------------------------------ #
    def check_requirements(self):
        """Check that all required tools are available."""
        required = [
            "python3", "canu", "nextDenovo", "pg_asm", "ipa",
            "flye", "hifiasm", "seqtk", "busco", "minimap2",
            "bwa", "samtools",
        ]
        missing = [c for c in required if not shutil.which(c)]
        # quast can be quast.py or quast
        if not shutil.which("quast.py") and not shutil.which("quast"):
            missing.append("quast.py/quast")
        if missing:
            self.log_warn(
                f"Missing tools in the active environment: {', '.join(missing)}"
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
    def write_run_metadata(self):
        """Write run metadata TSV."""
        rows = [
            ("field", "value"),
            ("pipeline", self.PIPELINE_NAME),
            ("run_id", self.run_id),
            ("date", self._ts()),
            ("host", socket.gethostname()),
            ("working_directory", self.working_dir),
            ("conda_env", os.environ.get("CONDA_DEFAULT_ENV", "N/A")),
            ("threads", str(self.threads)),
            ("genome_size", self.genomesize),
            ("fastq", self.fastq),
            ("platform", self.platform),
            ("motif", self.motif or "auto"),
            ("telomere_mode", self.telomere_mode),
            ("steps", ",".join(str(s) for s in self.steps)),
            ("os", plat.platform()),
            ("cpu_count", str(os.cpu_count() or "N/A")),
        ]
        with open(self.bench_run_tsv, "w") as f:
            for k, v in rows:
                f.write(f"{k}\t{v}\n")

    def init_benchmark_step_table(self):
        """Create header for step benchmark TSV."""
        if not os.path.exists(self.bench_step_tsv):
            with open(self.bench_step_tsv, "w") as f:
                f.write("run_id\tstep\tstep_name\tstart_time\tend_time\truntime_sec\tstatus\tlog_file\n")

    def append_step_benchmark(self, step, start_ts, end_ts, status, log_file):
        """Append a row to the step benchmark table."""
        runtime = (end_ts - start_ts).total_seconds()
        sname = STEP_NAMES.get(step, "Unknown")
        with open(self.bench_step_tsv, "a") as f:
            f.write(
                f"{self.run_id}\t{step}\t{sname}\t"
                f"{start_ts.strftime('%Y-%m-%d %H:%M:%S')}\t"
                f"{end_ts.strftime('%Y-%m-%d %H:%M:%S')}\t"
                f"{runtime:.0f}\t{status}\t{log_file}\n"
            )

    def write_benchmark_summary(self):
        """Write a human-readable benchmark summary."""
        lines = [
            f"Pipeline: {self.PIPELINE_NAME}",
            f"Run ID: {self.run_id}",
            f"Date: {self._ts()}",
            f"Benchmark table: {self.bench_step_tsv}",
            f"Run metadata: {self.bench_run_tsv}",
            "",
        ]
        if os.path.isfile(self.bench_step_tsv):
            counts = {}
            with open(self.bench_step_tsv) as f:
                next(f, None)  # skip header
                for row in f:
                    parts = row.strip().split("\t")
                    if len(parts) >= 7:
                        st = parts[6]
                        counts[st] = counts.get(st, 0) + 1
            for k, v in counts.items():
                lines.append(f"{k}: {v} step(s)")
        with open(self.bench_summary, "w") as f:
            f.write("\n".join(lines) + "\n")

    # ------------------------------------------------------------------ #
    # Main pipeline dispatch
    # ------------------------------------------------------------------ #
    def run(self):
        """Execute the pipeline."""
        self.log(f"Starting {self.PIPELINE_NAME}")
        self.log_info(f"Sequencing platform: {self.platform}")
        self.log_info(f"Telomere detection mode: {self.telomere_mode}")
        if self.motif:
            self.log_info(f"User telomere motif: {self.motif}")
        else:
            self.log_info("No user motif supplied; using auto-discovery + built-in families")
        if self.assembly_only:
            self.log_info("Assembly-only mode enabled")

        self.check_requirements()
        self.write_versions()
        self.write_run_metadata()
        self.init_benchmark_step_table()
        self.resolve_reference_fasta()
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
                status = "0"
            except SystemExit as e:
                status = str(e.code) if e.code else "1"
                end = datetime.now()
                sys.stdout, sys.stderr = old_stdout, old_stderr
                log_fh.close()
                self.append_step_benchmark(step_num, start, end, status, log_file)
                self.log_error(f"Step {step_num} failed. See {log_file}")
                self.write_benchmark_summary()
                sys.exit(int(status))
            except Exception as e:
                status = "1"
                end = datetime.now()
                sys.stdout, sys.stderr = old_stdout, old_stderr
                log_fh.close()
                self.append_step_benchmark(step_num, start, end, status, log_file)
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
            self.append_step_benchmark(step_num, start, end, status, log_file)
            self.log(f"===== STEP {step_num} END: {sname} ({(end - start).total_seconds():.0f}s) =====")

        self.write_benchmark_summary()
        self.log(f"{self.PIPELINE_NAME} completed successfully")
