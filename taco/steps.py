"""
Step functions for TACO genome assembly pipeline.

Each step is implemented as a function that takes a PipelineRunner instance
and performs the specified pipeline operation.
"""
import os
import sys
import glob
import csv
import re
import shutil
import subprocess
import json
import math
import tempfile
from pathlib import Path
from collections import defaultdict

from taco.telomere_detect import detect_telomeres, write_detection_outputs
from taco.clustering import cluster_and_select


def _parse_genome_size(size_str):
    """Parse genome size string like '12m', '500k', '1.5g' to integer bp.

    Returns 0 if the string cannot be parsed.
    """
    if not size_str:
        return 0
    size_str = size_str.strip().lower()
    multipliers = {'k': 1_000, 'm': 1_000_000, 'g': 1_000_000_000}
    try:
        if size_str[-1] in multipliers:
            return int(float(size_str[:-1]) * multipliers[size_str[-1]])
        return int(float(size_str))
    except (ValueError, IndexError):
        return 0


def rename_and_sort_fasta(runner, infa, outfa, prefix):
    """
    Read FASTA, sort contigs by length descending, rename as prefix_1, prefix_2, etc.
    """
    if not os.path.isfile(infa) or os.path.getsize(infa) == 0:
        runner.log_warn(f"rename_and_sort_fasta: input '{infa}' missing or empty; skipping.")
        return

    def read_fasta(fp):
        name = None
        seq = []
        for line in fp:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if name is not None:
            yield name, ''.join(seq)

    def wrap(s, w=60):
        return '\n'.join(s[i:i+w] for i in range(0, len(s), w)) if s else ''

    with open(infa, 'r') as f:
        recs = list(read_fasta(f))

    recs.sort(key=lambda x: len(x[1]), reverse=True)

    with open(outfa, 'w') as o:
        for i, (name, seq) in enumerate(recs, start=1):
            o.write(f">{prefix}_{i}\n")
            o.write(wrap(seq) + ("\n" if seq and not seq.endswith("\n") else ""))


def _extract_by_list(listfile, infa, outfa):
    """Extract sequences from FASTA by list of IDs."""
    if not os.path.isfile(listfile) or os.path.getsize(listfile) == 0:
        with open(outfa, 'w') as f:
            pass
        return

    with open(listfile) as f:
        ids = {line.strip() for line in f if line.strip()}

    with open(infa) as f, open(outfa, 'w') as o:
        name = None
        seq = []
        for line in f:
            if line.startswith('>'):
                if name is not None and name in ids:
                    o.write(f">{name}\n")
                    s = ''.join(seq)
                    for i in range(0, len(s), 60):
                        o.write(s[i:i+60] + "\n")
                name = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if name is not None and name in ids:
            o.write(f">{name}\n")
            s = ''.join(seq)
            for i in range(0, len(s), 60):
                o.write(s[i:i+60] + "\n")


def _read_fasta_records(path):
    """Read a FASTA file and yield (name, sequence) tuples."""
    name = None
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, ''.join(seq)


def _write_fasta(records, path, wrap=60):
    """Write list of (name, seq) to a FASTA file with line wrapping."""
    with open(path, 'w') as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i:i+wrap] + "\n")


def _fasta_sort_minlen(infa, outfa, prefix="contig", minlen=0):
    """
    Sort FASTA by contig length (descending), rename as prefix_1, prefix_2, ...,
    and filter out contigs shorter than minlen.

    Replaces: funannotate sort -i <in> -b <prefix> -o <out> --minlen <N>
    """
    if not os.path.isfile(infa) or os.path.getsize(infa) == 0:
        with open(outfa, 'w') as f:
            pass
        return 0

    recs = [(name, seq) for name, seq in _read_fasta_records(infa)
            if len(seq) >= minlen]
    recs.sort(key=lambda x: len(x[1]), reverse=True)

    renamed = [(f"{prefix}_{i}", seq) for i, (_, seq) in enumerate(recs, start=1)]
    _write_fasta(renamed, outfa)
    return len(renamed)


def _fasta_clean_contained(infa, outfa, pct_cov=30, exhaustive=True, runner=None):
    """
    Remove contigs that are mostly contained within larger contigs, using
    minimap2 self-alignment.  Keeps the longer contig when two overlap.

    Replaces: funannotate clean -i <in> -p <pct_cov> -o <out> [--exhaustive]

    Algorithm:
      1. Self-align with minimap2 -x asm5 -DP
      2. For each alignment, if the smaller query is covered >= pct_cov%
         by the larger target, mark the query for removal.
      3. Write surviving contigs to outfa.
      4. If exhaustive, repeat until no more removals.
    """
    if not os.path.isfile(infa) or os.path.getsize(infa) == 0:
        with open(outfa, 'w') as f:
            pass
        return

    if not shutil.which("minimap2"):
        if runner:
            runner.log_warn("minimap2 not found; skipping contig cleaning (funannotate clean replacement)")
        shutil.copy(infa, outfa)
        return

    current = infa
    round_num = 0
    max_rounds = 20  # safety limit for exhaustive mode

    while True:
        round_num += 1
        recs = {name: seq for name, seq in _read_fasta_records(current)}
        if len(recs) <= 1:
            break

        # Self-align with minimap2
        with tempfile.NamedTemporaryFile(suffix=".paf", delete=False, mode='w') as paf_tmp:
            paf_path = paf_tmp.name

        try:
            cmd = f"minimap2 -x asm5 -DP --no-long-join -r 500 -c {current} {current}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            with open(paf_path, 'w') as f:
                f.write(result.stdout)
        except Exception:
            break

        # Parse PAF to find contained contigs
        remove = set()
        for line in open(paf_path):
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            qname, qlen, qstart, qend = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
            tname, tlen = parts[5], int(parts[6])
            if qname == tname:
                continue  # skip self-hits

            q_cov = (qend - qstart) * 100.0 / qlen if qlen > 0 else 0

            # If the query is shorter/equal and sufficiently covered, mark for removal
            if qlen <= tlen and q_cov >= pct_cov:
                remove.add(qname)

        os.unlink(paf_path)

        if not remove:
            break

        if runner:
            runner.log(f"  Cleaning round {round_num}: removing {len(remove)} contained contig(s)")

        # Write survivors to temp file for next round
        survivors = [(n, s) for n, s in _read_fasta_records(current) if n not in remove]
        tmp_out = current + f".clean_round{round_num}.tmp"
        _write_fasta(survivors, tmp_out)

        # Clean up previous temp (if not the original input)
        if current != infa and os.path.isfile(current):
            os.unlink(current)
        current = tmp_out

        if not exhaustive:
            break
        if round_num >= max_rounds:
            break

    # Move final result to outfa
    if current != infa:
        shutil.move(current, outfa)
    else:
        shutil.copy(current, outfa)


def _build_busco_csv(runner):
    """Build assembly.busco.csv from BUSCO results."""
    out_csv = os.path.join("assemblies", "assembly.busco.csv")
    desired = ["canu", "reference", "flye", "ipa", "nextDenovo", "peregrine", "hifiasm"]

    def newest(paths):
        return max(paths, key=os.path.getmtime) if paths else None

    def find_run_root(base):
        candidates = []
        p1 = os.path.join("busco", base)
        p2 = os.path.join("busco", f"run_{base}")
        if os.path.isdir(p1):
            candidates.append(p1)
        if os.path.isdir(p2):
            candidates.append(p2)
        if not candidates:
            candidates = sorted(glob.glob(os.path.join("busco", f"{base}*"))) + \
                         sorted(glob.glob(os.path.join("busco", f"run_{base}*")))
        return newest(candidates)

    def read_counts_from_anywhere(run_root):
        if not run_root:
            return None
        fulls = glob.glob(os.path.join(run_root, "**", "full_table*.tsv"), recursive=True)
        if fulls:
            p = newest(fulls)
            with open(p, newline="") as f:
                lines = [ln.rstrip("\n") for ln in f]
            status_idx = None
            for ln in lines:
                if ln.startswith("#") and "Status" in ln:
                    hdr = ln.lstrip("#").strip().split("\t")
                    for i, h in enumerate(hdr):
                        if h.strip().lower() == "status":
                            status_idx = i
                    break
            S = D = F = M = 0
            n = 0
            for ln in lines:
                if not ln or ln.startswith("#"):
                    continue
                parts = ln.split("\t")
                st = parts[1].strip() if status_idx is None else (parts[status_idx].strip() if len(parts) > status_idx else "")
                n += 1
                if st == "Complete":
                    S += 1
                elif st == "Duplicated":
                    D += 1
                elif st == "Fragmented":
                    F += 1
                elif st == "Missing":
                    M += 1
            if n == 0:
                return None
            return {"S": S, "D": D, "F": F, "M": M, "n": n}
        sums = glob.glob(os.path.join(run_root, "**", "short_summary*.txt"), recursive=True)
        if sums:
            txt = open(newest(sums), "r", errors="ignore").read()
            m = re.search(r"C:(\d+(?:\.\d+)?)%.*?S:(\d+(?:\.\d+)?)%.*?D:(\d+(?:\.\d+)?)%.*?F:(\d+(?:\.\d+)?)%.*?M:(\d+(?:\.\d+)?)%.*?n:(\d+)", txt, re.S)
            if not m:
                return None
            Cpct, Spct, Dpct, Fpct, Mpct = map(float, m.groups()[:5])
            n = int(m.group(6))
            Ccnt = round(n * Cpct / 100.0)
            Mcnt = round(n * Mpct / 100.0)
            return {"Cpct": Cpct, "Spct": Spct, "Dpct": Dpct, "Fpct": Fpct, "Mpct": Mpct, "Ccnt": Ccnt, "Mcnt": Mcnt, "n": n}
        return None

    def safe_pct(x, n):
        return 0.0 if not n else 100.0 * float(x) / float(n)

    def fmt_pct(x):
        s = f"{x:.1f}"
        return s.rstrip('0').rstrip('.') if '.' in s else s

    results = {}
    for asm in desired:
        run_root = find_run_root(asm)
        counts = read_counts_from_anywhere(run_root)
        if counts and "S" in counts:
            S, D, F, M, n = counts["S"], counts["D"], counts["F"], counts["M"], counts["n"]
            C = S + D
            results[asm] = {
                "Cpct": fmt_pct(safe_pct(C, n)),
                "Spct": fmt_pct(safe_pct(S, n)),
                "Dpct": fmt_pct(safe_pct(D, n)),
                "Fpct": fmt_pct(safe_pct(F, n)),
                "Mpct": fmt_pct(safe_pct(M, n)),
                "Ccnt": str(C), "Mcnt": str(M), "n": str(n),
            }
        elif counts and "Cpct" in counts:
            results[asm] = {
                "Cpct": fmt_pct(counts["Cpct"]), "Spct": fmt_pct(counts["Spct"]),
                "Dpct": fmt_pct(counts["Dpct"]), "Fpct": fmt_pct(counts["Fpct"]),
                "Mpct": fmt_pct(counts["Mpct"]),
                "Ccnt": str(counts["Ccnt"]), "Mcnt": str(counts["Mcnt"]), "n": str(counts["n"]),
            }
        else:
            results[asm] = {"Cpct": "0", "Spct": "0", "Dpct": "0", "Fpct": "0", "Mpct": "0", "Ccnt": "0", "Mcnt": "0", "n": "0"}

    rows = []
    header = ["Metric"] + desired
    rows.append(header)
    rows.append(["BUSCO C (%)"] + [results[a]["Cpct"] for a in desired])
    rows.append(["BUSCO S (%)"] + [results[a]["Spct"] for a in desired])
    rows.append(["BUSCO D (%)"] + [results[a]["Dpct"] for a in desired])
    rows.append(["BUSCO F (%)"] + [results[a]["Fpct"] for a in desired])
    rows.append(["BUSCO M (%)"] + [results[a]["Mpct"] for a in desired])
    rows.append(["BUSCO C (count)"] + [results[a]["Ccnt"] for a in desired])
    rows.append(["BUSCO M (count)"] + [results[a]["Mcnt"] for a in desired])

    with open(out_csv, "w", newline="") as f:
        csv.writer(f).writerows(rows)

    runner.log(f"Wrote {out_csv}")


def _build_quast_csv(runner):
    """Build assembly.quast.csv from QUAST results."""
    out = os.path.join("assemblies", "assembly.quast.csv")
    treport = os.path.join("quast_out", "transposed_report.tsv")
    report = os.path.join("quast_out", "report.tsv")
    desired = ["canu", "reference", "flye", "ipa", "nextDenovo", "peregrine", "hifiasm"]

    def norm(s):
        return (s or "").lower().replace("-", "").replace("_", "").replace(".result", "").replace(".fasta", "").strip()

    def build_from_report(pth):
        with open(pth, newline="") as f:
            rows = list(csv.reader(f, delimiter="\t"))
        if not rows:
            return None
        header = rows[0]
        asm_headers = header[1:]
        n_asm = len(asm_headers)

        mapcol = {}
        used = set()
        for d in desired:
            nd = norm(d)
            idx = None
            for j, h in enumerate(asm_headers):
                if j in used:
                    continue
                if nd in (norm(h),) or nd == norm(h) or norm(h).startswith(nd) or nd in norm(h):
                    idx = j
                    break
            if idx is None and nd == "hifiasm":
                for j, h in enumerate(asm_headers):
                    nh = norm(h)
                    if "hifiasm" in nh:
                        idx = j
                        break
            mapcol[d] = idx
            if idx is not None:
                used.add(idx)

        out = [["Metric"] + desired]
        for r in rows[1:]:
            if not r:
                continue
            metric = r[0]
            r = r + [""] * (1 + n_asm - len(r))
            vals = []
            for d in desired:
                j = mapcol[d]
                vals.append(r[1 + j] if j is not None else "")
            out.append([metric] + vals)
        return out

    def build_from_transposed(pth):
        with open(pth, newline="") as f:
            rows = list(csv.reader(f, delimiter="\t"))
        if not rows:
            return None
        header = rows[0]
        metrics = header[1:]
        maprow = {d: None for d in desired}
        for i in range(1, len(rows)):
            asm_name = rows[i][0] if rows[i] else ""
            na = norm(asm_name)
            for d in desired:
                nd = norm(d)
                if maprow[d] is not None:
                    continue
                if nd == na or nd in na or na in nd:
                    maprow[d] = i
                elif nd == "hifiasm" and ("hifiasm" in na):
                    maprow[d] = i

        out = [["Metric"] + desired]
        for j, m in enumerate(metrics, start=1):
            line = [m]
            for d in desired:
                i = maprow[d]
                if i is None:
                    line.append("")
                else:
                    row = rows[i]
                    row = row + [""] * (j + 1 - len(row))
                    line.append(row[j] if j < len(row) else "")
            out.append(line)
        return out

    rows = None
    if os.path.exists(report):
        rows = build_from_report(report)
    elif os.path.exists(treport):
        rows = build_from_transposed(treport)
    else:
        runner.log_error("Could not find quast_out/report.tsv or quast_out/transposed_report.tsv")
        return

    if rows:
        with open(out, "w", newline="") as f:
            csv.writer(f).writerows(rows)
        runner.log(f"Wrote {out}")


def _build_assembly_info(runner):
    """Build assembly_info.csv from BUSCO, QUAST, and telomere metrics."""
    os.makedirs("assemblies", exist_ok=True)
    info_csv = "assemblies/assembly_info.csv"
    desired_header = "Metric,canu,reference,flye,ipa,nextDenovo,peregrine,hifiasm"

    with open(info_csv, "w") as f:
        f.write(desired_header + "\n")

    candidates = [
        "assemblies/assembly.busco.csv",
        "assemblies/assembly.quast.csv",
        "assemblies/assembly.telo.csv",
    ]

    for candidate in candidates:
        if os.path.isfile(candidate) and os.path.getsize(candidate) > 0:
            with open(candidate) as f:
                lines = [ln.rstrip('\r\n') for ln in f]

            if lines:
                header = lines[0].split(',')
                desired = desired_header.split(',')

                for line in lines[1:]:
                    if not line:
                        continue
                    parts = line.split(',')
                    out_row = [parts[0] if parts else ""]
                    for i, col in enumerate(desired[1:]):
                        try:
                            idx = header.index(col)
                            out_row.append(parts[idx] if idx < len(parts) else "")
                        except ValueError:
                            out_row.append("")

                    with open(info_csv, "a") as f:
                        f.write(",".join(out_row) + "\n")

    runner.log(f"Wrote {info_csv}")


def _assembler_skip(runner, step_num, name, reason):
    """Log a standard skip message and return."""
    runner.log_warn(f"Step {step_num}: {name} skipped — {reason}")


def step_01_canu(runner):
    """Step 1 - Assembly using HiCanu (non-fatal)."""
    runner.log("Step 1 - Assembly of the genome using HiCanu")

    if not shutil.which("canu"):
        _assembler_skip(runner, 1, "canu", "binary not found. Install a stable release from "
                        "https://github.com/marbl/canu/releases and place it on PATH.")
        return

    # Platform flag
    flag_map = {"pacbio-hifi": "-pacbio-hifi", "nanopore": "-nanopore", "pacbio": "-pacbio"}
    canu_flag = flag_map.get(runner.platform, "-pacbio-hifi")

    # Warn about development builds (broken Java in some conda-installed builds)
    ver = runner.log_version("canu", "canu")
    if ver and ("master" in ver or "changes" in ver):
        runner.log_warn(f"Canu appears to be a development build: {ver}")
        runner.log_warn("Development builds may have broken Java dependencies. "
                        "Install a stable release from https://github.com/marbl/canu/releases")

    cmd = (f"canu -p canu -d hicanu genomeSize={runner.genomesize} "
           f"maxThreads={runner.threads} {canu_flag} {runner.fastq}")
    result = runner.run_cmd(cmd, desc="Running canu", check=False)
    if result.returncode != 0:
        runner.log_warn("Step 1: canu failed. Skipping. Other assemblers will continue. "
                        "Check logs/step_1.log for details.")


def step_02_nextdenovo(runner):
    """Step 2 - Assembly using NextDenovo (non-fatal)."""
    runner.log("Step 2 - Assembly of the genome using NextDenovo")

    if not shutil.which("nextDenovo"):
        _assembler_skip(runner, 2, "nextDenovo", "binary not found. Install via: conda install -c bioconda nextdenovo")
        return

    runner.log_version("nextDenovo", "nextDenovo")
    cmd = f"nextDenovo run_{runner.project}.cfg"
    result = runner.run_cmd(cmd, desc="Running nextDenovo", check=False)
    if result.returncode != 0:
        runner.log_warn("Step 2: nextDenovo failed. Skipping. Check logs/step_2.log for details.")


def step_03_peregrine(runner):
    """Step 3 - Assembly using Peregrine (non-fatal; skipped for Nanopore)."""
    runner.log("Step 3 - Assembly of the genome using Peregrine")

    # Peregrine does not support Nanopore reads
    if runner.platform == "nanopore":
        _assembler_skip(runner, 3, "peregrine", "does not support Nanopore reads")
        return

    if not shutil.which("pg_asm"):
        _assembler_skip(runner, 3, "peregrine", "pg_asm not found. Install peregrine-2021 from "
                        "https://github.com/cschin/peregrine-2021")
        return

    runner.log_version("pg_asm", "pg_asm")
    cmd = f"pg_asm reads_{runner.project}.lst peregrine-2021"
    result = runner.run_cmd(cmd, desc="Running peregrine", check=False)
    if result.returncode != 0:
        runner.log_warn("Step 3: peregrine failed. Skipping. Check logs/step_3.log for details.")


def step_04_ipa(runner):
    """Step 4 - Assembly using IPA (non-fatal; HiFi only)."""
    runner.log("Step 4 - Assembly of the genome using IPA")

    # IPA only supports PacBio HiFi reads
    if runner.platform != "pacbio-hifi":
        _assembler_skip(runner, 4, "IPA", f"only supports PacBio HiFi reads (platform={runner.platform})")
        return

    if not shutil.which("ipa"):
        _assembler_skip(runner, 4, "IPA", "binary not found. Install via: conda install -c bioconda pbipa")
        return

    # IPA calls snakemake internally with --reason, which was removed in v8
    smk = shutil.which("snakemake")
    if smk:
        try:
            sv = subprocess.run(f"{smk} --version", shell=True, capture_output=True, text=True)
            smk_ver = sv.stdout.strip()
            if smk_ver and int(smk_ver.split(".")[0]) >= 8:
                runner.log_warn(f"Snakemake {smk_ver} detected — IPA requires Snakemake 7.x "
                                "(v8 removed the --reason flag that IPA depends on). "
                                "Fix: conda install 'snakemake>=7,<8'")
                _assembler_skip(runner, 4, "IPA", "incompatible Snakemake version")
                return
        except (ValueError, IndexError):
            pass  # couldn't parse version, let IPA try anyway

    runner.log_version("ipa", "ipa")
    if os.path.isdir("ipa"):
        runner.log("Removing existing IPA run directory")
        shutil.rmtree("ipa")
    cmd = f"ipa local --nthreads {runner.threads} --njobs 1 --run-dir ipa -i {runner.fastq}"
    result = runner.run_cmd(cmd, desc="Running ipa", check=False)
    if result.returncode != 0:
        runner.log_warn("Step 4: IPA failed. Skipping. Check logs/step_4.log for details.")


def step_05_flye(runner):
    """Step 5 - Assembly using Flye (non-fatal)."""
    runner.log("Step 5 - Assembly of the genome using Flye")

    if not shutil.which("flye"):
        _assembler_skip(runner, 5, "flye", "binary not found. Install via: conda install -c bioconda flye")
        return

    # Flye read-type flags per platform
    # --nano-hq: ONT reads basecalled with Guppy5+/Dorado (Q20+) — default for
    #   modern ONT data.  For older pre-Q20 ONT reads, users should set
    #   FLYE_ONT_FLAG=--nano-raw in the environment.
    # --pacbio-raw: PacBio CLR (continuous long reads, higher error rate)
    # --pacbio-hifi: PacBio CCS/HiFi reads
    ont_flag = os.environ.get("FLYE_ONT_FLAG", "--nano-hq")
    flag_map = {"pacbio-hifi": "--pacbio-hifi", "nanopore": ont_flag, "pacbio": "--pacbio-raw"}
    flye_flag = flag_map.get(runner.platform, "--pacbio-hifi")

    runner.log_version("flye", "flye")
    cmd = f"flye {flye_flag} {runner.fastq} --out-dir flye --threads {runner.threads}"
    result = runner.run_cmd(cmd, desc="Running flye", check=False)
    if result.returncode != 0:
        runner.log_warn("Step 5: flye failed. Skipping. Check logs/step_5.log for details.")


def step_06_hifiasm(runner):
    """Step 6 - Assembly using Hifiasm (non-fatal; HiFi only).

    hifiasm is designed for PacBio HiFi reads. It does not natively support
    ONT or PacBio CLR reads as primary input. For those platforms, use Flye,
    Canu, or NextDenovo instead.
    """
    runner.log("Step 6 - Assembly of the genome using Hifiasm")

    # hifiasm only supports PacBio HiFi as primary input
    if runner.platform == "nanopore":
        _assembler_skip(runner, 6, "hifiasm",
                        "does not support Nanopore reads as primary input")
        return
    if runner.platform == "pacbio":
        _assembler_skip(runner, 6, "hifiasm",
                        "does not support PacBio CLR reads as primary input")
        return

    if not shutil.which("hifiasm"):
        _assembler_skip(runner, 6, "hifiasm", "binary not found. Install via: conda install -c bioconda hifiasm")
        return

    os.makedirs("hifiasm", exist_ok=True)
    runner.log_version("hifiasm", "hifiasm")

    # HiFi-only: use default mode (diploid-aware); add -l0 for known haploids
    hifiasm_flags = ""

    cmd = f"cd hifiasm && hifiasm -o hifiasm.asm -t {runner.threads} {hifiasm_flags} {runner.fastq} 2> hifiasm.log"
    result = runner.run_cmd(cmd, desc="Running hifiasm", check=False)
    if result.returncode != 0:
        runner.log_warn("Step 6: hifiasm failed. Skipping. Check logs/step_6.log for details.")
        return

    if os.path.isfile("hifiasm/hifiasm.asm.bp.p_ctg.gfa"):
        cmd = "cd hifiasm && awk '/^S/{print \">\"$2; print $3}' hifiasm.asm.bp.p_ctg.gfa > hifiasm.fasta"
        result = runner.run_cmd(cmd, desc="Converting GFA to FASTA", check=False)
        if result.returncode != 0:
            runner.log_warn("Step 6: GFA-to-FASTA conversion failed. Check logs/step_6.log.")
    else:
        runner.log_warn("Step 6: hifiasm primary contig GFA not found (hifiasm.asm.bp.p_ctg.gfa). "
                        "Assembly may have failed or produced no output.")


def step_07_normalize(runner):
    """Step 7 - Copy and normalize all assemblies."""
    runner.log("Step 7 - Copy all assemblies")
    os.makedirs("assemblies", exist_ok=True)

    assembler_paths = [
        ("canu", "./hicanu/canu.contigs.fasta"),
        ("nextDenovo", "./NextDenovo/03.ctg_graph/nd.asm.fasta"),
        ("peregrine", "./peregrine-2021/asm_ctgs_m_p.fa"),
        ("ipa", "./ipa/assembly-results/final.p_ctg.fasta"),
        ("flye", "./flye/assembly.fasta"),
        ("hifiasm", "./hifiasm/hifiasm.fasta"),
    ]

    for prefix, src_path in assembler_paths:
        if os.path.isfile(src_path) and os.path.getsize(src_path) > 0:
            dest = f"./assemblies/{prefix}.result.fasta"
            shutil.copy(src_path, dest)
            tmp_renamed = f"assemblies/.{prefix}.renamed.tmp.fasta"
            rename_and_sort_fasta(runner, dest, tmp_renamed, prefix)
            shutil.move(tmp_renamed, dest)

    if runner.reference_fasta and os.path.isfile(runner.reference_fasta) and os.path.getsize(runner.reference_fasta) > 0:
        shutil.copy(runner.reference_fasta, "./assemblies/reference.result.fasta")


def step_08_busco(runner):
    """Step 8 - Run BUSCO on all assembled genomes."""
    runner.log("Step 8 - Run BUSCO on all assembled genomes (including reference)")
    os.makedirs("busco", exist_ok=True)
    os.makedirs("assemblies", exist_ok=True)

    assemblies = glob.glob("assemblies/*.result.fasta")
    if not assemblies:
        runner.log_error("No assemblies in ./assemblies for BUSCO.")
        raise RuntimeError("No assemblies found")

    lineage = runner.busco_lineage or "ascomycota_odb10"

    def has_busco_metrics(base):
        for d in [f"busco/{base}", f"busco/run_{base}"]:
            if os.path.isdir(d):
                files = glob.glob(f"{d}/**/short_summary*.txt", recursive=True) + \
                        glob.glob(f"{d}/**/full_table*.tsv", recursive=True)
                if files:
                    return True
        return False

    busco_available = shutil.which("busco") is not None

    for fasta_path in assemblies:
        base = os.path.basename(fasta_path).replace(".result.fasta", "").replace(".fasta", "")
        if has_busco_metrics(base):
            runner.log(f"Found existing BUSCO metrics for {base} → using busco/{base} (or legacy path)")
            continue

        if not busco_available:
            runner.log_error(f"Missing BUSCO metrics for '{base}' and 'busco' binary not available; cannot run.")
            raise RuntimeError("BUSCO not available")

        force_flag = ""
        if os.path.isdir(f"busco/{base}"):
            runner.log_info(f"Incomplete BUSCO dir detected for {base} → re-running with -f")
            force_flag = "-f"

        runner.log_info(f"BUSCO on {base} (lineage={lineage}, threads={runner.threads})")
        runner.log_version("busco", "busco")

        cwd = os.getcwd()
        os.chdir("busco")
        cmd = f"busco -i ../{fasta_path} -l {lineage} -m genome -c {runner.threads} -o {base} {force_flag}"
        try:
            runner.run_cmd(cmd, desc=f"Running BUSCO on {base}", check=False)
        finally:
            os.chdir(cwd)

    _build_busco_csv(runner)
    runner.log("Wrote assemblies/assembly.busco.csv")


def step_09_telomere(runner):
    """Step 9 - Hybrid telomere detection and scoring on all assemblies."""
    runner.log("Step 9 - Hybrid telomere detection and scoring")
    os.makedirs("assemblies", exist_ok=True)

    existing_assemblies = glob.glob("assemblies/*.result.fasta")
    if not existing_assemblies:
        runner.log_info("No per-assembly *.result.fasta files found; generating from assembler outputs")
        pairs = [
            ("canu", "./hicanu/canu.contigs.fasta"),
            ("nextDenovo", "./NextDenovo/03.ctg_graph/nd.asm.fasta"),
            ("peregrine", "./peregrine-2021/asm_ctgs_m_p.fa"),
            ("ipa", "./ipa/assembly-results/final.p_ctg.fasta"),
            ("flye", "./flye/assembly.fasta"),
            ("hifiasm", "./hifiasm/hifiasm.fasta"),
        ]
        for name, src in pairs:
            if os.path.isfile(src) and os.path.getsize(src) > 0:
                shutil.copy(src, f"./assemblies/{name}.result.fasta")
        if runner.reference_fasta and os.path.isfile(runner.reference_fasta) and os.path.getsize(runner.reference_fasta) > 0:
            shutil.copy(runner.reference_fasta, "./assemblies/reference.result.fasta")

    existing_assemblies = glob.glob("assemblies/*.result.fasta")
    if not existing_assemblies:
        runner.log_error("No assemblies found in ./assemblies. Run step 7 first or supply --reference.")
        raise RuntimeError("No assemblies found")

    cols = ["canu", "reference", "flye", "ipa", "nextDenovo", "peregrine", "hifiasm"]
    tdouble = {}
    tsingle = {}
    tsupported = {}

    for fasta_path in sorted(existing_assemblies):
        asm = os.path.basename(fasta_path).replace(".result.fasta", "").replace(".fasta", "")
        out_prefix = fasta_path.replace(".result.fasta", "").replace(".fasta", "")

        runner.log(f"Running hybrid telomere detection on {asm}")
        try:
            results = detect_telomeres(
                fasta_path,
                mode=runner.telomere_mode,
                user_motif=runner.motif,
                end_window=runner.telo_end_window,
                score_window=runner.telo_score_window,
                kmer_min=runner.telo_kmer_min,
                kmer_max=runner.telo_kmer_max,
                threads=runner.threads,
                taxon=getattr(runner, 'taxon', 'other'),
            )
            write_detection_outputs(results, fasta_path, out_prefix)
        except Exception as e:
            runner.log_warn(f"Hybrid detection failed for {asm}: {e}; skipping")
            results = []

        # Count classifications
        counts = defaultdict(int)
        for r in results:
            counts[r["classification"]] += 1

        tdouble[asm] = counts.get("strict_t2t", 0)
        tsingle[asm] = counts.get("single_tel_strong", 0)
        tsupported[asm] = counts.get("telomere_supported", 0)
        runner.log(f"  {asm}: {tdouble[asm]} strict_t2t, {tsingle[asm]} single_tel_strong, "
                   f"{tsupported[asm]} telomere_supported")

    # Build CSV matrices — metric names must match Step 14/16 for final_result.csv merge
    matrix_csv = "assemblies/assembly.telo.csv"
    with open(matrix_csv, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["Metric"] + cols
        writer.writerow(header)
        writer.writerow(["Telomere strict T2T contigs"] + [tdouble.get(c, 0) for c in cols])
        writer.writerow(["Telomere single-end strong contigs"] + [tsingle.get(c, 0) for c in cols])
        writer.writerow(["Telomere-supported contigs"] + [tdouble.get(c, 0) + tsingle.get(c, 0) + tsupported.get(c, 0) for c in cols])
    runner.log(f"Wrote {os.path.basename(matrix_csv)}")

    total_csv = "assemblies/total_telo.csv"
    with open(total_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Assembler", "Telomere strict T2T contigs", "Telomere single-end strong contigs", "Telomere-supported contigs"])
        for c in cols:
            writer.writerow([c, tdouble.get(c, 0), tsingle.get(c, 0), tdouble.get(c, 0) + tsingle.get(c, 0) + tsupported.get(c, 0)])
    runner.log(f"Wrote {os.path.basename(total_csv)}")


def _validate_quickmerge_t2t(merged_fasta, input_fastas, runner,
                             max_length_ratio=1.3):
    """Validate quickmerge output: only keep contigs that became strict_t2t
    and pass a length sanity check.

    Quickmerge can join two single-end telomere contigs from different assemblers
    to create a complete T2T contig (left-tel from assembler A + right-tel from
    assembler B).  However, it can also create chimeras (joining contigs from
    different chromosomes).

    Validation criteria:
    1. Merged contig must classify as strict_t2t (telomere on BOTH ends) — this
       proves the join actually recovered a missing telomere end.
    2. Merged contig length must be ≤ max_length_ratio × max(input contig lengths).
       A legitimate join overlaps substantially in the middle, so merged length
       should be close to max(input_lengths).  A chimera concatenates two
       chromosomes, so length ≈ sum(input_lengths), which is ~2× and fails.

    Returns list of (name, sequence) tuples for validated T2T contigs.
    """
    if not os.path.isfile(merged_fasta) or os.path.getsize(merged_fasta) == 0:
        return []

    # Get max input contig length across all input FASTAs for this pair
    max_input_len = 0
    for fa in input_fastas:
        if os.path.isfile(fa):
            for _, seq in _read_fasta_records(fa):
                max_input_len = max(max_input_len, len(seq))

    if max_input_len == 0:
        return []

    length_threshold = int(max_input_len * max_length_ratio)

    # Run telomere detection on merged output
    try:
        merged_results = detect_telomeres(
            merged_fasta,
            mode=runner.telomere_mode,
            user_motif=runner.motif,
            end_window=runner.telo_end_window,
            score_window=runner.telo_score_window,
            kmer_min=runner.telo_kmer_min,
            kmer_max=runner.telo_kmer_max,
            threads=runner.threads,
            taxon=getattr(runner, 'taxon', 'other'),
        )
    except Exception as e:
        runner.log_warn(f"Telomere detection on merged output failed: {e}")
        return []

    # Find strict_t2t contigs that pass length filter
    t2t_names = set()
    for r in merged_results:
        if r["classification"] == "strict_t2t":
            if r["length"] <= length_threshold:
                t2t_names.add(r["contig"])
            else:
                runner.log_warn(
                    f"Rejected merged T2T '{r['contig']}' ({r['length']:,} bp) — "
                    f"exceeds {length_threshold:,} bp ({max_length_ratio}× max input "
                    f"{max_input_len:,} bp), likely chimera"
                )

    if not t2t_names:
        return []

    # Extract validated sequences
    validated = []
    for name, seq in _read_fasta_records(merged_fasta):
        if name in t2t_names:
            validated.append((name, seq))

    return validated


def step_10_telomere_pool(runner):
    """Step 10 - Build optimized telomere contig pool.

    TELOMERE-AWARE VALIDATED MERGE STRATEGY:

    The pool is built from original assembler .telo.fasta files PLUS validated
    quickmerge outputs.  Quickmerge can join two single-end telomere contigs
    (left-tel from assembler A + right-tel from assembler B) to create a
    complete T2T contig — this is valuable for recovering missing telomere ends.

    However, quickmerge can also create chimeras (joining contigs from different
    chromosomes).  To prevent chimeras from entering the pool, each quickmerge
    output is validated:
    1. Telomere detection must classify the merged contig as strict_t2t
       (telomere on BOTH ends — proof the join recovered a missing end)
    2. Merged length must be ≤ 1.3× max(input contig lengths) — a legitimate
       join overlaps in the middle (match-mismatch-match pattern), so the
       merged length ≈ max(inputs).  A chimera ≈ sum(inputs), which is ~2×.

    Only validated T2T contigs from quickmerge enter the pool.  The existing
    clustering then deduplicates, preferring the best representative per
    chromosome group.
    """
    runner.log("Step 10 - Build optimized telomere contig pool")
    os.makedirs("assemblies", exist_ok=True)

    fasta_files = glob.glob("assemblies/*.telo.fasta")
    if not fasta_files:
        runner.log_info("No per-assembly *.telo.fasta files found; attempting to continue")
        fasta_files = glob.glob("assemblies/*.result.fasta")
        if not fasta_files:
            runner.log_error("Still no FASTA files after generation.")
            raise RuntimeError("No FASTA files found")

    # ===== Telomere-aware validated quickmerge =====
    # Run quickmerge pairwise, then validate: only accept merged contigs that
    # (a) classify as strict_t2t and (b) pass length sanity check.
    validated_t2t_from_merge = []  # list of (name, seq) tuples
    if len(fasta_files) > 1 and shutil.which("merge_wrapper.py"):
        for old_f in glob.glob("merged_*.fasta"):
            os.remove(old_f)

        for i in range(len(fasta_files)):
            for j in range(i + 1, len(fasta_files)):
                file1 = fasta_files[i]
                file2 = fasta_files[j]
                base1 = os.path.basename(file1).replace(".telo.fasta", "").replace(".fasta", "")
                base2 = os.path.basename(file2).replace(".telo.fasta", "").replace(".fasta", "")

                prefix = f"merged_{base1}_{base2}"
                cmd = f"merge_wrapper.py -l 1000000 {file1} {file2} --prefix {prefix}"
                runner.run_cmd(cmd, desc=f"Merging {base1} and {base2}", check=False)

                # Validate: only keep merged contigs that are strict_t2t + sane length
                merged_fa = f"{prefix}.fasta"
                if not os.path.isfile(merged_fa):
                    # quickmerge may output with different naming; check alternatives
                    candidates = glob.glob(f"{prefix}*.fasta")
                    if candidates:
                        merged_fa = candidates[0]

                if os.path.isfile(merged_fa) and os.path.getsize(merged_fa) > 0:
                    new_t2t = _validate_quickmerge_t2t(
                        merged_fa, [file1, file2], runner,
                        max_length_ratio=1.3,
                    )
                    if new_t2t:
                        runner.log(
                            f"Validated {len(new_t2t)} T2T contigs from "
                            f"quickmerge({base1} × {base2})"
                        )
                        validated_t2t_from_merge.extend(new_t2t)
                    else:
                        runner.log_info(
                            f"No validated T2T from quickmerge({base1} × {base2}) — "
                            f"merged contigs either not T2T or failed length check"
                        )

        if validated_t2t_from_merge:
            runner.log(f"Total validated T2T contigs from quickmerge: {len(validated_t2t_from_merge)}")
        else:
            runner.log_info("No validated T2T contigs recovered from any quickmerge pair")
    elif len(fasta_files) <= 1:
        runner.log_info("Only one telo FASTA; skipping quickmerge")
    else:
        runner.log_info("merge_wrapper.py not found; skipping quickmerge")

    # ===== Build pool: original .telo.fasta + validated quickmerge T2T =====
    with open("allmerged_telo.fasta", "w") as out:
        for f in sorted(fasta_files):
            with open(f) as inp:
                out.write(inp.read())
        # Append validated T2T contigs from quickmerge (these are proven to have
        # telomere on both ends and pass length sanity)
        for name, seq in validated_t2t_from_merge:
            out.write(f">{name}_qm_validated\n")
            for k in range(0, len(seq), 60):
                out.write(seq[k:k+60] + "\n")

    n_qm = len(validated_t2t_from_merge)
    runner.log(f"Pool: {len(fasta_files)} assembler .telo.fasta files + {n_qm} validated quickmerge T2T")

    n_sorted = _fasta_sort_minlen("allmerged_telo.fasta", "allmerged_telo_sort.fasta",
                                   prefix="contig", minlen=500)
    runner.log(f"Sorted telomere contigs: {n_sorted} contigs >= 500 bp")

    # ===== Hybrid telomere detection on merged pool =====
    runner.log("Running hybrid telomere detection on merged pool")
    pool_fasta = "allmerged_telo_sort.fasta"
    try:
        pool_results = detect_telomeres(
            pool_fasta,
            mode=runner.telomere_mode,
            user_motif=runner.motif,
            end_window=runner.telo_end_window,
            score_window=runner.telo_score_window,
            kmer_min=runner.telo_kmer_min,
            kmer_max=runner.telo_kmer_max,
            threads=runner.threads,
            taxon=getattr(runner, 'taxon', 'other'),
        )
        write_detection_outputs(pool_results, pool_fasta, "allmerged")
    except Exception as e:
        runner.log_warn(f"Hybrid detection on merged pool failed: {e}")
        pool_results = []

    # Build tier ID lists from classification results
    t2t_ids = sorted(set(r["contig"] for r in pool_results if r["classification"] == "strict_t2t"))
    single_ids = sorted(set(r["contig"] for r in pool_results if r["classification"] == "single_tel_strong"))
    supported_ids = sorted(set(r["contig"] for r in pool_results
                               if r["classification"] in ("strict_t2t", "single_tel_strong", "telomere_supported")))

    # Write .list files
    for fname, ids in [("t2t.list", t2t_ids), ("single_tel.list", single_ids),
                       ("telomere_supported.list", supported_ids)]:
        with open(fname, "w") as f:
            for x in ids:
                f.write(x + "\n")

    # Read merged pool sequences for FASTA extraction
    pool_seqs = dict(_read_fasta_records(pool_fasta)) if os.path.isfile(pool_fasta) else {}

    # Write .fasta files directly from classification
    for fname, ids in [("t2t.fasta", t2t_ids), ("single_tel.fasta", single_ids),
                       ("telomere_supported.fasta", supported_ids)]:
        recs = [(n, pool_seqs[n]) for n in ids if n in pool_seqs]
        _write_fasta(recs, fname)

    runner.log(f"Pool classification: {len(t2t_ids)} t2t, {len(single_ids)} single, {len(supported_ids)} supported")

    # ===== Build assembler BUSCO quality map =====
    # Contigs from assemblers with lower BUSCO D (duplication) are preferred.
    # This prevents contigs from high-duplication assemblers (e.g. canu D=74.9%)
    # from winning the clustering representative selection over contigs from
    # cleaner assemblers (e.g. peregrine D=10%).
    asm_busco_d = {}  # assembler_name -> BUSCO D percentage
    info_csv = os.path.join("assemblies", "assembly_info.csv")
    if os.path.isfile(info_csv):
        try:
            with open(info_csv) as f:
                reader = csv.reader(f)
                header = next(reader, [])
                asm_names = [h.strip() for h in header[1:]]
                for row in reader:
                    if not row:
                        continue
                    metric = row[0].strip().lower()
                    if "busco d (%)" in metric:
                        for i, asm in enumerate(asm_names):
                            try:
                                asm_busco_d[asm.lower()] = float(row[i + 1])
                            except (IndexError, ValueError):
                                pass
                        break
        except Exception as e:
            runner.log_warn(f"Could not read BUSCO D from assembly_info.csv: {e}")

    if asm_busco_d:
        runner.log(f"Assembler BUSCO D quality map: {asm_busco_d}")

    # Build per-contig telomere score map for clustering representative selection.
    # Score = telomere_score * quality_weight, where quality_weight penalises
    # contigs from high-duplication assemblers.
    # Contigs are renamed to "contig_NNN" in the pool, but the original assembler
    # .telo.fasta files have prefixes that indicate which assembler produced them.
    # We track the source assembler for each contig by noting which input file
    # contributed it.
    contig_source_asm = {}  # contig_name -> assembler_name
    for f_path in sorted(fasta_files):
        basename = os.path.basename(f_path).replace(".telo.fasta", "").replace(".result.fasta", "").replace(".fasta", "")
        asm_key = basename.lower()
        if os.path.isfile(f_path):
            for line in open(f_path):
                if line.startswith(">"):
                    # After sort/rename, contigs get sequential "contig_NNN" names.
                    # We can't directly map those back.  Instead we'll count contigs
                    # per assembler and assign quality based on position in the
                    # concatenated pool.
                    pass  # filled below via position tracking

    # Track source assembler by position in concatenated pool
    contig_order = []  # list of (assembler_key, original_name)
    for f_path in sorted(fasta_files):
        basename = os.path.basename(f_path).replace(".telo.fasta", "").replace(".result.fasta", "").replace(".fasta", "")
        asm_key = basename.lower()
        if os.path.isfile(f_path):
            for line in open(f_path):
                if line.startswith(">"):
                    orig_name = line[1:].strip().split()[0]
                    contig_order.append((asm_key, orig_name))

    # After funannotate sort / _fasta_sort_minlen, contigs are renamed to
    # "contig_001", "contig_002", etc. sorted by length descending.
    # We need to map pool contig names back to their source assembler.
    # Read the sorted pool and match by sequence identity to the originals.
    # SIMPLER APPROACH: read pool_fasta, for each contig check which assembler
    # .telo.fasta contains a sequence with the same length (approximate match).
    # EVEN SIMPLER: read all original seqs with lengths, map pool contigs by length.
    orig_seq_lengths = {}  # (asm_key, length) -> asm_key
    for f_path in sorted(fasta_files):
        basename = os.path.basename(f_path).replace(".telo.fasta", "").replace(".result.fasta", "").replace(".fasta", "")
        asm_key = basename.lower()
        if os.path.isfile(f_path):
            name = None
            seqlen = 0
            for line in open(f_path):
                if line.startswith(">"):
                    if name is not None:
                        orig_seq_lengths[(asm_key, seqlen)] = asm_key
                    name = line[1:].strip().split()[0]
                    seqlen = 0
                else:
                    seqlen += len(line.strip())
            if name is not None:
                orig_seq_lengths[(asm_key, seqlen)] = asm_key

    # Map pool contigs to assembler by matching sequence length
    pool_contig_asm = {}  # pool_contig_name -> asm_key
    if os.path.isfile(pool_fasta):
        pool_seqs_for_map = dict(_read_fasta_records(pool_fasta))
        # Build length->asm lookup (group by length)
        len_to_asm = {}
        for (asm_key, slen), ak in orig_seq_lengths.items():
            len_to_asm.setdefault(slen, []).append(ak)
        for cname, cseq in pool_seqs_for_map.items():
            clen = len(cseq)
            candidates = len_to_asm.get(clen, [])
            if candidates:
                pool_contig_asm[cname] = candidates[0]  # first match

    if pool_contig_asm:
        # Count per assembler
        asm_counts = {}
        for cn, ak in pool_contig_asm.items():
            asm_counts[ak] = asm_counts.get(ak, 0) + 1
        runner.log(f"Pool contig source mapping: {asm_counts}")

    telo_scores = {}
    for r in pool_results:
        ls = r.get("left_score", 0.0) or 0.0
        rs = r.get("right_score", 0.0) or 0.0
        base_score = max(ls, rs)

        # Apply quality weight: penalise contigs from high-D assemblers
        cname = r["contig"]
        asm_key = pool_contig_asm.get(cname, "")
        busco_d = asm_busco_d.get(asm_key, 0.0)
        # Quality weight: 1.0 for D=0%, 0.25 for D=100%
        quality_weight = max(0.25, 1.0 - busco_d / 133.0)
        telo_scores[cname] = base_score * quality_weight

    # ===== Redundancy reduction for T2T contigs via minimap2 clustering =====
    if os.path.isfile("t2t.fasta") and os.path.getsize("t2t.fasta") > 0 and shutil.which("minimap2"):
        runner.log("Selecting best T2T representative per chromosome via minimap2 clustering")
        cmd = f"minimap2 -x asm20 -D -P -t {runner.threads} t2t.fasta t2t.fasta"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        with open("t2t.self.paf", "w") as f:
            f.write(result.stdout)
        cluster_and_select("t2t.fasta", "t2t.self.paf", "t2t_best.fasta",
                           "t2t_cluster_summary.tsv", label="T2T",
                           min_identity=0.95, min_coverage=0.85,
                           scores=telo_scores)
        shutil.copy("t2t_best.fasta", "t2t.fasta")
    else:
        runner.log_info("T2T clustering skipped (no T2T contigs or minimap2 not available)")

    # ===== Redundancy reduction for single-end contigs via minimap2 clustering =====
    if os.path.isfile("single_tel.fasta") and os.path.getsize("single_tel.fasta") > 0 and shutil.which("minimap2"):
        runner.log("Selecting best single-end telomeric representatives via minimap2 clustering")
        cmd = f"minimap2 -x asm20 -D -P -t {runner.threads} single_tel.fasta single_tel.fasta"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        with open("single_tel.self.paf", "w") as f:
            f.write(result.stdout)
        cluster_and_select("single_tel.fasta", "single_tel.self.paf", "single_tel_best.fasta",
                           "telomere_cluster_summary.tsv", label="single-end",
                           min_identity=0.95, min_coverage=0.85,
                           scores=telo_scores)
    elif os.path.isfile("single_tel.fasta") and os.path.getsize("single_tel.fasta") > 0:
        shutil.copy("single_tel.fasta", "single_tel_best.fasta")
    else:
        with open("single_tel_best.fasta", "w") as f:
            pass

    # Combine outputs
    with open("telomere_supported_best.fasta", "w") as out:
        if os.path.isfile("t2t.fasta") and os.path.getsize("t2t.fasta") > 0:
            with open("t2t.fasta") as f:
                out.write(f.read())
        if os.path.isfile("single_tel_best.fasta") and os.path.getsize("single_tel_best.fasta") > 0:
            with open("single_tel_best.fasta") as f:
                out.write(f.read())

    # Clean contained/duplicate contigs (replaces funannotate clean).
    # IMPORTANT: T2T contigs are NOT cleaned — the minimap2 clustering already
    # produced one best representative per chromosome.  Running clean_contained
    # on T2T can incorrectly remove T2T contigs that share repetitive elements
    # with longer contigs on different chromosomes (16→15 loss observed).
    if os.path.isfile("t2t.fasta") and os.path.getsize("t2t.fasta") > 0:
        shutil.copy("t2t.fasta", "t2t_clean.fasta")
        runner.log("T2T contigs preserved without clean_contained (clustering already deduplicated)")
    else:
        with open("t2t_clean.fasta", "w") as f:
            pass

    for src, dst in [("single_tel_best.fasta", "single_tel_best_clean.fasta"),
                     ("telomere_supported_best.fasta", "telomere_supported_best_clean.fasta")]:
        if os.path.isfile(src) and os.path.getsize(src) > 0:
            runner.log(f"Cleaning contained contigs from {src}")
            _fasta_clean_contained(src, dst, pct_cov=30, exhaustive=True, runner=runner)
        else:
            with open(dst, "w") as f:
                pass

    # Update references
    for src, dst in [("single_tel_best.fasta", "single_tel.fasta"),
                     ("single_tel_best_clean.fasta", "single_tel_clean.fasta"),
                     ("telomere_supported_best.fasta", "telomere_supported.fasta"),
                     ("telomere_supported_best_clean.fasta", "telomere_supported_clean.fasta")]:
        if os.path.isfile(src):
            shutil.copy(src, dst)

    # Count contigs
    strict_t2t_n = sum(1 for _ in open("t2t.fasta") if _.startswith(">")) if os.path.isfile("t2t.fasta") else 0
    single_tel_n = sum(1 for _ in open("single_tel_best.fasta") if _.startswith(">")) if os.path.isfile("single_tel_best.fasta") else 0
    tel_supported_n = sum(1 for _ in open("telomere_supported_best.fasta") if _.startswith(">")) if os.path.isfile("telomere_supported_best.fasta") else 0

    with open("telomere_support_summary.csv", "w") as f:
        f.write(f"strict_t2t_contigs,{strict_t2t_n}\n")
        f.write(f"single_telomere_best_contigs,{single_tel_n}\n")
        f.write(f"telomere_supported_best_contigs,{tel_supported_n}\n")

    runner.log("Telomere support summary:")
    with open("telomere_support_summary.csv") as f:
        runner.log(f.read())

    # Set protected telomere mode
    if os.path.isfile("t2t_clean.fasta") and os.path.getsize("t2t_clean.fasta") > 0:
        shutil.copy("t2t_clean.fasta", "protected_telomere_contigs.fasta")
        with open("protected_telomere_mode.txt", "w") as f:
            f.write("strict_t2t\n")
        runner.log("Protected contigs mode: strict_t2t")
    elif os.path.isfile("single_tel_best_clean.fasta") and os.path.getsize("single_tel_best_clean.fasta") > 0:
        shutil.copy("single_tel_best_clean.fasta", "protected_telomere_contigs.fasta")
        with open("protected_telomere_mode.txt", "w") as f:
            f.write("single_tel_best\n")
        runner.log("Protected contigs mode: single_tel_best")
    elif os.path.isfile("telomere_supported_best_clean.fasta") and os.path.getsize("telomere_supported_best_clean.fasta") > 0:
        shutil.copy("telomere_supported_best_clean.fasta", "protected_telomere_contigs.fasta")
        with open("protected_telomere_mode.txt", "w") as f:
            f.write("telomere_supported_best\n")
        runner.log("Protected contigs mode: telomere_supported_best")
    else:
        with open("protected_telomere_contigs.fasta", "w") as f:
            pass
        with open("protected_telomere_mode.txt", "w") as f:
            f.write("none\n")
        runner.log_warn("No telomere-supported contigs found")


def step_11_quast(runner):
    """Step 11 - QUAST metrics for all assemblies."""
    runner.log("Step 11 - QUAST metrics for all assemblies")
    os.makedirs("quast_out", exist_ok=True)
    os.makedirs("assemblies", exist_ok=True)

    assemblies = glob.glob("assemblies/*.result.fasta")
    if not assemblies:
        runner.log_error("No assemblies to run QUAST.")
        raise RuntimeError("No assemblies found")

    quast_bin = shutil.which("quast.py") or shutil.which("quast")
    if quast_bin:
        asm_str = " ".join(assemblies)
        cmd = f"{quast_bin} {asm_str} --threads {runner.threads} -o quast_out"
        runner.run_cmd(cmd, desc="Running QUAST", check=False)

    _build_quast_csv(runner)
    runner.log("Wrote assemblies/assembly.quast.csv")


def _run_merqury_preselection(runner):
    """Run Merqury on all assembler outputs for pre-selection scoring."""
    db = runner.merqury_db
    if not db:
        for cand in ["reads.meryl", "meryl/reads.meryl", "merqury/reads.meryl"] + glob.glob("*.meryl"):
            if os.path.isdir(cand):
                db = cand
                break

    if db and os.path.isdir(db) and shutil.which("merqury.sh"):
        runner.log_info(f"Running Merqury pre-selection using database: {db}")
        runner.log_version("merqury.sh", "merqury.sh")

        assemblers = ["canu", "reference", "flye", "ipa", "nextDenovo", "peregrine", "hifiasm"]
        for asm in assemblers:
            asm_fa = f"assemblies/{asm}.result.fasta"
            if not os.path.isfile(asm_fa) or os.path.getsize(asm_fa) == 0:
                continue
            qv_file = os.path.join("merqury", f"{asm}.qv")
            comp_file = os.path.join("merqury", f"{asm}.completeness.stats")
            if not os.path.isfile(qv_file) or not os.path.isfile(comp_file):
                cmd = f"merqury.sh {db} {asm_fa} merqury/{asm}"
                runner.run_cmd(cmd, desc=f"Merqury on {asm}", check=False)
    else:
        runner.log_warn("Merqury requested but merqury.sh or a valid .meryl database was not found; skipping.")


def _write_merqury_csv():
    """Write assemblies/assembly.merqury.csv from Merqury output files."""
    assemblers = ["canu", "reference", "flye", "ipa", "nextDenovo", "peregrine", "hifiasm"]
    rows = [["Metric"] + assemblers, ["Merqury QV"], ["Merqury completeness (%)"]]

    def parse_first_float(path):
        if not os.path.exists(path):
            return ""
        txt = open(path, "r", errors="ignore").read()
        for pat in [r'(?i)qv[^0-9]*([0-9]+(?:\.[0-9]+)?)',
                    r'(?i)completeness[^0-9]*([0-9]+(?:\.[0-9]+)?)',
                    r'([0-9]+(?:\.[0-9]+)?)']:
            m = re.search(pat, txt)
            if m:
                return m.group(1)
        return ""

    for asm in assemblers:
        rows[1].append(parse_first_float(os.path.join("merqury", f"{asm}.qv")))
        rows[2].append(parse_first_float(os.path.join("merqury", f"{asm}.completeness.stats")))

    with open("assemblies/assembly.merqury.csv", "w", newline="") as f:
        csv.writer(f).writerows(rows)


def _auto_select_backbone(runner):
    """Auto-select best backbone assembler using smart scoring or N50 mode.

    Reads assembly_info.csv (transposed format: Metric,canu,flye,...) and
    applies the scoring formula matching TACO.sh Step 12B.
    """
    info_csv = "assemblies/assembly_info.csv"
    mode = runner.auto_mode or "smart"
    debug_tsv = "assemblies/selection_debug.tsv"
    decision_txt = "assemblies/selection_decision.txt"

    if not os.path.isfile(info_csv) or os.path.getsize(info_csv) == 0:
        runner.log_warn("assembly_info.csv missing or empty; cannot auto-select assembler.")
        return None

    with open(info_csv, newline="") as f:
        rows = list(csv.reader(f))
    if not rows:
        return None

    header = [h.strip() for h in rows[0]]

    # Build row lookup: metric_name -> row
    metric_rows = {}
    for r in rows[1:]:
        if r and r[0]:
            metric_rows[r[0].strip().lower()] = r

    def get_val(metric_name_lower, col_idx, default=0.0):
        row = metric_rows.get(metric_name_lower)
        if row is None:
            return default
        try:
            return float(row[col_idx])
        except (IndexError, ValueError):
            return default

    best_name = None
    best_score = None
    records = []

    for idx, asm in enumerate(header[1:], start=1):
        asm = asm.strip()
        if not asm:
            continue
        fa = f"assemblies/{asm}.result.fasta"
        if not os.path.isfile(fa) or os.path.getsize(fa) == 0:
            continue

        busco_s = get_val("busco s (%)", idx)
        busco_d = get_val("busco d (%)", idx)
        contigs = get_val("# contigs", idx)
        n50 = get_val("n50", idx)
        t2t = get_val("telomere strict t2t contigs", idx)
        single = get_val("telomere single-end strong contigs", idx)
        merqury_qv = get_val("merqury qv", idx)
        merqury_comp = get_val("merqury completeness (%)", idx)

        if mode == "n50":
            if n50 <= 0:
                continue
            score = n50
        else:
            if contigs <= 0 or n50 <= 0:
                continue

            # Taxon-aware scoring weights
            # - Fungi (small genomes): stronger telomere rescue weight,
            #   high BUSCO_D penalty (duplicated assemblies look falsely good)
            # - Plants: stronger fragmentation penalty, moderate BUSCO_D
            #   (polyploidy can inflate D naturally)
            # - Vertebrate/animal: stronger contiguity preference,
            #   conservative telomere weight (large repeat-rich chromosomes)
            # - Insect/other: balanced defaults
            taxon = getattr(runner, 'taxon', 'other')
            w_busco_s = 1000
            w_t2t = 300
            w_single = 150
            w_contigs = 30
            w_n50 = 150
            w_busco_d = 500

            if taxon == "fungal":
                w_busco_d = 600  # penalize duplication more strongly
                w_t2t = 350      # telomere rescue is strong for small genomes
            elif taxon == "plant":
                w_contigs = 50   # penalize fragmentation more
                w_busco_d = 300  # relax D penalty (polyploidy inflates D)
                w_t2t = 200      # reduce telomere weight (false signals from repeats)
            elif taxon in ("vertebrate", "animal"):
                w_contigs = 40   # penalize fragmentation moderately
                w_n50 = 200      # stronger contiguity preference
                w_t2t = 200      # reduce telomere weight (ITRs common)

            score = (
                busco_s * w_busco_s
                + t2t * w_t2t
                + single * w_single
                + merqury_comp * 200
                + merqury_qv * 20
                - contigs * w_contigs
                + math.log10(n50) * w_n50
                - busco_d * w_busco_d
            )

        records.append({
            "assembler": asm, "busco_s": busco_s, "busco_d": busco_d,
            "t2t": t2t, "single_tel": single, "merqury_qv": merqury_qv,
            "merqury_comp": merqury_comp, "contigs": contigs,
            "n50": n50, "score": score,
        })

        runner.log(f"  {asm}: BUSCO_S={busco_s} BUSCO_D={busco_d} T2T={t2t} single={single} "
                   f"MerquryQV={merqury_qv} MerquryComp={merqury_comp} "
                   f"contigs={contigs} N50={n50} score={score:.1f}")

        if best_score is None or score > best_score:
            best_name = asm
            best_score = score

    # Write debug TSV
    with open(debug_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["assembler", "busco_s", "busco_d", "t2t", "single_tel",
                     "merqury_qv", "merqury_comp", "contigs", "n50", "score"])
        for r in records:
            w.writerow([r["assembler"], r["busco_s"], r["busco_d"], r["t2t"], r["single_tel"],
                        r["merqury_qv"], r["merqury_comp"], r["contigs"], r["n50"], r["score"]])

    # Write decision file
    with open(decision_txt, "w") as f:
        f.write(f"auto_mode\t{mode}\n")
        f.write(f"selected_assembler\t{best_name or ''}\n")
        f.write(f"selected_score\t{best_score if best_score is not None else ''}\n")
        f.write("score_formula\tBUSCO_S*1000 + T2T*300 + single*150 + MerquryComp*200 + MerquryQV*20 - contigs*30 + log10(N50)*150 - BUSCO_D*500\n")

    return best_name


def _clean_backbone_headers(asm_fa, out_fa):
    """Sanitize backbone FASTA headers (remove special chars, deduplicate names)."""
    seen = {}
    recs = []
    for name, seq in _read_fasta_records(asm_fa):
        if not seq:
            continue
        base = name.split()[0]
        base = "".join(c if c.isalnum() or c in "._-" else "_" for c in base)
        if base in seen:
            seen[base] += 1
            new = f"{base}_{seen[base]}"
        else:
            seen[base] = 1
            new = base
        recs.append((new, seq))
    _write_fasta(recs, out_fa)


def _parse_paf_best_hits(paf_path):
    """Parse PAF and return best (coverage, identity) per query contig.

    Returns dict: {qname: (best_cov, best_ident)}
    """
    best = {}
    if os.path.isfile(paf_path):
        with open(paf_path) as f:
            for ln in f:
                if not ln.strip() or ln.startswith("#"):
                    continue
                p = ln.rstrip().split("\t")
                if len(p) < 12:
                    continue
                qname, qlen = p[0], int(p[1])
                qs, qe = int(p[2]), int(p[3])
                matches, alnlen = int(p[9]), int(p[10])
                if qlen <= 0 or alnlen <= 0:
                    continue
                ident = matches / max(1, alnlen)
                cov = (qe - qs) / max(1, qlen)
                cur = best.get(qname, (0.0, 0.0))
                if cov > cur[0] or (abs(cov - cur[0]) < 1e-12 and ident > cur[1]):
                    best[qname] = (cov, ident)
    return best


def _filter_redundant_to_protected(paf_path, backbone_fa, out_fa, cov_thr=0.95, id_thr=0.95):
    """Remove backbone contigs that are redundant to the protected telomere pool.

    Uses PAF from minimap2 alignment of backbone vs protected pool.
    Keeps backbone contigs that are NOT covered >= cov_thr with identity >= id_thr.
    """
    best = _parse_paf_best_hits(paf_path)
    drop = {q for q, (cov, ident) in best.items() if cov >= cov_thr and ident >= id_thr}

    recs = [(n, s) for n, s in _read_fasta_records(backbone_fa) if n not in drop]
    _write_fasta(recs, out_fa)
    return len(drop)


def _filter_fragments_to_protected(paf_path, backbone_fa, out_fa,
                                   frag_cov_thr=0.50, frag_id_thr=0.90):
    """Second-pass fragment removal: drop backbone contigs that are partial
    copies of protected T2T chromosomes.

    The strict 95%/95% pass catches near-identical contigs.  This pass catches
    backbone FRAGMENTS — contigs where ≥50% of their length aligns to a T2T
    contig at ≥90% identity.  These are partial chromosome copies that carry
    duplicated gene content but are too divergent for the strict filter.

    For a genome where all chromosomes already have T2T representatives in the
    protected pool, any backbone contig that substantially aligns to T2T is
    redundant by definition — it represents a worse version of an already-
    complete chromosome.
    """
    best = _parse_paf_best_hits(paf_path)
    drop = {q for q, (cov, ident) in best.items()
            if cov >= frag_cov_thr and ident >= frag_id_thr}

    recs = [(n, s) for n, s in _read_fasta_records(backbone_fa) if n not in drop]
    _write_fasta(recs, out_fa)
    return len(drop)


def _name_dedup_fasta(protected_fa, input_fa, output_fa):
    """Remove contigs from input_fa whose names match protected_fa IDs."""
    protected_ids = set()
    if os.path.isfile(protected_fa):
        for name, _ in _read_fasta_records(protected_fa):
            protected_ids.add(name)

    recs = [(n, s) for n, s in _read_fasta_records(input_fa) if n not in protected_ids]
    _write_fasta(recs, output_fa)
    return len(protected_ids)


def _classify_contigs_telomere_status(fasta_path, runner):
    """Run telomere detection on contigs and return sets of telomere-bearing IDs.

    Returns (t2t_ids, telo_ids) where:
      - t2t_ids: contigs classified as strict_t2t
      - telo_ids: contigs with ANY telomere signal (strict_t2t, single_tel_strong,
        or telomere_supported)
    """
    t2t_ids = set()
    telo_ids = set()

    if not os.path.isfile(fasta_path) or os.path.getsize(fasta_path) == 0:
        return t2t_ids, telo_ids

    try:
        results = detect_telomeres(
            fasta_path,
            mode=runner.telomere_mode,
            user_motif=runner.motif,
            end_window=runner.telo_end_window,
            score_window=runner.telo_score_window,
            kmer_min=runner.telo_kmer_min,
            kmer_max=runner.telo_kmer_max,
            threads=runner.threads,
            taxon=getattr(runner, 'taxon', 'other'),
        )
        for r in results:
            cls = r.get("classification", "")
            name = r.get("contig", "")
            if cls == "strict_t2t":
                t2t_ids.add(name)
                telo_ids.add(name)
            elif cls in ("single_tel_strong", "telomere_supported"):
                telo_ids.add(name)
    except Exception:
        pass

    return t2t_ids, telo_ids


def _self_dedup_non_telomeric(fasta_path, telo_ids, out_path, threads,
                               cov_thr=0.80, id_thr=0.90):
    """Remove non-telomeric backbone contigs that are redundant to each other.

    When two non-telomeric contigs overlap at >= cov_thr / id_thr, the shorter
    one is removed.  Telomere-bearing contigs (in telo_ids) are always kept.
    """
    if not os.path.isfile(fasta_path) or os.path.getsize(fasta_path) == 0:
        _write_fasta([], out_path)
        return 0

    if not shutil.which("minimap2"):
        shutil.copy(fasta_path, out_path)
        return 0

    # Self-align
    paf_path = fasta_path + ".self.paf"
    cmd = f"minimap2 -x asm20 -D -P -t {threads} {fasta_path} {fasta_path}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    with open(paf_path, "w") as f:
        f.write(result.stdout)

    # Parse PAF: for each pair of non-telomeric contigs, mark shorter as drop
    seq_lens = {}
    for n, s in _read_fasta_records(fasta_path):
        seq_lens[n] = len(s)

    drop = set()
    with open(paf_path) as f:
        for ln in f:
            if not ln.strip():
                continue
            p = ln.rstrip().split("\t")
            if len(p) < 12:
                continue
            qname, qlen = p[0], int(p[1])
            qs, qe = int(p[2]), int(p[3])
            tname, tlen = p[5], int(p[6])
            ts, te = int(p[7]), int(p[8])
            matches, alnlen = int(p[9]), int(p[10])

            if qname == tname:
                continue
            if qname in telo_ids or tname in telo_ids:
                continue  # never remove telomere-bearing contigs
            if alnlen <= 0:
                continue

            ident = matches / max(1, alnlen)
            cov_q = (qe - qs) / max(1, qlen)

            if cov_q >= cov_thr and ident >= id_thr:
                # Shorter contig is the redundant one
                if qlen <= tlen:
                    drop.add(qname)
                else:
                    drop.add(tname)

    recs = [(n, s) for n, s in _read_fasta_records(fasta_path) if n not in drop]
    _write_fasta(recs, out_path)

    # Cleanup
    if os.path.isfile(paf_path):
        os.remove(paf_path)

    return len(drop)


def _telomere_rescue(single_fa, backbone_fa, paf_path, out_ids, out_fa,
                     min_ident=0.90, min_cov=0.80, min_ext=500):
    """Replace backbone contigs with telomeric single-end contigs when they extend the backbone.

    For each backbone contig that has a single-end telomeric contig aligned to it with
    sufficient identity, coverage, and terminal extension, replace the backbone contig
    with the telomeric contig.
    """
    single = dict(_read_fasta_records(single_fa))
    backbone = dict(_read_fasta_records(backbone_fa))

    best = {}  # tname -> ((ext, ident, cov, len, qname), qname)
    if os.path.isfile(paf_path):
        for ln in open(paf_path):
            if not ln.strip() or ln.startswith("#"):
                continue
            p = ln.rstrip().split("\t")
            if len(p) < 12:
                continue
            qname, qlen = p[0], int(p[1])
            qs, qe = int(p[2]), int(p[3])
            tname = p[5]
            matches, alnlen = int(p[9]), int(p[10])
            if qname not in single or tname not in backbone:
                continue
            if alnlen <= 0 or qlen <= 0:
                continue
            ident = matches / max(1, alnlen)
            cov_q = (qe - qs) / max(1, qlen)
            ext = max(qs, qlen - qe)  # terminal extension
            if ident < min_ident or cov_q < min_cov or ext < min_ext:
                continue
            score = (ext, ident, cov_q, len(single[qname]), qname)
            cur = best.get(tname)
            if cur is None or score > cur[0]:
                best[tname] = (score, qname)

    replaced = set()
    with open(out_ids, "w") as ids:
        for tname, (_, qname) in sorted(best.items()):
            replaced.add(tname)
            ids.write(f"{tname}\t{qname}\n")

    used_single = set(qname for _, qname in best.values())
    recs = []
    for qname in sorted(used_single):
        if qname in single:
            recs.append((qname, single[qname]))
    for tname in backbone:
        if tname not in replaced:
            recs.append((tname, backbone[tname]))

    _write_fasta(recs, out_fa)
    return len(replaced)


def _parse_paf_rescue_hits(paf_path, donor_seqs, backbone_seqs, terminal_window=500):
    """Parse PAF for telomere rescue candidate screening with detailed metrics.

    Returns list of dicts with per-hit metrics for structural screening.
    """
    hits = []
    if not os.path.isfile(paf_path):
        return hits

    with open(paf_path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue
            p = ln.rstrip().split("\t")
            if len(p) < 12:
                continue
            qname, qlen = p[0], int(p[1])
            qs, qe = int(p[2]), int(p[3])
            strand = p[4]
            tname, tlen = p[5], int(p[6])
            ts, te = int(p[7]), int(p[8])
            matches, alnlen = int(p[9]), int(p[10])

            if qname not in donor_seqs or tname not in backbone_seqs:
                continue
            if alnlen <= 0 or qlen <= 0 or tlen <= 0:
                continue

            ident = matches / max(1, alnlen)
            cov_backbone = (te - ts) / max(1, tlen)
            cov_donor = (qe - qs) / max(1, qlen)
            ext = max(qs, qlen - qe)  # donor extension beyond alignment
            len_gain = qlen - tlen
            touches_left = ts <= terminal_window
            touches_right = (tlen - te) <= terminal_window

            hits.append({
                "donor": qname,
                "backbone": tname,
                "donor_len": qlen,
                "backbone_len": tlen,
                "ident": ident,
                "aligned_bp": alnlen,
                "cov_backbone": cov_backbone,
                "cov_donor": cov_donor,
                "ext": ext,
                "len_gain": len_gain,
                "touches_left": touches_left,
                "touches_right": touches_right,
                "strand": strand,
            })
    return hits


def _screen_rescue_candidates(hits, min_ident=0.85, min_aligned_bp=8000,
                               min_cov_backbone=0.60, min_cov_donor=0.50,
                               min_ext=1000):
    """Split hits into obvious rejections and plausible candidates.

    Returns (rejected, candidates) where each is a list of dicts with
    added 'reject_reason' or 'structural_score' fields.
    """
    rejected = []
    candidates = []

    # Keep only the best hit per (donor, backbone) pair
    best_per_pair = {}
    for h in hits:
        key = (h["donor"], h["backbone"])
        cur = best_per_pair.get(key)
        if cur is None or h["cov_backbone"] > cur["cov_backbone"]:
            best_per_pair[key] = h

    for h in best_per_pair.values():
        reasons = []
        if h["ident"] < min_ident:
            reasons.append(f"ident={h['ident']:.4f}<{min_ident}")
        if h["aligned_bp"] < min_aligned_bp:
            reasons.append(f"aligned_bp={h['aligned_bp']}<{min_aligned_bp}")
        if h["cov_backbone"] < min_cov_backbone:
            reasons.append(f"cov_bb={h['cov_backbone']:.4f}<{min_cov_backbone}")
        if h["cov_donor"] < min_cov_donor:
            reasons.append(f"cov_donor={h['cov_donor']:.4f}<{min_cov_donor}")
        if h["ext"] < min_ext:
            reasons.append(f"ext={h['ext']}<{min_ext}")
        if not h["touches_left"] and not h["touches_right"]:
            reasons.append("no_terminal_touch")

        if reasons:
            h["reject_reason"] = "; ".join(reasons)
            rejected.append(h)
        else:
            # Composite structural score for ranking
            score = (
                h["cov_backbone"] * 3.0 +
                h["cov_donor"] * 2.0 +
                h["ident"] * 2.0 +
                min(h["ext"] / 5000.0, 1.0) * 1.5 +
                (1.0 if h["touches_left"] or h["touches_right"] else 0.0) +
                (0.5 if h["touches_left"] and h["touches_right"] else 0.0) +
                max(min(h["len_gain"] / 10000.0, 0.5), -0.3)
            )
            h["structural_score"] = score
            candidates.append(h)

    # Rank candidates by structural score descending
    candidates.sort(key=lambda x: x["structural_score"], reverse=True)
    for i, c in enumerate(candidates):
        c["candidate_rank"] = i + 1

    return rejected, candidates


def _write_rescue_debug_tsv(hits, path):
    """Write all evaluated hits to debug TSV."""
    cols = ["donor", "backbone", "ident", "aligned_bp", "cov_backbone",
            "cov_donor", "ext", "len_gain", "touches_left", "touches_right"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for h in hits:
            w.writerow([
                h["donor"], h["backbone"],
                f"{h['ident']:.4f}", h["aligned_bp"],
                f"{h['cov_backbone']:.4f}", f"{h['cov_donor']:.4f}",
                h["ext"], h["len_gain"],
                h["touches_left"], h["touches_right"],
            ])


def _write_candidates_tsv(candidates, path):
    """Write plausible rescue candidates to TSV."""
    cols = ["candidate_rank", "backbone", "donor", "ident", "aligned_bp",
            "cov_backbone", "cov_donor", "ext", "len_gain",
            "touches_left", "touches_right", "structural_score"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for c in candidates:
            w.writerow([
                c["candidate_rank"], c["backbone"], c["donor"],
                f"{c['ident']:.4f}", c["aligned_bp"],
                f"{c['cov_backbone']:.4f}", f"{c['cov_donor']:.4f}",
                c["ext"], c["len_gain"],
                c["touches_left"], c["touches_right"],
                f"{c['structural_score']:.4f}",
            ])


def _run_busco_trial(trial_fa, lineage, threads, trial_label, out_dir):
    """Run BUSCO on a trial assembly and return parsed metrics dict or None."""
    busco_bin = shutil.which("busco")
    if not busco_bin:
        return None

    trial_out = os.path.join(out_dir, trial_label)
    if os.path.isdir(trial_out):
        shutil.rmtree(trial_out)

    cmd = (f"busco -o {trial_label} -i {trial_fa} -l {lineage} "
           f"-m genome -c {threads} --out_path {out_dir} --offline")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Parse BUSCO summary
    for root_dir in [trial_out, os.path.join(out_dir, f"run_{trial_label}")]:
        if not os.path.isdir(root_dir):
            continue
        # Try JSON
        for jf in sorted(glob.glob(os.path.join(root_dir, "**", "short_summary*.json"),
                                    recursive=True), key=os.path.getmtime, reverse=True):
            try:
                with open(jf) as f:
                    data = json.load(f)
                per = data.get("percentages") or {}
                if per:
                    return {k: float(per.get(k, 0)) for k in "CSDFM"}
            except Exception:
                continue
        # Try TXT
        for tf in sorted(glob.glob(os.path.join(root_dir, "**", "short_summary*.txt"),
                                    recursive=True), key=os.path.getmtime, reverse=True):
            try:
                txt = open(tf, errors="ignore").read()
                m = re.search(r"C:(\d+(?:\.\d+)?)%.*?S:(\d+(?:\.\d+)?)%.*?"
                              r"D:(\d+(?:\.\d+)?)%.*?F:(\d+(?:\.\d+)?)%.*?"
                              r"M:(\d+(?:\.\d+)?)%", txt, re.S)
                if m:
                    return {"C": float(m.group(1)), "S": float(m.group(2)),
                            "D": float(m.group(3)), "F": float(m.group(4)),
                            "M": float(m.group(5))}
            except Exception:
                continue
    return None


def _run_purge_dups(runner, input_fa, output_fa):
    """Run purge_dups on the assembly to remove haplotigs and overlaps.

    purge_dups is designed for haploid/primary assembly deduplication.
    For polyploid genomes (e.g. plants with --taxon plant), exercise caution
    as purge_dups may incorrectly collapse homeologous sequences.

    Returns True if purge_dups ran successfully, False otherwise.
    """
    # Check all required tools
    for tool in ["purge_dups", "pbcstat", "calcuts", "split_fa", "get_seqs"]:
        if not shutil.which(tool):
            runner.log_warn(f"{tool} not found; skipping purge_dups pipeline")
            shutil.copy(input_fa, output_fa)
            return False

    # Taxon-aware warning for polyploid-prone taxa
    taxon = getattr(runner, 'taxon', 'other')
    if taxon == "plant":
        runner.log_warn("purge_dups may collapse homeologous sequences in polyploid plants. "
                        "Use --no-purge-dups if this is a polyploid genome.")

    runner.log_info("Running purge_dups for haplotig/overlap cleanup")
    pd_dir = "assemblies/purge_dups_work"
    if os.path.isdir(pd_dir):
        shutil.rmtree(pd_dir)
    os.makedirs(pd_dir, exist_ok=True)

    # Step 1: Align reads to assembly with platform-appropriate preset
    mm2_preset = "map-hifi"
    if runner.platform == "nanopore":
        mm2_preset = "map-ont"
    elif runner.platform == "pacbio":
        mm2_preset = "map-pb"

    paf = os.path.join(pd_dir, "reads_vs_asm.paf.gz")
    cmd = (f"minimap2 -x {mm2_preset} -t {runner.threads} {input_fa} {runner.fastq} "
           f"| gzip -c > {paf}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if not os.path.isfile(paf) or os.path.getsize(paf) == 0:
        runner.log_warn("purge_dups: read alignment produced no output; skipping")
        shutil.copy(input_fa, output_fa)
        return False

    # Step 2: Calculate coverage stats
    stat_file = os.path.join(pd_dir, "PB.stat")
    cutoff_file = os.path.join(pd_dir, "cutoffs")
    cmd = f"pbcstat {paf} -O {pd_dir}"
    subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if not os.path.isfile(stat_file) or os.path.getsize(stat_file) == 0:
        runner.log_warn("purge_dups pbcstat failed; skipping")
        shutil.copy(input_fa, output_fa)
        return False

    cmd = f"calcuts {stat_file} > {cutoff_file}"
    subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if not os.path.isfile(cutoff_file) or os.path.getsize(cutoff_file) == 0:
        runner.log_warn("purge_dups calcuts failed; skipping")
        shutil.copy(input_fa, output_fa)
        return False

    # Step 3: Self-alignment for duplicate detection
    split_fa = os.path.join(pd_dir, "asm.split.fa")
    cmd = f"split_fa {input_fa} > {split_fa}"
    subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if not os.path.isfile(split_fa) or os.path.getsize(split_fa) == 0:
        runner.log_warn("purge_dups split_fa failed; skipping")
        shutil.copy(input_fa, output_fa)
        return False

    self_paf = os.path.join(pd_dir, "self_aln.paf.gz")
    cmd = (f"minimap2 -x asm5 -t {runner.threads} -DP {split_fa} {split_fa} "
           f"| gzip -c > {self_paf}")
    subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Step 4: Purge duplicates
    # Note: -2 flag enables two-round purging (more aggressive).
    # Use single-round for fungi/small genomes; two-round for vertebrates/plants.
    base_cov = os.path.join(pd_dir, "PB.base.cov")
    dups_bed = os.path.join(pd_dir, "dups.bed")
    purge_flag = ""
    if taxon in ("vertebrate", "animal", "plant"):
        purge_flag = "-2"  # two-round for larger, more complex genomes
    cmd = f"purge_dups {purge_flag} -c {base_cov} -T {cutoff_file} {self_paf} > {dups_bed}"
    subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if not os.path.isfile(dups_bed) or os.path.getsize(dups_bed) == 0:
        runner.log_warn("purge_dups produced no dups.bed; skipping")
        shutil.copy(input_fa, output_fa)
        return False

    # Step 5: Get purged sequences
    # get_seqs writes purged.fa and hap.fa in its working directory.
    # Use absolute path for input FASTA since we set cwd=pd_dir.
    abs_input = os.path.abspath(input_fa)
    cmd = f"get_seqs -e {dups_bed} {abs_input}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                            cwd=pd_dir)

    purged_fa = os.path.join(pd_dir, "purged.fa")
    hap_fa = os.path.join(pd_dir, "hap.fa")

    if os.path.isfile(purged_fa) and os.path.getsize(purged_fa) > 0:
        before_recs = list(_read_fasta_records(input_fa))
        after_recs = list(_read_fasta_records(purged_fa))
        before_bp = sum(len(s) for _, s in before_recs)
        after_bp = sum(len(s) for _, s in after_recs)
        hap_n = len(list(_read_fasta_records(hap_fa))) if os.path.isfile(hap_fa) else 0
        runner.log(f"purge_dups: {len(before_recs)} → {len(after_recs)} contigs "
                   f"({hap_n} haplotigs removed), "
                   f"{before_bp:,} → {after_bp:,} bp")
        shutil.copy(purged_fa, output_fa)
        return True
    else:
        runner.log_warn("purge_dups produced no purged.fa; keeping original assembly")
        shutil.copy(input_fa, output_fa)
        return False


def _run_polishing(runner, input_fa, output_fa):
    """Run platform-appropriate polishing on the assembly.

    Platform selection logic:
      - HiFi: Skip polishing by default (HiFi reads are ~Q40+, polishing has
        minimal benefit and can introduce errors). If NextPolish2 is installed,
        run it as an optional refinement.
      - ONT:  Medaka (preferred, neural-network-based) → Racon (fallback)
      - CLR:  Racon (standard for CLR polishing)

    Returns True if polishing ran successfully.
    """
    platform = runner.platform or "pacbio-hifi"

    if platform == "pacbio-hifi":
        # HiFi assemblies are already high-accuracy (~Q40+).
        # Polishing is optional and may not improve quality.
        np2_bin = shutil.which("nextPolish2") or shutil.which("nextpolish2")
        if np2_bin:
            runner.log_info("Polishing HiFi assembly with NextPolish2 (optional refinement)")
            polish_dir = "assemblies/polish_work"
            os.makedirs(polish_dir, exist_ok=True)

            bam = os.path.join(polish_dir, "reads.sorted.bam")
            cmd = (f"minimap2 -ax map-hifi -t {runner.threads} {input_fa} {runner.fastq} "
                   f"| samtools sort -@ {runner.threads} -o {bam}")
            subprocess.run(cmd, shell=True, capture_output=True, text=True)
            subprocess.run(f"samtools index {bam}", shell=True,
                           capture_output=True, text=True)

            cmd = f"{np2_bin} {input_fa} {bam} > {output_fa}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if os.path.isfile(output_fa) and os.path.getsize(output_fa) > 0:
                runner.log("NextPolish2 polishing complete")
                return True
            else:
                runner.log_warn("NextPolish2 produced no output")
        else:
            runner.log_info("HiFi assembly: polishing skipped (HiFi reads are already "
                            "high-accuracy; install NextPolish2 for optional refinement)")
        shutil.copy(input_fa, output_fa)
        return False

    if platform == "nanopore":
        # Medaka is preferred for ONT — neural-network-based, better accuracy
        medaka_bin = shutil.which("medaka_consensus") or shutil.which("medaka")
        if medaka_bin:
            runner.log_info("Polishing ONT assembly with Medaka")
            polish_dir = "assemblies/polish_work"
            os.makedirs(polish_dir, exist_ok=True)
            medaka_out = os.path.join(polish_dir, "medaka_out")

            # Medaka model can be overridden via MEDAKA_MODEL env var.
            # Default: r1041_e82_400bps_sup_v4.3.0 (R10.4.1 chemistry, SUP basecalling)
            # For older R9.4.1 data, set MEDAKA_MODEL=r941_min_sup_g507
            medaka_model = os.environ.get("MEDAKA_MODEL", "r1041_e82_400bps_sup_v4.3.0")
            cmd = (f"medaka_consensus -i {runner.fastq} -d {input_fa} "
                   f"-o {medaka_out} -t {runner.threads} -m {medaka_model}")
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            medaka_fa = os.path.join(medaka_out, "consensus.fasta")
            if os.path.isfile(medaka_fa) and os.path.getsize(medaka_fa) > 0:
                shutil.copy(medaka_fa, output_fa)
                runner.log("Medaka polishing complete")
                return True
            else:
                runner.log_warn("Medaka produced no output; falling back to Racon")

    # Racon path: default for CLR, fallback for ONT
    racon_bin = shutil.which("racon")
    if not racon_bin:
        runner.log_warn("racon not found; skipping polishing")
        shutil.copy(input_fa, output_fa)
        return False

    runner.log_info(f"Polishing with Racon ({platform} reads)")
    mm2_preset = "map-ont" if platform == "nanopore" else "map-pb"

    polish_dir = "assemblies/polish_work"
    os.makedirs(polish_dir, exist_ok=True)
    overlap_paf = os.path.join(polish_dir, "reads_vs_asm.paf")
    cmd = (f"minimap2 -x {mm2_preset} -t {runner.threads} {input_fa} {runner.fastq} "
           f"> {overlap_paf}")
    subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Racon outputs polished FASTA to stdout
    cmd = f"racon -t {runner.threads} {runner.fastq} {overlap_paf} {input_fa}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.stdout.strip():
        with open(output_fa, "w") as f:
            f.write(result.stdout)
        if os.path.getsize(output_fa) > 0:
            runner.log("Racon polishing complete")
            return True

    runner.log_warn("Racon produced no output; keeping unpolished assembly")
    shutil.copy(input_fa, output_fa)
    return False


def step_12_refine(runner):
    """Step 12 - T2T-first telomere-aware backbone refinement with BUSCO trial validation.

    v1.2.0 workflow (T2T-first philosophy):
      12A  Merqury pre-selection (optional)
      12B  Auto-select backbone assembler
      12C  Prepare cleaned backbone + chimera safety check
      12D  T2T-first foundation building:
           12D1  Strict dedup — remove backbone contigs 95%/95% redundant to T2T pool
           12D2  Fragment removal — remove backbone fragments 50%/90% to T2T pool
           12D3  Classify backbone contigs by telomere status
           12D4  Aggressive non-telomeric dedup — remove non-telo backbone contigs
                 70%/85% redundant to T2T pool (more aggressive for contigs
                 that lack telomere support)
           12D5  Self-dedup non-telomeric backbone — remove smaller of overlapping
                 non-telo backbone contigs (80%/90%)
      12E  Telomere rescue with donor verification:
           only accept donors with verified telomere signal
      12F  BUSCO trial validation for plausible candidates only
      12G  Final combine (T2T foundation + telomere-rescued backbone gap-fill)
      12H  purge_dups default cleanup
      12I  Platform-aware polishing (Medaka/NextPolish2/Racon)
      12J  Genome-size-aware pruning (never prunes telomere-bearing contigs)
    """
    runner.log("Step 12 - Telomere-aware backbone refinement (v1.2.0)")
    os.makedirs("merqury", exist_ok=True)
    os.makedirs("assemblies", exist_ok=True)

    # ---- 12A. Optional Merqury pre-selection ----
    if runner.merqury_enable:
        _run_merqury_preselection(runner)
    else:
        runner.log_info("Merqury disabled for this run")

    _write_merqury_csv()
    _build_assembly_info(runner)

    # ---- 12B. Auto-select backbone assembler ----
    assembler = runner.assembler
    if not assembler:
        runner.log_info("Selection criteria: BUSCO + telomere + Merqury + contiguity + N50")
        assembler = _auto_select_backbone(runner)
        if assembler:
            runner.log_info(f"Auto-selected assembler: {assembler}")
        else:
            runner.log_warn("Auto-selection failed; using first available assembler as fallback")

    if not assembler:
        for asm in ["canu", "flye", "nextDenovo", "peregrine", "ipa", "hifiasm", "reference"]:
            if os.path.isfile(f"assemblies/{asm}.result.fasta") and \
               os.path.getsize(f"assemblies/{asm}.result.fasta") > 0:
                assembler = asm
                break

    if not assembler:
        runner.log_error("No valid assembler selected")
        raise RuntimeError("No assembler selected")

    asm_fa = f"assemblies/{assembler}.result.fasta"
    if not os.path.isfile(asm_fa) or os.path.getsize(asm_fa) == 0:
        runner.log_error(f"Selected assembler '{assembler}' has no valid FASTA at {asm_fa}")
        raise RuntimeError("Invalid assembler selection")

    # Read protected telomere pool from Step 10
    protected_fasta = ""
    protected_mode = "none"
    if os.path.isfile("protected_telomere_contigs.fasta") and \
       os.path.getsize("protected_telomere_contigs.fasta") > 0:
        protected_fasta = "protected_telomere_contigs.fasta"
        if os.path.isfile("protected_telomere_mode.txt"):
            with open("protected_telomere_mode.txt") as f:
                protected_mode = f.read().strip()

    runner.log_info(f"Protected mode: {protected_mode}")
    runner.log_info(f"Selected backbone: {asm_fa}")

    # ---- 12C. Prepare cleaned backbone + chimera safety ----
    backbone_fa = "assemblies/backbone.clean.fa"
    _clean_backbone_headers(asm_fa, backbone_fa)

    if protected_fasta and os.path.isfile(protected_fasta) and os.path.getsize(protected_fasta) > 0:
        max_individual_len = 0
        for telo_fa in glob.glob("assemblies/*.telo.fasta"):
            for _, seq in _read_fasta_records(telo_fa):
                max_individual_len = max(max_individual_len, len(seq))
        if max_individual_len == 0:
            for _, seq in _read_fasta_records(backbone_fa):
                max_individual_len = max(max_individual_len, len(seq))

        if max_individual_len > 0:
            chimera_threshold = int(max_individual_len * 1.5)
            protected_recs = list(_read_fasta_records(protected_fasta))
            n_before = len(protected_recs)
            filtered_recs = [(n, s) for n, s in protected_recs if len(s) <= chimera_threshold]
            n_removed = n_before - len(filtered_recs)
            if n_removed > 0:
                runner.log_warn(f"Chimera safety: removed {n_removed} protected contigs > "
                                f"{chimera_threshold:,} bp")
                _write_fasta(filtered_recs, protected_fasta)
            else:
                runner.log(f"Chimera safety: all {n_before} protected contigs pass "
                           f"(threshold {chimera_threshold:,} bp)")

    # ---- 12D. Strict protected-pool dedup (95%/95%) ----
    if protected_fasta and os.path.isfile(protected_fasta) and os.path.getsize(protected_fasta) > 0:
        runner.log_info(f"Using protected telomere contigs from {protected_fasta}")
        shutil.copy(protected_fasta, "assemblies/protected.telomere.fa")

        if shutil.which("minimap2"):
            runner.log_version("minimap2", "minimap2")
            paf = "assemblies/backbone_vs_protected.paf"
            cmd = f"minimap2 -x asm20 -t {runner.threads} assemblies/protected.telomere.fa {backbone_fa}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            with open(paf, "w") as f:
                f.write(result.stdout)

            n_dropped = _filter_redundant_to_protected(
                paf, backbone_fa, "assemblies/backbone.filtered.fa",
                cov_thr=float(os.environ.get("PROTECT_COV", "0.95")),
                id_thr=float(os.environ.get("PROTECT_ID", "0.95")),
            )
            runner.log(f"Strict dedup: removed {n_dropped} backbone contigs "
                       f"redundant to protected pool (95%/95%)")

            # Fragment removal pass (50%/90%)
            frag_paf = "assemblies/backbone_frag_vs_protected.paf"
            cmd = (f"minimap2 -x asm20 -t {runner.threads} "
                   f"assemblies/protected.telomere.fa assemblies/backbone.filtered.fa")
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            with open(frag_paf, "w") as f:
                f.write(result.stdout)

            n_frag = _filter_fragments_to_protected(
                frag_paf, "assemblies/backbone.filtered.fa",
                "assemblies/backbone.filtered.defrag.fa",
                frag_cov_thr=float(os.environ.get("FRAG_COV", "0.50")),
                frag_id_thr=float(os.environ.get("FRAG_ID", "0.90")),
            )
            if n_frag > 0:
                runner.log(f"Fragment removal: removed {n_frag} backbone fragments (50%/90%)")
                shutil.copy("assemblies/backbone.filtered.defrag.fa",
                            "assemblies/backbone.filtered.fa")
            else:
                runner.log("Fragment removal: no additional fragments removed")
        else:
            runner.log_warn("minimap2 not found; cannot filter redundant backbone contigs")
            shutil.copy(backbone_fa, "assemblies/backbone.filtered.fa")

        _name_dedup_fasta("assemblies/protected.telomere.fa",
                          "assemblies/backbone.filtered.fa",
                          "assemblies/backbone.filtered.nodup.fa")
    else:
        runner.log_warn("No protected telomere contigs; keeping backbone for rescue only")
        with open("assemblies/protected.telomere.fa", "w") as f:
            pass
        shutil.copy(backbone_fa, "assemblies/backbone.filtered.nodup.fa")

    # ---- 12D3. Classify remaining backbone contigs by telomere status ----
    backbone_nodup_fa = "assemblies/backbone.filtered.nodup.fa"
    runner.log_info("Classifying backbone contigs by telomere status")
    bb_t2t_ids, bb_telo_ids = _classify_contigs_telomere_status(
        backbone_nodup_fa, runner)
    runner.log(f"Backbone telomere census: {len(bb_t2t_ids)} T2T, "
               f"{len(bb_telo_ids)} any-telomere, out of "
               f"{sum(1 for _ in _read_fasta_records(backbone_nodup_fa))} total")

    # ---- 12D4. Aggressive non-telomeric dedup against T2T pool ----
    # Non-telomeric backbone contigs that overlap T2T pool at 70%/85% are
    # redundant partial copies that don't add telomere information.
    # Telomere-bearing backbone contigs are kept regardless (they contribute
    # chromosome-end information even if they overlap a T2T contig).
    if protected_fasta and os.path.isfile(protected_fasta) and \
       os.path.getsize(protected_fasta) > 0 and shutil.which("minimap2"):
        aggr_paf = "assemblies/backbone_aggr_vs_protected.paf"
        cmd = (f"minimap2 -x asm20 -t {runner.threads} "
               f"assemblies/protected.telomere.fa {backbone_nodup_fa}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        with open(aggr_paf, "w") as f:
            f.write(result.stdout)

        # Parse PAF: only drop NON-telomeric contigs that hit T2T at >= 70%/85%
        aggr_best = _parse_paf_best_hits(aggr_paf)
        aggr_cov = float(os.environ.get("AGGR_NONTELO_COV", "0.70"))
        aggr_id = float(os.environ.get("AGGR_NONTELO_ID", "0.85"))
        aggr_drop = {q for q, (cov, ident) in aggr_best.items()
                     if cov >= aggr_cov and ident >= aggr_id
                     and q not in bb_telo_ids}

        if aggr_drop:
            recs = [(n, s) for n, s in _read_fasta_records(backbone_nodup_fa)
                    if n not in aggr_drop]
            _write_fasta(recs, backbone_nodup_fa)
            runner.log(f"Aggressive non-telo dedup: removed {len(aggr_drop)} "
                       f"non-telomeric backbone contigs redundant to T2T pool "
                       f"({aggr_cov:.0%}/{aggr_id:.0%})")
        else:
            runner.log("Aggressive non-telo dedup: no additional contigs removed")

    # ---- 12D5. Self-dedup of non-telomeric backbone contigs ----
    # When two non-telomeric backbone contigs overlap at 80%/90%,
    # remove the shorter one.  Telomere-bearing contigs are never removed.
    n_self_dedup = _self_dedup_non_telomeric(
        backbone_nodup_fa, bb_telo_ids,
        "assemblies/backbone.filtered.selfdedup.fa",
        runner.threads,
        cov_thr=float(os.environ.get("SELFDEDUP_COV", "0.80")),
        id_thr=float(os.environ.get("SELFDEDUP_ID", "0.90")),
    )
    if n_self_dedup > 0:
        runner.log(f"Non-telo self-dedup: removed {n_self_dedup} redundant "
                   f"non-telomeric backbone contigs")
        shutil.copy("assemblies/backbone.filtered.selfdedup.fa",
                    backbone_nodup_fa)
    else:
        runner.log("Non-telo self-dedup: no additional contigs removed")

    # ---- 12E. Telomere rescue candidate screening ----
    # T2T-first philosophy: only accept rescue donors that bring telomere
    # evidence to the assembly.  Non-telomeric donors are rejected because
    # they don't improve chromosome-end completeness.
    single_tel_src = ""
    for candidate in ["single_tel_best_clean.fasta", "single_tel_best.fasta",
                      "telomere_supported_best.fasta", "single_tel_clean.fasta",
                      "single_tel.fasta"]:
        if os.path.isfile(candidate) and os.path.getsize(candidate) > 0:
            single_tel_src = candidate
            break

    working_backbone = backbone_nodup_fa

    if single_tel_src and shutil.which("minimap2"):
        runner.log_info(f"Screening telomere rescue candidates from {single_tel_src}")

        rescue_paf = "assemblies/single_tel_vs_backbone.paf"
        cmd = (f"minimap2 -x asm20 -t {runner.threads} "
               f"{working_backbone} {single_tel_src}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        with open(rescue_paf, "w") as f:
            f.write(result.stdout)

        donor_seqs = dict(_read_fasta_records(single_tel_src))
        backbone_seqs = dict(_read_fasta_records(working_backbone))

        all_hits = _parse_paf_rescue_hits(rescue_paf, donor_seqs, backbone_seqs)
        runner.log(f"Telomere rescue: {len(all_hits)} raw alignment hits")

        # Write debug TSV with all hits
        _write_rescue_debug_tsv(all_hits, "assemblies/single_tel.replaced.debug.tsv")

        # Structural screening
        rejected, candidates = _screen_rescue_candidates(
            all_hits,
            min_ident=float(os.environ.get("RESCUE_MIN_IDENT", "0.85")),
            min_aligned_bp=int(os.environ.get("RESCUE_MIN_ALN_BP", "8000")),
            min_cov_backbone=float(os.environ.get("RESCUE_MIN_COV_BB", "0.60")),
            min_cov_donor=float(os.environ.get("RESCUE_MIN_COV_DONOR", "0.50")),
            min_ext=int(os.environ.get("RESCUE_MIN_EXT", "1000")),
        )
        runner.log(f"Structural screening: {len(rejected)} rejected, "
                   f"{len(candidates)} plausible candidates")

        # ---- Donor telomere verification ----
        # T2T-first philosophy: only accept donors that carry telomere signal.
        # Write donor FASTA for classification, then filter candidates.
        donor_telo_verified = set()
        if candidates:
            donor_check_fa = "assemblies/rescue_donors_check.fa"
            unique_donors = {c["donor"] for c in candidates}
            _write_fasta([(n, donor_seqs[n]) for n in unique_donors if n in donor_seqs],
                         donor_check_fa)
            _, donor_telo_set = _classify_contigs_telomere_status(
                donor_check_fa, runner)
            donor_telo_verified = donor_telo_set

            # Move non-telomeric donors to rejected
            telo_verified_candidates = []
            for c in candidates:
                if c["donor"] in donor_telo_verified:
                    telo_verified_candidates.append(c)
                else:
                    c["reject_reason"] = "donor_no_telomere_signal"
                    rejected.append(c)

            n_donor_reject = len(candidates) - len(telo_verified_candidates)
            if n_donor_reject > 0:
                runner.log(f"Donor telomere verification: rejected {n_donor_reject} "
                           f"donors lacking telomere signal")
            candidates = telo_verified_candidates

            # Re-rank after filtering
            for i, c in enumerate(candidates):
                c["candidate_rank"] = i + 1

        runner.log(f"After donor verification: {len(candidates)} telomere-verified candidates")

        # Write candidate TSV
        _write_candidates_tsv(candidates, "assemblies/single_tel.candidates.tsv")

        # Write rejection summary
        with open("assemblies/rescue_rejection_summary.txt", "w") as f:
            f.write(f"Total hits evaluated: {len(all_hits)}\n")
            f.write(f"Obvious rejections: {len(rejected)}\n")
            f.write(f"Plausible candidates: {len(candidates)}\n\n")
            for r in rejected:
                f.write(f"REJECT {r['donor']} → {r['backbone']}: {r['reject_reason']}\n")

        # ---- 12F. BUSCO trial validation for plausible candidates ----
        max_accepted = int(os.environ.get("STEP12_MAX_ACCEPTED", "20"))
        min_bp_ratio = float(os.environ.get("STEP12_MIN_BP_RATIO", "0.90"))

        # Taxon-aware BUSCO thresholds:
        # - Haploid/fungal: strict (2% C-drop max)
        # - Plant: relaxed (polyploidy inflates D naturally, allow 4% C-drop)
        # - Vertebrate/animal: moderate (3% C-drop)
        taxon = getattr(runner, 'taxon', 'other')
        default_c_drop = "2.0"
        default_m_rise = "0.3"
        if taxon == "plant":
            default_c_drop = "4.0"
            default_m_rise = "1.0"
        elif taxon in ("vertebrate", "animal"):
            default_c_drop = "3.0"
            default_m_rise = "0.5"

        max_busco_c_drop = float(os.environ.get("STEP12_MAX_BUSCO_C_DROP", default_c_drop))
        max_busco_m_rise = float(os.environ.get("STEP12_MAX_BUSCO_M_RISE", default_m_rise))

        lineage = runner.busco_lineage or "ascomycota_odb10"
        busco_available = shutil.which("busco") is not None and runner.run_busco

        trial_dir = "assemblies/rescue_trials"
        os.makedirs(trial_dir, exist_ok=True)

        # Get baseline metrics
        baseline_recs = list(_read_fasta_records(working_backbone))
        baseline_bp = sum(len(s) for _, s in baseline_recs)
        baseline_contigs = len(baseline_recs)
        baseline_busco = None

        if busco_available and candidates:
            runner.log_info("Computing baseline BUSCO for trial validation")
            baseline_busco = _run_busco_trial(working_backbone, lineage,
                                               runner.threads, "baseline", trial_dir)
            if baseline_busco:
                runner.log(f"Baseline BUSCO: C={baseline_busco['C']:.1f}% "
                           f"M={baseline_busco['M']:.1f}%")
            else:
                runner.log_warn("Baseline BUSCO parse failed; using structural checks only")

        accepted_count = 0
        trial_results = []
        current_backbone = dict(backbone_seqs)  # working copy

        for cand in candidates:
            if accepted_count >= max_accepted:
                runner.log_info(f"Reached max accepted rescues ({max_accepted})")
                break

            backbone_name = cand["backbone"]
            donor_name = cand["donor"]

            if backbone_name not in current_backbone:
                continue
            if donor_name not in donor_seqs:
                continue

            # Build trial assembly: replace backbone contig with donor
            trial_recs = []
            for n, s in current_backbone.items():
                if n == backbone_name:
                    trial_recs.append((donor_name, donor_seqs[donor_name]))
                else:
                    trial_recs.append((n, s))

            trial_fa = os.path.join(trial_dir, f"trial_{cand['candidate_rank']}.fa")
            _write_fasta(trial_recs, trial_fa)

            trial_bp = sum(len(s) for _, s in trial_recs)
            trial_contigs = len(trial_recs)

            # Size check
            if baseline_bp > 0 and trial_bp < baseline_bp * min_bp_ratio:
                reason = "bp_drop"
                trial_results.append({
                    "candidate_rank": cand["candidate_rank"],
                    "backbone": backbone_name, "donor": donor_name,
                    "accepted": False, "reason": reason,
                    "trial_bp": trial_bp, "trial_contigs": trial_contigs,
                    "trial_supported": "", "baseline_busco_c": "",
                    "trial_busco_c": "", "baseline_busco_m": "",
                    "trial_busco_m": "", "busco_c_delta": "", "busco_m_delta": "",
                })
                runner.log(f"  Candidate {cand['candidate_rank']} "
                           f"({donor_name} → {backbone_name}): REJECTED (bp_drop)")
                continue

            # BUSCO trial if available
            trial_busco = None
            if busco_available and baseline_busco:
                trial_label = f"trial_{cand['candidate_rank']}"
                trial_busco = _run_busco_trial(trial_fa, lineage,
                                                runner.threads, trial_label, trial_dir)

            accepted = True
            reason = "accepted"

            if trial_busco and baseline_busco:
                c_delta = trial_busco["C"] - baseline_busco["C"]
                m_delta = trial_busco["M"] - baseline_busco["M"]

                if c_delta < -max_busco_c_drop:
                    accepted = False
                    reason = "busco_c_drop_gt_2"
                elif m_delta > max_busco_m_rise:
                    accepted = False
                    reason = "busco_m_rise"
            elif busco_available and baseline_busco and trial_busco is None:
                # BUSCO ran but failed to parse
                accepted = False
                reason = "trial_busco_failed"

            trial_results.append({
                "candidate_rank": cand["candidate_rank"],
                "backbone": backbone_name, "donor": donor_name,
                "accepted": accepted, "reason": reason,
                "trial_bp": trial_bp, "trial_contigs": trial_contigs,
                "trial_supported": "",
                "baseline_busco_c": f"{baseline_busco['C']:.1f}" if baseline_busco else "",
                "trial_busco_c": f"{trial_busco['C']:.1f}" if trial_busco else "",
                "baseline_busco_m": f"{baseline_busco['M']:.1f}" if baseline_busco else "",
                "trial_busco_m": f"{trial_busco['M']:.1f}" if trial_busco else "",
                "busco_c_delta": f"{c_delta:.1f}" if (trial_busco and baseline_busco) else "",
                "busco_m_delta": f"{m_delta:.1f}" if (trial_busco and baseline_busco) else "",
            })

            if accepted:
                accepted_count += 1
                # Update working backbone
                del current_backbone[backbone_name]
                current_backbone[donor_name] = donor_seqs[donor_name]
                # Update baseline BUSCO for next iteration
                if trial_busco:
                    baseline_busco = trial_busco
                baseline_bp = trial_bp
                runner.log(f"  Candidate {cand['candidate_rank']} "
                           f"({donor_name} → {backbone_name}): ACCEPTED")
            else:
                runner.log(f"  Candidate {cand['candidate_rank']} "
                           f"({donor_name} → {backbone_name}): REJECTED ({reason})")

        runner.log(f"Telomere rescue: {accepted_count} replacements accepted")

        # Write trial summary
        trial_cols = ["candidate_rank", "backbone", "donor", "accepted", "reason",
                      "trial_bp", "trial_contigs", "trial_supported",
                      "baseline_busco_c", "trial_busco_c", "baseline_busco_m",
                      "trial_busco_m", "busco_c_delta", "busco_m_delta"]
        with open("assemblies/rescue_trial_summary.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(trial_cols)
            for tr in trial_results:
                w.writerow([tr.get(c, "") for c in trial_cols])

        # Write rescued backbone
        rescued_recs = list(current_backbone.items())
        _write_fasta(rescued_recs, "assemblies/backbone.telomere_rescued.fa")

        # Write replaced IDs for reporting
        with open("assemblies/single_tel.replaced.ids", "w") as f:
            for tr in trial_results:
                if tr["accepted"]:
                    f.write(f"{tr['backbone']}\t{tr['donor']}\n")
    else:
        if not single_tel_src:
            runner.log_info("No single-end telomeric contigs available for rescue")
        elif not shutil.which("minimap2"):
            runner.log_warn("minimap2 not found; skipping telomere rescue")

        shutil.copy(working_backbone, "assemblies/backbone.telomere_rescued.fa")

        # Write empty output files for consistency
        for p in ["assemblies/single_tel.replaced.debug.tsv",
                   "assemblies/single_tel.candidates.tsv",
                   "assemblies/rescue_rejection_summary.txt",
                   "assemblies/rescue_trial_summary.tsv"]:
            with open(p, "w") as f:
                pass
        with open("assemblies/single_tel.replaced.ids", "w") as f:
            pass

    # ---- 12F2. Post-rescue dedup against protected set ----
    if os.path.isfile("assemblies/protected.telomere.fa") and \
       os.path.getsize("assemblies/protected.telomere.fa") > 0 and \
       os.path.isfile("assemblies/backbone.telomere_rescued.fa") and \
       shutil.which("minimap2"):
        dedup_paf = "assemblies/rescued_vs_protected.paf"
        cmd = (f"minimap2 -x asm20 -t {runner.threads} "
               f"assemblies/protected.telomere.fa assemblies/backbone.telomere_rescued.fa")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        with open(dedup_paf, "w") as f:
            f.write(result.stdout)

        _filter_redundant_to_protected(
            dedup_paf, "assemblies/backbone.telomere_rescued.fa",
            "assemblies/backbone.telomere_rescued.dedup.fa",
            cov_thr=0.95, id_thr=0.95,
        )
    else:
        shutil.copy("assemblies/backbone.telomere_rescued.fa",
                     "assemblies/backbone.telomere_rescued.dedup.fa")

    # ---- 12G. Final combine ----
    prot = "assemblies/protected.telomere.fa"
    rescued = "assemblies/backbone.telomere_rescued.dedup.fa"
    raw_out = "assemblies/final_merge.raw.fasta"

    if os.path.isfile(prot) and os.path.getsize(prot) > 0 and \
       os.path.isfile(rescued) and os.path.getsize(rescued) > 0:
        with open(raw_out, "w") as out:
            with open(prot) as f:
                out.write(f.read())
            with open(rescued) as f:
                out.write(f.read())
    elif os.path.isfile(prot) and os.path.getsize(prot) > 0:
        shutil.copy(prot, raw_out)
    elif os.path.isfile(rescued) and os.path.getsize(rescued) > 0:
        shutil.copy(rescued, raw_out)
    else:
        shutil.copy(asm_fa, raw_out)

    runner.log(f"Built assemblies/final_merge.raw.fasta (mode: {protected_mode})")

    # ---- 12H. purge_dups default cleanup ----
    if not getattr(runner, 'no_purge_dups', False):
        purged_fa = "assemblies/final_merge.purged.fasta"
        _run_purge_dups(runner, raw_out, purged_fa)
        if os.path.isfile(purged_fa) and os.path.getsize(purged_fa) > 0:
            shutil.copy(purged_fa, raw_out)
    else:
        runner.log_info("purge_dups skipped (--no-purge-dups)")

    # ---- 12I. Platform-aware polishing ----
    if not getattr(runner, 'no_polish', False):
        polished_fa = "assemblies/final_merge.polished.fasta"
        _run_polishing(runner, raw_out, polished_fa)
        if os.path.isfile(polished_fa) and os.path.getsize(polished_fa) > 0:
            shutil.copy(polished_fa, raw_out)
    else:
        runner.log_info("Polishing skipped (--no-polish)")

    # ---- 12J. Genome-size-aware pruning (safety net) ----
    # T2T-first: never prune telomere-bearing contigs.  Classify the final
    # assembly and protect anything with telomere signal, not just the
    # original T2T pool.
    expected_size = _parse_genome_size(runner.genomesize)
    if expected_size > 0:
        all_recs = list(_read_fasta_records(raw_out))
        total_len = sum(len(s) for _, s in all_recs)
        budget = int(expected_size * 1.15)

        if total_len > budget:
            runner.log_warn(f"Assembly size {total_len:,} bp exceeds "
                            f"{expected_size:,} bp + 15% = {budget:,} bp")

            # Build expanded protected set: original T2T pool + any contig
            # in the final assembly that carries telomere signal
            protected_names = set()
            if os.path.isfile(prot):
                for name, _ in _read_fasta_records(prot):
                    protected_names.add(name)

            # Classify final assembly contigs for telomere status
            final_t2t, final_telo = _classify_contigs_telomere_status(
                raw_out, runner)
            protected_names |= final_telo  # protect ALL telomere-bearing contigs

            keep_always = [(n, s) for n, s in all_recs if n in protected_names]
            removable = [(n, s) for n, s in all_recs if n not in protected_names]
            removable.sort(key=lambda x: len(x[1]))

            runner.log(f"Pruning protection: {len(keep_always)} telomere-bearing "
                       f"contigs protected, {len(removable)} non-telomeric eligible")

            keep_len = sum(len(s) for _, s in keep_always)
            kept_backbone = []
            running_len = keep_len
            for name, seq in reversed(removable):
                if running_len + len(seq) <= budget:
                    kept_backbone.append((name, seq))
                    running_len += len(seq)
                else:
                    runner.log(f"  Pruning '{name}' ({len(seq):,} bp, non-telomeric)")

            n_pruned = len(removable) - len(kept_backbone)
            if n_pruned > 0:
                final_recs = keep_always + kept_backbone
                _write_fasta(final_recs, raw_out)
                new_total = sum(len(s) for _, s in final_recs)
                runner.log(f"Genome-size pruning: removed {n_pruned} non-telomeric "
                           f"fragments ({total_len:,} → {new_total:,} bp)")
            else:
                runner.log("Genome-size pruning: no contigs removed")
        else:
            runner.log(f"Assembly size {total_len:,} bp within budget "
                       f"({budget:,} bp) — no pruning needed")

    n_sorted = _fasta_sort_minlen(raw_out, f"merged_{assembler}_sort.fa",
                                   prefix="contig", minlen=500)
    runner.log(f"Sorted final merged assembly: {n_sorted} contigs >= 500 bp")

    shutil.copy(f"merged_{assembler}_sort.fa", "assemblies/final.merged.fasta")
    runner.log("Wrote assemblies/final.merged.fasta")


def step_13_busco_final(runner):
    """Step 13 - BUSCO analysis on final merged assembly."""
    runner.log("Step 13 - BUSCO analysis")
    os.makedirs("busco", exist_ok=True)
    os.makedirs("assemblies", exist_ok=True)

    final_fa = "assemblies/final.merged.fasta"
    if not os.path.isfile(final_fa) or os.path.getsize(final_fa) == 0:
        runner.log_error(f"Final merged FASTA '{final_fa}' not found.")
        raise RuntimeError("No final assembly")

    lineage = runner.busco_lineage or "ascomycota_odb10"

    def has_final_metrics():
        for d in ["busco/final", "busco/run_final"]:
            if os.path.isdir(d):
                files = glob.glob(f"{d}/**/short_summary*.json", recursive=True) + \
                        glob.glob(f"{d}/**/short_summary*.txt", recursive=True) + \
                        glob.glob(f"{d}/**/full_table*.tsv", recursive=True)
                if files:
                    return True
        return False

    if has_final_metrics():
        runner.log("Found existing BUSCO metrics for final")
    else:
        busco_bin = shutil.which("busco")
        if not busco_bin:
            runner.log_error("Missing final BUSCO metrics and 'busco' is not available")
            raise RuntimeError("BUSCO not available")

        runner.log_info(f"BUSCO on final (lineage={lineage}, threads={runner.threads})")
        runner.log_version("busco", "busco")

        cwd = os.getcwd()
        os.chdir("busco")
        cmd = f"busco -o final -i ../{final_fa} -l {lineage} -m genome -c {runner.threads}"
        try:
            runner.run_cmd(cmd, desc="Running BUSCO on final", check=False)
        finally:
            os.chdir(cwd)

    # Parse BUSCO results from busco/final or busco/run_final
    def _parse_busco_final():
        """Parse BUSCO metrics from final run outputs."""
        roots = []
        for d in ("busco/final", "busco/run_final"):
            if os.path.isdir(d):
                roots.append(d)
        if not roots:
            return None

        for root in sorted(roots, key=os.path.getmtime, reverse=True):
            # Try JSON first
            js_files = sorted(glob.glob(os.path.join(root, "**", "short_summary*.json"), recursive=True),
                              key=os.path.getmtime, reverse=True)
            for jf in js_files:
                try:
                    with open(jf) as f:
                        data = json.load(f)
                    n = data.get("n") or (data.get("lineage_dataset") or {}).get("n")
                    per = data.get("percentages") or {}
                    if n and per:
                        return {k: float(per.get(k, 0)) for k in "CSDFM"} | {"n": int(n)}
                except Exception:
                    continue

            # Try TXT
            txt_files = sorted(glob.glob(os.path.join(root, "**", "short_summary*.txt"), recursive=True),
                               key=os.path.getmtime, reverse=True)
            for tf in txt_files:
                try:
                    txt = open(tf, errors="ignore").read()
                    m = re.search(r"C:(\d+(?:\.\d+)?)%.*?S:(\d+(?:\.\d+)?)%.*?D:(\d+(?:\.\d+)?)%.*?F:(\d+(?:\.\d+)?)%.*?M:(\d+(?:\.\d+)?)%.*?n:(\d+)", txt, re.S)
                    if m:
                        return {"C": float(m.group(1)), "S": float(m.group(2)),
                                "D": float(m.group(3)), "F": float(m.group(4)),
                                "M": float(m.group(5)), "n": int(m.group(6))}
                except Exception:
                    continue

            # Try full_table TSV
            ft_files = sorted(glob.glob(os.path.join(root, "**", "full_table*.tsv"), recursive=True),
                              key=os.path.getmtime, reverse=True)
            for ft in ft_files:
                try:
                    S = D = F = M = n = 0
                    for ln in open(ft):
                        if not ln or ln.startswith("#"):
                            continue
                        parts = ln.split("\t")
                        st = parts[1].strip() if len(parts) > 1 else ""
                        n += 1
                        if st == "Complete": S += 1
                        elif st == "Duplicated": D += 1
                        elif st == "Fragmented": F += 1
                        elif st == "Missing": M += 1
                    if n > 0:
                        pct = lambda x: 100.0 * float(x) / float(n) if n else 0.0
                        return {"C": pct(S + D), "S": pct(S), "D": pct(D),
                                "F": pct(F), "M": pct(M), "n": n}
                except Exception:
                    continue
        return None

    metrics = _parse_busco_final()
    def fmt_pct(x):
        s = f"{float(x):.1f}"
        return s.rstrip('0').rstrip('.') if '.' in s else s

    if metrics:
        rows = [
            ["Metric", "merged"],
            ["BUSCO C (%)", fmt_pct(metrics["C"])],
            ["BUSCO S (%)", fmt_pct(metrics["S"])],
            ["BUSCO D (%)", fmt_pct(metrics["D"])],
            ["BUSCO F (%)", fmt_pct(metrics["F"])],
            ["BUSCO M (%)", fmt_pct(metrics["M"])],
        ]
    else:
        runner.log_warn("Could not parse BUSCO outputs for final assembly")
        rows = [
            ["Metric", "merged"],
            ["BUSCO C (%)", ""],
            ["BUSCO S (%)", ""],
            ["BUSCO D (%)", ""],
            ["BUSCO F (%)", ""],
            ["BUSCO M (%)", ""],
        ]

    with open("assemblies/merged.busco.csv", "w", newline="") as f:
        csv.writer(f).writerows(rows)
    runner.log("Wrote assemblies/merged.busco.csv")


def step_14_telomere_final(runner):
    """Step 14 - Telomere analysis on final assembly (hybrid detection)."""
    runner.log("Step 14 - Telomere analysis (final assembly, hybrid detection)")
    os.makedirs("assemblies", exist_ok=True)

    final_fa = "assemblies/final.merged.fasta"
    if not os.path.isfile(final_fa) or os.path.getsize(final_fa) == 0:
        runner.log_error(f"Final merged FASTA '{final_fa}' not found.")
        raise RuntimeError("No final assembly")

    # Run hybrid telomere detection on final assembly
    runner.log(f"Running hybrid telomere detection on final assembly (mode: {runner.telomere_mode})")
    try:
        results = detect_telomeres(
            final_fa,
            mode=runner.telomere_mode,
            user_motif=runner.motif,
            end_window=runner.telo_end_window,
            score_window=runner.telo_score_window,
            kmer_min=runner.telo_kmer_min,
            kmer_max=runner.telo_kmer_max,
            threads=runner.threads,
            taxon=getattr(runner, 'taxon', 'other'),
        )
        write_detection_outputs(results, final_fa, "assemblies/final")
    except Exception as e:
        runner.log_warn(f"Hybrid detection on final assembly failed: {e}")
        results = []

    # Count tiers
    counts = defaultdict(int)
    for r in results:
        counts[r["classification"]] += 1
    t2t = counts.get("strict_t2t", 0)
    single = counts.get("single_tel_strong", 0)
    supported = counts.get("telomere_supported", 0)

    runner.log(f"Final assembly telomeres: {t2t} strict_t2t, {single} single_tel_strong, "
               f"{supported} telomere_supported")

    # Write merged.telo.csv with tier-based counts
    rows = [
        ["Metric", "merged"],
        ["Telomere strict T2T contigs", str(t2t)],
        ["Telomere single-end strong contigs", str(single)],
        ["Telomere-supported contigs", str(t2t + single + supported)],
    ]
    with open("assemblies/merged.telo.csv", "w", newline="") as f:
        csv.writer(f).writerows(rows)
    runner.log("Wrote assemblies/merged.telo.csv")


def step_15_quast_final(runner):
    """Step 15 - QUAST metrics for final assembly."""
    runner.log("Step 15 - QUAST metrics for final assembly")
    os.makedirs("quast_final", exist_ok=True)
    os.makedirs("assemblies", exist_ok=True)

    final_fa = "assemblies/final.merged.fasta"
    if not os.path.isfile(final_fa) or os.path.getsize(final_fa) == 0:
        runner.log_error(f"Final merged FASTA '{final_fa}' not found.")
        raise RuntimeError("No final assembly")

    quast_bin = shutil.which("quast.py") or shutil.which("quast")
    if quast_bin:
        cmd = f"{quast_bin} {final_fa} --labels final --threads {runner.threads} -o quast_final"
        runner.run_cmd(cmd, desc="Running QUAST on final assembly", check=False)

    # Parse QUAST output into merged.quast.csv
    quast_report = os.path.join("quast_final", "report.tsv")
    quast_treport = os.path.join("quast_final", "transposed_report.tsv")
    rows = [["Metric", "merged"]]

    if os.path.exists(quast_report):
        with open(quast_report) as f:
            for line in f:
                parts = line.rstrip('\r\n').split('\t')
                if len(parts) >= 2:
                    rows.append([parts[0], parts[1]])
    elif os.path.exists(quast_treport):
        with open(quast_treport) as f:
            treport_rows = list(csv.reader(f, delimiter='\t'))
        if len(treport_rows) >= 2:
            metrics = treport_rows[0][1:]
            vals = treport_rows[1][1:] if len(treport_rows) > 1 else []
            for i, m in enumerate(metrics):
                v = vals[i] if i < len(vals) else ""
                rows.append([m, v])

    with open("assemblies/merged.quast.csv", "w", newline="") as f:
        csv.writer(f).writerows(rows)
    runner.log("Wrote assemblies/merged.quast.csv")


def step_16_final_report(runner):
    """Step 16 - Generate final comparison report.

    Runs Merqury on final assembly (if enabled), then combines per-assembler
    metrics (assembly_info.csv) with final-merged metrics into
    final_results/final_result.csv.
    """
    runner.log("Step 16 - Final assembly comparison")
    os.makedirs("final_results", exist_ok=True)
    os.makedirs("merqury", exist_ok=True)

    # ---- Run Merqury on final assembly ----
    if runner.merqury_enable:
        db = runner.merqury_db
        if not db:
            for cand in ["reads.meryl", "meryl/reads.meryl", "merqury/reads.meryl"] + glob.glob("*.meryl"):
                if os.path.isdir(cand):
                    db = cand
                    break
        if db and os.path.isdir(db) and shutil.which("merqury.sh"):
            runner.log_info("Running Merqury on final assembly")
            runner.log_version("merqury.sh", "merqury.sh")
            cmd = f"merqury.sh {db} assemblies/final.merged.fasta merqury/final"
            runner.run_cmd(cmd, desc="Merqury on final assembly", check=False)
        else:
            runner.log_warn("Merqury requested for final assembly but merqury.sh or .meryl database was not found; skipping.")
    else:
        runner.log_info("Merqury disabled for final assembly comparison")

    # Write merged Merqury CSV
    def _parse_first_float(path):
        if not os.path.exists(path):
            return ""
        txt = open(path, "r", errors="ignore").read()
        for pat in [r'(?i)qv[^0-9]*([0-9]+(?:\.[0-9]+)?)',
                    r'(?i)completeness[^0-9]*([0-9]+(?:\.[0-9]+)?)',
                    r'([0-9]+(?:\.[0-9]+)?)']:
            m = re.search(pat, txt)
            if m:
                return m.group(1)
        return ""

    qv = _parse_first_float("merqury/final.qv")
    comp = _parse_first_float("merqury/final.completeness.stats")
    with open("assemblies/merged.merqury.csv", "w", newline="") as f:
        csv.writer(f).writerows([
            ["Metric", "merged"],
            ["Merqury QV", qv],
            ["Merqury completeness (%)", comp],
        ])

    # ---- Build comprehensive final report ----
    # Load all merged metric files
    from collections import OrderedDict
    final_map = OrderedDict()

    for _, path in [("telo", "assemblies/merged.telo.csv"),
                    ("busco", "assemblies/merged.busco.csv"),
                    ("quast", "assemblies/merged.quast.csv"),
                    ("merqury", "assemblies/merged.merqury.csv")]:
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            continue
        with open(path, newline="") as f:
            r = csv.reader(f)
            hdr = next(r, None)
            if not hdr or len(hdr) < 2:
                continue
            idx_metric, idx_merged = 0, 1
            for i, h in enumerate(hdr):
                name = (h or "").strip().lower()
                if name == "metric":
                    idx_metric = i
                if name == "merged":
                    idx_merged = i
            for row in r:
                if not row:
                    continue
                m_name = (row[idx_metric] if idx_metric < len(row) else "").strip()
                v = (row[idx_merged] if idx_merged < len(row) else "").strip()
                if m_name:
                    final_map[m_name] = v

    # Add telomere tier counts from score TSV if available
    _telo_scores = "assemblies/final.telomere_end_scores.tsv"
    if os.path.isfile(_telo_scores):
        _t2t = _single = _supported = 0
        with open(_telo_scores) as tf:
            for tr in csv.DictReader(tf, delimiter="\t"):
                _tier = tr.get("tier", "none")
                if _tier == "strict_t2t":
                    _t2t += 1
                elif _tier == "single_tel_strong":
                    _single += 1
                elif _tier == "telomere_supported":
                    _supported += 1
        final_map["Telomere strict T2T contigs"] = str(_t2t)
        final_map["Telomere single-end strong contigs"] = str(_single)
        final_map["Telomere-supported contigs"] = str(_t2t + _single + _supported)

    # Add selection metadata
    decision_file = "assemblies/selection_decision.txt"
    if os.path.isfile(decision_file):
        with open(decision_file) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                k, *rest = line.split("\t", 1)
                v = rest[0] if rest else ""
                if k == "selected_score":
                    final_map["Selection score"] = v
                elif k == "selected_assembler":
                    final_map["Selected assembler"] = v
                elif k == "auto_mode":
                    final_map["Auto-selection mode"] = v
                elif k == "score_formula":
                    final_map["Score formula"] = v

    # Protected mode
    if os.path.isfile("protected_telomere_mode.txt"):
        with open("protected_telomere_mode.txt") as f:
            final_map["Protected telomere mode"] = f.read().strip()

    # Pool stats
    def _count_fasta(path):
        try:
            with open(path) as f:
                return str(sum(1 for ln in f if ln.startswith(">")))
        except Exception:
            return ""

    final_map["Step10 strict T2T pool contigs"] = _count_fasta("t2t_clean.fasta")
    final_map["Step10 best single-end telomere pool contigs"] = _count_fasta("single_tel_best_clean.fasta")
    final_map["Step10 optimized telomere-supported pool contigs"] = _count_fasta("telomere_supported_best_clean.fasta")

    # Rescue counts
    try:
        with open("assemblies/single_tel.replaced.ids") as f:
            final_map["Step12 rescued telomere replacements"] = str(sum(1 for _ in f))
    except Exception:
        final_map["Step12 rescued telomere replacements"] = ""

    # ---- Combine assembly_info.csv with merged metrics ----
    info_csv = "assemblies/assembly_info.csv"
    if os.path.isfile(info_csv) and os.path.getsize(info_csv) > 0:
        with open(info_csv, newline="") as f:
            r = csv.reader(f)
            hdr = next(r, None)
            body = [row for row in r if row]

        if hdr:
            hdr_out = list(hdr)
            if "merged" not in [h.lower() for h in hdr]:
                hdr_out.append("merged")

            out_rows = []
            for row in body:
                base_row = row[:]
                while len(base_row) < len(hdr_out) - 1:
                    base_row.append("")
                metric = (base_row[0] if base_row else "").strip()
                out_rows.append(base_row + [final_map.get(metric, "")])

            present_metrics = set(r[0].strip() for r in body if r and r[0].strip())
            for m, v in final_map.items():
                if m not in present_metrics:
                    r = [""] * len(hdr_out)
                    r[0] = m
                    r[-1] = v
                    out_rows.append(r)

            with open("final_results/final_result.csv", "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(hdr_out)
                w.writerows(out_rows)
        else:
            with open("final_results/final_result.csv", "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["Metric", "merged"])
                for m, v in final_map.items():
                    w.writerow([m, v])
    else:
        with open("final_results/final_result.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Metric", "merged"])
            for m, v in final_map.items():
                w.writerow([m, v])

    runner.log("Wrote final_results/final_result.csv")

    # Also copy assembly_info for easy access
    if os.path.isfile("assemblies/assembly_info.csv"):
        shutil.copy("assemblies/assembly_info.csv", "final_results/assembly_info.csv")


def step_17_cleanup(runner):
    """Step 17 - Cleanup temporary files into structured output folders."""
    runner.log("Step 17 - Cleanup temporary files")
    os.makedirs("temp/merge/fasta", exist_ok=True)
    os.makedirs("temp/merge/param", exist_ok=True)
    os.makedirs("temp/busco", exist_ok=True)
    os.makedirs("temp/log", exist_ok=True)
    os.makedirs("final_results", exist_ok=True)

    def _safe_move(src, dst_dir):
        """Move src into dst_dir, overwriting any existing file of the same name."""
        dst = os.path.join(dst_dir, os.path.basename(src))
        if os.path.exists(dst):
            os.remove(dst)
        shutil.move(src, dst_dir)

    # Move merge intermediates
    for pattern in ["aln_summary_merged*.tsv", "anchor_summary_merged_*.txt"]:
        for f in glob.glob(pattern):
            try:
                _safe_move(f, "temp/merge/")
            except Exception:
                pass
    for pattern in [".merged_*", "merged_*.delta", "merged_*.coords", "merged_*.snps",
                    "merged_*.delta.*", "merged_*.crunch", "merged_*.filter",
                    "merged_*.qdiff", "merged_*.rdiff", "merged_*.mcoords"]:
        for f in glob.glob(pattern):
            try:
                _safe_move(f, "temp/merge/")
            except Exception:
                pass

    # Move BUSCO logs
    if os.path.isdir("busco"):
        for f in glob.glob("busco/**/*.log", recursive=True):
            try:
                _safe_move(f, "temp/busco/")
            except Exception:
                pass
    for f in glob.glob("busco_*.log") + glob.glob("*busco*.log"):
        try:
            _safe_move(f, "temp/busco/")
        except Exception:
            pass

    # Move merged FASTA intermediates
    for f in glob.glob("merged_*.fasta") + glob.glob("merged_*.fa"):
        try:
            _safe_move(f, "temp/merge/fasta/")
        except Exception:
            pass

    # Move param files
    for f in glob.glob("param_summary_merged_*.txt"):
        try:
            _safe_move(f, "temp/merge/param/")
        except Exception:
            pass

    # Move misc logs
    for f in glob.glob("*.log"):
        if "busco" in f.lower():
            continue
        try:
            shutil.move(f, "temp/log/")
        except Exception:
            pass

    # Move key results to final_results/
    for src in [
        "assemblies/final.merged.fasta",
        "assemblies/selection_debug.tsv",
        "assemblies/selection_decision.txt",
        "assemblies/merged.merqury.csv",
        "assemblies/merged.busco.csv",
        "assemblies/merged.quast.csv",
        "assemblies/merged.telo.csv",
        "telomere_cluster_summary.tsv",
        "telomere_support_summary.csv",
        "protected_telomere_mode.txt",
    ]:
        if os.path.isfile(src):
            try:
                shutil.move(src, "final_results/")
            except Exception:
                pass

    runner.log("Cleanup complete.")


def step_18_assembly_only(runner):
    """Step 18 - Assembly-only comparison summary (with optional Merqury)."""
    runner.log("Step 18 - Assembly-only comparison summary")
    os.makedirs("assemblies", exist_ok=True)
    os.makedirs("merqury", exist_ok=True)
    os.makedirs("final_results", exist_ok=True)

    # Run Merqury on all assemblies if enabled
    if runner.merqury_enable:
        _run_merqury_preselection(runner)
    else:
        runner.log_info("Merqury disabled for assembly-only comparison")

    # Write Merqury CSV
    _write_merqury_csv()

    # Create assembly info (combines BUSCO + QUAST + telomere + Merqury)
    _build_assembly_info(runner)

    if os.path.isfile("assemblies/assembly_info.csv"):
        shutil.copy("assemblies/assembly_info.csv", "final_results/assembly_only_result.csv")
        runner.log("Wrote final_results/assembly_only_result.csv")

    runner.log("Wrote assemblies/assembly_info.csv")


STEP_FUNCTIONS = {
    1: step_01_canu,
    2: step_02_nextdenovo,
    3: step_03_peregrine,
    4: step_04_ipa,
    5: step_05_flye,
    6: step_06_hifiasm,
    7: step_07_normalize,
    8: step_08_busco,
    9: step_09_telomere,
    10: step_10_telomere_pool,
    11: step_11_quast,
    12: step_12_refine,
    13: step_13_busco_final,
    14: step_14_telomere_final,
    15: step_15_quast_final,
    16: step_16_final_report,
    17: step_17_cleanup,
    18: step_18_assembly_only,
}
