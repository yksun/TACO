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

    flag_map = {"pacbio-hifi": "--pacbio-hifi", "nanopore": "--nano-hq", "pacbio": "--pacbio-raw"}
    flye_flag = flag_map.get(runner.platform, "--pacbio-hifi")

    runner.log_version("flye", "flye")
    cmd = f"flye {flye_flag} {runner.fastq} --out-dir flye --threads {runner.threads}"
    result = runner.run_cmd(cmd, desc="Running flye", check=False)
    if result.returncode != 0:
        runner.log_warn("Step 5: flye failed. Skipping. Check logs/step_5.log for details.")


def step_06_hifiasm(runner):
    """Step 6 - Assembly using Hifiasm (non-fatal)."""
    runner.log("Step 6 - Assembly of the genome using Hifiasm")

    if not shutil.which("hifiasm"):
        _assembler_skip(runner, 6, "hifiasm", "binary not found. Install via: conda install -c bioconda hifiasm")
        return

    os.makedirs("hifiasm", exist_ok=True)
    runner.log_version("hifiasm", "hifiasm")

    # Platform-specific flags for hifiasm
    platform_flag = ""
    if runner.platform == "nanopore":
        platform_flag = "--ont"
    elif runner.platform == "pacbio":
        platform_flag = ""  # hifiasm uses default mode for CLR (less common)

    cmd = f"cd hifiasm && hifiasm -o hifiasm.asm -t {runner.threads} {platform_flag} {runner.fastq} 2> hifiasm.log"
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
            # Core formula matches TACO.sh:
            #   BUSCO_S*1000 + T2T*300 + single*150 + MerquryComp*200
            #   + MerquryQV*20 - contigs*30 + log10(N50)*150
            # Enhancement: explicit BUSCO_D penalty (-500 per % duplication)
            # prevents highly duplicated assemblies from being selected even
            # when they have many telomeric contigs.
            score = (
                busco_s * 1000
                + t2t * 300
                + single * 150
                + merqury_comp * 200
                + merqury_qv * 20
                - contigs * 30
                + math.log10(n50) * 150
                - busco_d * 500
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


def step_12_refine(runner):
    """Step 12 - Final assembly refinement with optimized telomere-end replacement."""
    runner.log("Step 12 - Final assembly refinement with optimized telomere-end replacement")
    os.makedirs("merqury", exist_ok=True)
    os.makedirs("assemblies", exist_ok=True)

    # ---- 12A. Optional Merqury pre-selection ----
    if runner.merqury_enable:
        _run_merqury_preselection(runner)
    else:
        runner.log_info("Merqury disabled for this run")

    # Always write Merqury summary CSV
    _write_merqury_csv()

    # Build assembly info (combines BUSCO + QUAST + telomere + Merqury)
    _build_assembly_info(runner)

    # ---- 12B. Auto-select backbone assembler ----
    assembler = runner.assembler
    if not assembler:
        runner.log_info(f"Selection criteria: BUSCO + telomere + Merqury + contiguity + N50")
        assembler = _auto_select_backbone(runner)
        if assembler:
            runner.log_info(f"Auto-selected assembler: {assembler}")
        else:
            runner.log_warn("Auto-selection failed; using first available assembler as fallback")

    # Fallback: pick first available assembler
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
        runner.log_error(f"Selected assembler '{assembler}' does not have valid FASTA at {asm_fa}")
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

    runner.log_info(f"Protected mode for final refinement: {protected_mode}")
    runner.log_info(f"Selected backbone assembly: {asm_fa}")

    # ---- 12C. Prepare cleaned backbone ----
    backbone_fa = "assemblies/backbone.clean.fa"
    _clean_backbone_headers(asm_fa, backbone_fa)

    # ---- 12C2. Chimera safety check on protected pool ----
    # Even with Step 10's validated-merge approach, run a soft sanity check.
    # Use the max individual contig length across ALL assembler .telo.fasta
    # files as reference.  A legitimate T2T contig should be ≤ 1.5× this max
    # (allowing for assembly variation and gap-resolution extensions).
    # A chimera joining two chromosomes would be ~2× and gets caught.
    if protected_fasta and os.path.isfile(protected_fasta) and os.path.getsize(protected_fasta) > 0:
        # Find max contig length across all original assembler outputs
        max_individual_len = 0
        for telo_fa in glob.glob("assemblies/*.telo.fasta"):
            for _, seq in _read_fasta_records(telo_fa):
                max_individual_len = max(max_individual_len, len(seq))

        # Fallback to backbone if no telo files found
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
                runner.log_warn(
                    f"Chimera safety: removed {n_removed} protected contigs > "
                    f"{chimera_threshold:,} bp (1.5× max individual contig "
                    f"{max_individual_len:,} bp) — likely chimeras"
                )
                _write_fasta(filtered_recs, protected_fasta)
            else:
                runner.log(f"Chimera safety: all {n_before} protected contigs pass "
                           f"(threshold {chimera_threshold:,} bp)")

    # ---- 12D. Remove backbone contigs redundant to protected telomere pool ----
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

            # Use TACO.sh default thresholds (95%/95%).  Only drop backbone
            # contigs that are nearly identical to protected T2T contigs.
            # Lower thresholds (60%/90%) were tried but WORSENED duplication
            # by dropping clean backbone contigs and replacing them with dirty
            # pool contigs from high-D assemblers.
            # Pass 1 — strict: drop near-identical contigs (95%/95%)
            n_dropped = _filter_redundant_to_protected(
                paf, backbone_fa, "assemblies/backbone.filtered.fa",
                cov_thr=float(os.environ.get("PROTECT_COV", "0.95")),
                id_thr=float(os.environ.get("PROTECT_ID", "0.95")),
            )
            runner.log(f"Pass 1 (strict): removed {n_dropped} backbone contigs "
                       f"redundant to protected pool (95%/95%)")

            # Pass 2 — minimap2-based fragment removal: drop backbone
            # fragments that partially overlap T2T chromosomes (≥50%
            # coverage at ≥90% identity).  Redundans runs later on the
            # full combined assembly (12G2) where it can see both T2T
            # and backbone contigs for smarter redundancy detection.
            frag_paf = "assemblies/backbone_frag_vs_protected.paf"
            cmd = (f"minimap2 -x asm20 -t {runner.threads} "
                   f"assemblies/protected.telomere.fa assemblies/backbone.filtered.fa")
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            with open(frag_paf, "w") as f:
                f.write(result.stdout)

            n_frag = _filter_fragments_to_protected(
                frag_paf,
                "assemblies/backbone.filtered.fa",
                "assemblies/backbone.filtered.defrag.fa",
                frag_cov_thr=float(os.environ.get("FRAG_COV", "0.50")),
                frag_id_thr=float(os.environ.get("FRAG_ID", "0.90")),
            )
            if n_frag > 0:
                runner.log(f"Pass 2 (fragment): removed {n_frag} backbone "
                           f"fragments partially covered by T2T contigs (50%/90%)")
                shutil.copy("assemblies/backbone.filtered.defrag.fa",
                            "assemblies/backbone.filtered.fa")
            else:
                runner.log("Pass 2 (fragment): no additional fragments removed")
        else:
            runner.log_warn("minimap2 not found; cannot filter redundant backbone contigs.")
            shutil.copy(backbone_fa, "assemblies/backbone.filtered.fa")

        # Name-based dedup: remove backbone contigs whose names match protected IDs
        _name_dedup_fasta("assemblies/protected.telomere.fa",
                          "assemblies/backbone.filtered.fa",
                          "assemblies/backbone.filtered.nodup.fa")
    else:
        runner.log_warn("No protected telomere contigs found; keeping backbone for telomere rescue only")
        with open("assemblies/protected.telomere.fa", "w") as f:
            pass
        shutil.copy(backbone_fa, "assemblies/backbone.filtered.nodup.fa")

    # ---- 12E. Telomere rescue using optimized best single-end pool ----
    single_tel_src = ""
    for candidate in ["single_tel_best_clean.fasta", "single_tel_clean.fasta", "single_tel.fasta"]:
        if os.path.isfile(candidate) and os.path.getsize(candidate) > 0:
            single_tel_src = candidate
            break

    if single_tel_src:
        runner.log_info(f"Running telomere rescue using {single_tel_src}")
        shutil.copy("assemblies/backbone.filtered.nodup.fa", "assemblies/backbone.telomere_rescued.fa")

        if shutil.which("minimap2"):
            rescue_paf = "assemblies/single_tel_vs_backbone.paf"
            cmd = (f"minimap2 -x asm20 -t {runner.threads} "
                   f"assemblies/backbone.telomere_rescued.fa {single_tel_src}")
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            with open(rescue_paf, "w") as f:
                f.write(result.stdout)

            n_rescued = _telomere_rescue(
                single_tel_src,
                "assemblies/backbone.telomere_rescued.fa",
                rescue_paf,
                "assemblies/single_tel.replaced.ids",
                "assemblies/backbone.telomere_rescued.next.fa",
            )
            runner.log(f"Telomere rescue: {n_rescued} backbone contigs replaced with telomeric versions")

            if os.path.isfile("assemblies/backbone.telomere_rescued.next.fa") and \
               os.path.getsize("assemblies/backbone.telomere_rescued.next.fa") > 0:
                shutil.move("assemblies/backbone.telomere_rescued.next.fa",
                            "assemblies/backbone.telomere_rescued.fa")
        else:
            runner.log_warn("minimap2 not found; skipping telomere rescue.")
    else:
        runner.log_info("No single-end telomeric contigs available for telomere rescue")
        shutil.copy("assemblies/backbone.filtered.nodup.fa", "assemblies/backbone.telomere_rescued.fa")

    # ---- 12F. Remove rescued contigs redundant to protected set ----
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
            dedup_paf,
            "assemblies/backbone.telomere_rescued.fa",
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

    # ---- 12G2. Redundans reduction + scaffolding + gap closing ----
    # Run Redundans on the FULL combined assembly (T2T + backbone) so it can
    # see both the protected chromosomes and the surviving backbone fragments.
    # Running on backbone alone (as was tried earlier) produces "Nothing
    # reduced!" because the redundancy is between backbone fragments and T2T
    # contigs, not among the backbone fragments themselves.
    #
    # Redundans performs three steps in order:
    #   1. Reduction — detect and remove heterozygous/duplicate contigs
    #   2. Scaffolding — join fragments using long reads (+ reference if given)
    #   3. Gap closing — fill gaps created during scaffolding
    #
    # If --reference was provided, it is passed via -r for reference-guided
    # scaffolding.  Without --reference, Redundans uses only the long reads.
    if shutil.which("redundans.py") and os.path.isfile(raw_out) and \
       os.path.getsize(raw_out) > 0:
        runner.log_info("Running Redundans (reduction + scaffolding + gap closing) "
                        "on full combined assembly")
        runner.log_version("redundans.py", "redundans.py")

        redundans_full_dir = "assemblies/redundans_full"
        if os.path.isdir(redundans_full_dir):
            shutil.rmtree(redundans_full_dir)

        # Determine minimap2 preset based on platform
        mm2_preset = "map-hifi"
        if hasattr(runner, 'platform') and runner.platform:
            if "nano" in runner.platform.lower():
                mm2_preset = "map-ont"
            elif "clr" in runner.platform.lower() or runner.platform.lower() == "pacbio":
                mm2_preset = "map-pb"

        # Build reference flag
        ref_flag = ""
        if runner.reference_fasta and os.path.isfile(runner.reference_fasta) \
           and os.path.getsize(runner.reference_fasta) > 0:
            ref_flag = f"-r {runner.reference_fasta} "
            runner.log_info(f"Reference-guided mode using --reference: "
                            f"{runner.reference_fasta}")
        else:
            runner.log_info("De novo mode (no --reference provided, using long reads only)")

        # Count input contigs for logging
        before_recs = list(_read_fasta_records(raw_out))
        before_n = len(before_recs)
        before_len = sum(len(s) for _, s in before_recs)

        red_full_cmd = (
            f"redundans.py "
            f"-f {raw_out} "
            f"-l {runner.fastq} "
            f"{ref_flag}"
            f"-o {redundans_full_dir} "
            f"-t {runner.threads} "
            f"--minimap2reduce "
            f"-p {mm2_preset} "
            f"--identity {os.environ.get('RED_IDENTITY', '0.51')} "
            f"--overlap {os.environ.get('RED_OVERLAP', '0.80')} "
            f"--minLength 200"
        )
        result = subprocess.run(red_full_cmd, shell=True,
                                capture_output=True, text=True)
        runner.log(f"  redundans stdout: {result.stdout.strip()}")
        if result.stderr.strip():
            runner.log(f"  redundans stderr: {result.stderr.strip()[:500]}")

        # Redundans output priority: scaffolds.filled.fa > scaffolds.fa > contigs.reduced.fa
        final_red_fa = None
        for candidate in [
            os.path.join(redundans_full_dir, "scaffolds.filled.fa"),
            os.path.join(redundans_full_dir, "scaffolds.fa"),
            os.path.join(redundans_full_dir, "contigs.reduced.fa"),
        ]:
            if os.path.isfile(candidate) and os.path.getsize(candidate) > 0:
                final_red_fa = candidate
                break

        if final_red_fa:
            after_recs = list(_read_fasta_records(final_red_fa))
            after_n = len(after_recs)
            after_len = sum(len(s) for _, s in after_recs)
            runner.log(f"Redundans result: {before_n} → {after_n} contigs/scaffolds, "
                       f"{before_len:,} → {after_len:,} bp "
                       f"(from {os.path.basename(final_red_fa)})")
            shutil.copy(final_red_fa, raw_out)
        else:
            runner.log_warn("Redundans produced no output; keeping pre-Redundans assembly")
    elif not shutil.which("redundans.py"):
        runner.log_info("redundans.py not found; skipping Redundans reduction/scaffolding/gap closing")

    # ---- 12H. Genome-size-aware pruning ----
    # If total assembly size exceeds the expected genome size by >15%,
    # the smallest non-T2T backbone contigs (most likely redundant fragments)
    # are removed until the assembly is within budget.
    # This catches fragments that passed both the strict and relaxed filters
    # but still inflate assembly size and BUSCO D.
    expected_size = _parse_genome_size(runner.genomesize)
    if expected_size > 0:
        all_recs = list(_read_fasta_records(raw_out))
        total_len = sum(len(s) for _, s in all_recs)
        budget = int(expected_size * 1.15)  # 15% tolerance

        if total_len > budget:
            runner.log_warn(
                f"Assembly size {total_len:,} bp exceeds expected "
                f"{expected_size:,} bp + 15% = {budget:,} bp"
            )
            # Identify protected T2T contig names (they must not be removed)
            protected_names = set()
            if os.path.isfile(prot):
                for name, _ in _read_fasta_records(prot):
                    protected_names.add(name)

            # Split into protected (keep unconditionally) and removable
            keep_always = [(n, s) for n, s in all_recs if n in protected_names]
            removable = [(n, s) for n, s in all_recs if n not in protected_names]

            # Sort removable by length ascending — remove smallest first
            # (smallest fragments are most likely to be redundant junk)
            removable.sort(key=lambda x: len(x[1]))

            keep_len = sum(len(s) for _, s in keep_always)
            kept_backbone = []
            running_len = keep_len
            # Walk from largest to smallest, adding back until budget reached
            for name, seq in reversed(removable):
                if running_len + len(seq) <= budget:
                    kept_backbone.append((name, seq))
                    running_len += len(seq)
                else:
                    runner.log(f"  Pruning '{name}' ({len(seq):,} bp) — "
                               f"would exceed genome-size budget")

            n_pruned = len(removable) - len(kept_backbone)
            if n_pruned > 0:
                final_recs = keep_always + kept_backbone
                _write_fasta(final_recs, raw_out)
                new_total = sum(len(s) for _, s in final_recs)
                runner.log(f"Genome-size pruning: removed {n_pruned} small backbone "
                           f"fragments ({total_len:,} → {new_total:,} bp)")
            else:
                runner.log("Genome-size pruning: no contigs removed")
        else:
            runner.log(f"Assembly size {total_len:,} bp within budget "
                       f"({budget:,} bp) — no pruning needed")

    n_sorted = _fasta_sort_minlen(raw_out,
                                   f"merged_{assembler}_sort.fa",
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

    # Move merge intermediates
    for pattern in ["aln_summary_merged*.tsv", "anchor_summary_merged_*.txt"]:
        for f in glob.glob(pattern):
            shutil.move(f, "temp/merge/")
    for pattern in [".merged_*", "merged_*.delta", "merged_*.coords", "merged_*.snps",
                    "merged_*.delta.*", "merged_*.crunch", "merged_*.filter",
                    "merged_*.qdiff", "merged_*.rdiff", "merged_*.mcoords"]:
        for f in glob.glob(pattern):
            shutil.move(f, "temp/merge/")

    # Move BUSCO logs
    if os.path.isdir("busco"):
        for f in glob.glob("busco/**/*.log", recursive=True):
            try:
                shutil.move(f, "temp/busco/")
            except Exception:
                pass
    for f in glob.glob("busco_*.log") + glob.glob("*busco*.log"):
        try:
            shutil.move(f, "temp/busco/")
        except Exception:
            pass

    # Move merged FASTA intermediates
    for f in glob.glob("merged_*.fasta") + glob.glob("merged_*.fa"):
        try:
            shutil.move(f, "temp/merge/fasta/")
        except Exception:
            pass

    # Move param files
    for f in glob.glob("param_summary_merged_*.txt"):
        try:
            shutil.move(f, "temp/merge/param/")
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
