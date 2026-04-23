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

    Returns:
        int: number of contigs written (if name_map_path is not set)
        OR (int, dict): (count, {new_name -> old_name}) when used by step_12
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


def _fasta_sort_minlen_with_map(infa, outfa, prefix="contig", minlen=0):
    """Like _fasta_sort_minlen, but also returns a name mapping dict.

    Returns (count, name_map) where name_map = {new_name: old_name}.
    """
    if not os.path.isfile(infa) or os.path.getsize(infa) == 0:
        with open(outfa, 'w') as f:
            pass
        return 0, {}

    recs = [(name, seq) for name, seq in _read_fasta_records(infa)
            if len(seq) >= minlen]
    recs.sort(key=lambda x: len(x[1]), reverse=True)

    name_map = {}
    renamed = []
    for i, (old_name, seq) in enumerate(recs, start=1):
        new_name = f"{prefix}_{i}"
        renamed.append((new_name, seq))
        name_map[new_name] = old_name
    _write_fasta(renamed, outfa)
    return len(renamed), name_map


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

    # Run canu with stderr/stdout captured so we can report the actual error.
    # Canu dev builds from bioconda frequently fail with broken Java (JLI_StringDup).
    runner.log(f"Running canu")
    runner.log(f"$ {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        # Log captured stderr/stdout for diagnosis
        err_snippet = (result.stderr or "").strip()
        out_snippet = (result.stdout or "").strip()
        if err_snippet:
            # Show first ~20 lines of stderr
            err_lines = err_snippet.split("\n")[:20]
            for line in err_lines:
                runner.log_warn(f"  canu stderr: {line}")
        if out_snippet and not err_snippet:
            out_lines = out_snippet.split("\n")[:10]
            for line in out_lines:
                runner.log_warn(f"  canu stdout: {line}")

        # Detect common failure modes and give actionable guidance
        combined = (err_snippet + out_snippet).lower()
        if "jli_stringdup" in combined or "undefined symbol" in combined:
            runner.log_warn("This is the known bioconda canu Java runtime bug. "
                            "Fix: install a stable canu binary from "
                            "https://github.com/marbl/canu/releases "
                            "(not from conda) and place it on PATH.")
        elif "java" in combined and ("error" in combined or "exception" in combined):
            runner.log_warn("Canu Java runtime error. Try: "
                            "conda install -c conda-forge 'openjdk>=17' "
                            "or install canu from GitHub releases.")
        elif "no reads" in combined or "no input" in combined:
            runner.log_warn("Canu could not find input reads. "
                            f"Check that {runner.fastq} exists and is readable.")

        runner.log_warn("Step 1: canu failed. Skipping. Other assemblers will continue.")


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
    # Exclude leftover final assembly files from previous runs — these should
    # not be treated as assembler inputs for the pool.
    _exclude_prefixes = {"final", "final_merge", "backbone"}
    fasta_files = [f for f in fasta_files
                   if os.path.basename(f).replace(".telo.fasta", "")
                   not in _exclude_prefixes]
    if not fasta_files:
        runner.log_info("No per-assembly *.telo.fasta files found; attempting to continue")
        fasta_files = glob.glob("assemblies/*.result.fasta")
        fasta_files = [f for f in fasta_files
                       if os.path.basename(f).replace(".result.fasta", "")
                       not in _exclude_prefixes]
        if not fasta_files:
            runner.log_error("Still no FASTA files after generation.")
            raise RuntimeError("No FASTA files found")

    # ===== Telomere-aware validated quickmerge =====
    # Run quickmerge pairwise, then validate: only accept merged contigs that
    # (a) classify as strict_t2t and (b) pass length sanity check.
    validated_t2t_from_merge = []  # list of (name, seq) tuples
    qm_pair_map = {}  # quickmerge_contig_name -> (assembler1, assembler2)
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
                        for tname, tseq in new_t2t:
                            qm_pair_map[tname] = (base1, base2)
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

    n_sorted, pool_name_map = _fasta_sort_minlen_with_map(
        "allmerged_telo.fasta", "allmerged_telo_sort.fasta",
        prefix="contig", minlen=500)
    runner.log(f"Sorted telomere contigs: {n_sorted} contigs >= 500 bp")

    # Build comprehensive provenance map:
    # pool_contig_name → (assembler, original_name, source_type)
    #
    # The concatenation order is tracked in contig_order (assembler_key, orig_name).
    # Quickmerge contigs are appended after with suffix "_qm_validated".
    # pool_name_map gives {new_sorted_name: old_concat_name}.
    # We need: old_concat_name → (assembler, original_name, source_type).
    concat_provenance = {}  # old_concat_name → (asm, orig, source_type)
    for f_path in sorted(fasta_files):
        basename = os.path.basename(f_path).replace(".telo.fasta", "").replace(
            ".result.fasta", "").replace(".fasta", "")
        asm_key = basename.lower()
        if os.path.isfile(f_path):
            for cname, _ in _read_fasta_records(f_path):
                concat_provenance[cname] = (asm_key, cname, "assembler")
    # ---- Region-level provenance for quickmerge contigs ----
    # For each validated quickmerge T2T, align it against both source assembler
    # FASTAs to determine which assembler contributed which region.
    # qm_regions: qm_contig_name -> list of (start, end, assembler, source_contig)
    qm_regions = {}
    if validated_t2t_from_merge and qm_pair_map and shutil.which("minimap2"):
        # Write all validated QM contigs to a temporary FASTA
        qm_tmp_fa = "assemblies/_qm_validated_tmp.fa"
        _write_fasta(validated_t2t_from_merge, qm_tmp_fa)

        for qm_name, (asm1, asm2) in qm_pair_map.items():
            regions = []
            qm_seq_len = 0
            for tname, tseq in validated_t2t_from_merge:
                if tname == qm_name:
                    qm_seq_len = len(tseq)
                    break
            if qm_seq_len == 0:
                continue

            # Write just this one QM contig
            qm_single_fa = "assemblies/_qm_single_tmp.fa"
            for tname, tseq in validated_t2t_from_merge:
                if tname == qm_name:
                    _write_fasta([(tname, tseq)], qm_single_fa)
                    break

            # Align against each source assembler
            for asm_name in [asm1, asm2]:
                asm_fa = f"assemblies/{asm_name}.telo.fasta"
                if not os.path.isfile(asm_fa):
                    asm_fa = f"assemblies/{asm_name}.result.fasta"
                if not os.path.isfile(asm_fa):
                    continue

                paf_tmp = "assemblies/_qm_region_tmp.paf"
                cmd = (f"minimap2 -x asm20 -t {runner.threads} "
                       f"{qm_single_fa} {asm_fa}")
                result = subprocess.run(cmd, shell=True, capture_output=True,
                                        text=True)
                with open(paf_tmp, "w") as f:
                    f.write(result.stdout)

                # Parse PAF: find which target regions of the QM contig each
                # source assembler covers.  PAF target = qm_contig (reference),
                # query = assembler contig.
                if os.path.isfile(paf_tmp):
                    with open(paf_tmp) as f:
                        for ln in f:
                            p = ln.rstrip().split("\t")
                            if len(p) < 12:
                                continue
                            src_contig = p[0]
                            tstart, tend = int(p[7]), int(p[8])
                            matches, alnlen = int(p[9]), int(p[10])
                            if alnlen <= 0:
                                continue
                            ident = matches / max(1, alnlen)
                            if ident >= 0.85 and (tend - tstart) >= 1000:
                                regions.append((tstart, tend, asm_name,
                                                src_contig, ident))

            # Merge overlapping regions per assembler and pick best coverage
            if regions:
                # Sort by start position
                regions.sort(key=lambda x: x[0])
                # Deduplicate: for overlapping regions from different assemblers,
                # keep both (they show the merge boundary).  Just store all.
                qm_regions[qm_name] = regions

        # Cleanup temp files
        for tmp in ["assemblies/_qm_validated_tmp.fa",
                     "assemblies/_qm_single_tmp.fa",
                     "assemblies/_qm_region_tmp.paf"]:
            if os.path.isfile(tmp):
                os.remove(tmp)

    for name, seq in validated_t2t_from_merge:
        qm_name = f"{name}_qm_validated"
        asm1, asm2 = qm_pair_map.get(name, ("unknown", "unknown"))
        concat_provenance[qm_name] = (
            "quickmerge", name, "quickmerge",
            asm1, asm2,
            qm_regions.get(name, [])
        )

    # Now map pool sorted names to provenance
    # Extended format: (asm, orig_name, source_type, [asm1, asm2, regions])
    pool_provenance = {}
    for new_name, old_name in pool_name_map.items():
        if old_name in concat_provenance:
            prov = concat_provenance[old_name]
            if len(prov) == 3:
                # Regular assembler contig: (asm, orig, source_type)
                pool_provenance[new_name] = prov
            else:
                # Quickmerge contig: (asm, orig, source_type, asm1, asm2, regions)
                pool_provenance[new_name] = prov
        else:
            pool_provenance[new_name] = ("unknown", old_name, "unknown")

    # Save provenance TSV for Step 12 GFF generation
    # Extended format with quickmerge pair info and region detail
    prov_tsv = "pool_contig_provenance.tsv"
    with open(prov_tsv, "w") as f:
        f.write("pool_name\tassembler\toriginal_name\tsource_type\t"
                "qm_assembler1\tqm_assembler2\tqm_regions\n")
        for pname in sorted(pool_provenance.keys(),
                            key=lambda x: int(x.split("_")[1]) if "_" in x else 0):
            prov = pool_provenance[pname]
            asm, orig, stype = prov[0], prov[1], prov[2]
            if len(prov) >= 6 and stype == "quickmerge":
                asm1, asm2, regions = prov[3], prov[4], prov[5]
                # Format regions as "start-end:asm:contig;..."
                reg_str = ";".join(
                    f"{s}-{e}:{a}:{c}"
                    for s, e, a, c, _ in regions
                ) if regions else ""
                f.write(f"{pname}\t{asm}\t{orig}\t{stype}\t{asm1}\t{asm2}\t{reg_str}\n")
            else:
                f.write(f"{pname}\t{asm}\t{orig}\t{stype}\t\t\t\n")
    runner.log(f"Saved pool contig provenance: {prov_tsv} ({len(pool_provenance)} entries)")

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
        missing = [n for n in ids if n not in pool_seqs]
        if missing:
            runner.log_warn(f"Name mismatch writing {fname}: {len(missing)}/{len(ids)} "
                            f"classified IDs not found in pool_seqs "
                            f"(first 3: {missing[:3]})")
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
    # Map pool contigs to source assembler using the provenance map
    # (accurately tracked through concatenation and sort in Step 10).
    pool_contig_asm = {}  # pool_contig_name -> asm_key
    for pname, (asm, orig, stype) in pool_provenance.items():
        pool_contig_asm[pname] = asm

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

    # Count contigs from final output files (after clean_contained dedup)
    strict_t2t_n = sum(1 for _ in open("t2t.fasta") if _.startswith(">")) if os.path.isfile("t2t.fasta") else 0
    single_tel_n = sum(1 for _ in open("single_tel.fasta") if _.startswith(">")) if os.path.isfile("single_tel.fasta") else 0
    tel_supported_n = sum(1 for _ in open("telomere_supported.fasta") if _.startswith(">")) if os.path.isfile("telomere_supported.fasta") else 0

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

    # Taxon-aware scoring weights (defined outside loop so they're accessible
    # for the decision file regardless of which assemblers were scored).
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
        # Note: MerquryComp and MerquryQV contribute 0 when Merqury is not enabled
        has_merqury = any(r["merqury_qv"] > 0 or r["merqury_comp"] > 0 for r in records)
        if has_merqury:
            f.write(f"score_formula\tBUSCO_S*{w_busco_s} + T2T*{w_t2t} + single*{w_single} + MerquryComp*200 + MerquryQV*20 - contigs*{w_contigs} + log10(N50)*{w_n50} - BUSCO_D*{w_busco_d} (taxon={taxon})\n")
        else:
            f.write(f"score_formula\tBUSCO_S*{w_busco_s} + T2T*{w_t2t} + single*{w_single} - contigs*{w_contigs} + log10(N50)*{w_n50} - BUSCO_D*{w_busco_d} (taxon={taxon}, Merqury not available)\n")

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


def _filter_redundant_to_protected(paf_path, backbone_fa, out_fa, cov_thr=0.95, id_thr=0.95,
                                    logger=None):
    """Remove backbone contigs that are redundant to the protected telomere pool.

    Uses PAF from minimap2 alignment of backbone vs protected pool.
    Keeps backbone contigs that are NOT covered >= cov_thr with identity >= id_thr.

    Returns (n_dropped, drop_details) where drop_details is a list of
    (backbone_name, backbone_len, cov, ident) tuples for logging.
    """
    best = _parse_paf_best_hits(paf_path)
    drop = {}  # name -> (cov, ident)
    for q, (cov, ident) in best.items():
        if cov >= cov_thr and ident >= id_thr:
            drop[q] = (cov, ident)

    drop_details = []
    recs = []
    for n, s in _read_fasta_records(backbone_fa):
        if n in drop:
            cov, ident = drop[n]
            drop_details.append((n, len(s), cov, ident))
            if logger:
                logger(f"  Dedup drop: {n} ({len(s):,} bp, "
                       f"cov={cov:.3f}, ident={ident:.3f})")
        else:
            recs.append((n, s))
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


def _write_provenance_gff(final_fa, gff_path, name_map, protected_ids,
                          replaced_map, backbone_assembler, pool_asm_map,
                          pool_provenance_map=None, runner=None):
    """Write a GFF3 file documenting the provenance of each contig in the final assembly.

    Each contig gets a GFF3 record (type=contig) spanning its full length, with attributes:
      - source_assembler: which assembler originally produced this contig
      - role: backbone, upgrade_donor, novel_t2t
      - original_name: contig name before sort/rename (the pool contig name)
      - assembler_contig: the ORIGINAL assembler contig name (before Step 10 pool renaming)
      - source_type: assembler or quickmerge
      - replacement_class: (for donors) upgrade_tier2_to_t2t, fill_missing_end, etc.
      - replaced_contig: (for donors) which backbone contig was replaced
      - qm_assembler1, qm_assembler2: (for quickmerge) the two source assemblers
      - description: human-readable provenance summary
      - length: contig length in bp

    For quickmerge contigs, child records (type=region) are emitted showing which
    assembler contributed each genomic region, with attributes:
      - Parent: the parent contig ID
      - source_assembler: which assembler contributed this region
      - assembler_contig: the original contig name from that assembler
      - description: human-readable region summary (e.g., "Region 1-500000 from canu contig 'tig00001'")

    Args:
        final_fa: path to final renamed FASTA
        gff_path: output GFF3 path
        name_map: {new_contig_name: old_contig_name} from _fasta_sort_minlen_with_map
        protected_ids: set of contig names from the protected T2T pool
        replaced_map: {donor_name: (backbone_name, replacement_class)}
        backbone_assembler: name of the selected backbone assembler
        pool_asm_map: {contig_name: assembler_key}
        pool_provenance_map: {pool_name: (asm, orig, source_type[, asm1, asm2, regions])}
            Normal entries are 3-tuples; quickmerge entries are 6-tuples with
            source assembler pair and region-level mapping.
        runner: pipeline runner for logging
    """
    if not os.path.isfile(final_fa) or os.path.getsize(final_fa) == 0:
        return

    if pool_provenance_map is None:
        pool_provenance_map = {}

    donor_names = set(replaced_map.keys())

    lines = ["##gff-version 3"]
    lines.append(f"# TACO provenance annotation for {os.path.basename(final_fa)}")
    lines.append(f"# Backbone assembler: {backbone_assembler}")
    lines.append(f"# Generated by TACO v1.2.0")
    lines.append(f"# Columns: seqid source type start end score strand phase attributes")
    lines.append("#")

    n_bb = n_upgrade = n_novel = 0

    for new_name, seq in _read_fasta_records(final_fa):
        length = len(seq)
        old_name = name_map.get(new_name, new_name)

        # purge_dups renames contigs by appending _0, _1, etc. suffixes
        # (e.g., "contig_14" becomes "contig_14_0" or "peregrine_3" becomes
        # "peregrine_3_1").  Strip this suffix for provenance/role lookups.
        lookup_name = old_name
        if lookup_name not in pool_provenance_map and \
           lookup_name not in protected_ids and \
           lookup_name not in donor_names:
            # Try stripping trailing _N suffix (purge_dups artifact)
            m = re.match(r'^(.+)_(\d+)$', lookup_name)
            if m:
                base = m.group(1)
                if base in pool_provenance_map or \
                   base in protected_ids or \
                   base in donor_names:
                    lookup_name = base

        # Look up full provenance from Step 10 TSV
        # Quickmerge entries are 6-tuples: (asm, orig, source_type, asm1, asm2, regions)
        # Normal entries are 3-tuples: (asm, orig, source_type)
        prov = pool_provenance_map.get(lookup_name)
        qm_asm1 = qm_asm2 = ""
        qm_regions = []
        if prov and len(prov) >= 6:
            orig_asm, orig_contig, source_type = prov[0], prov[1], prov[2]
            qm_asm1, qm_asm2 = prov[3], prov[4]
            qm_regions = prov[5] if prov[5] else []
        elif prov:
            orig_asm, orig_contig, source_type = prov[0], prov[1], prov[2]
        else:
            orig_asm = "unknown"
            orig_contig = lookup_name
            source_type = "unknown"

        # Determine role, source assembler, and build description
        # Use lookup_name (purge_dups suffix stripped) for all set/dict lookups
        if lookup_name in donor_names:
            bb_name, repl_class = replaced_map[lookup_name]
            role = "upgrade_donor"
            source_asm = pool_asm_map.get(lookup_name, orig_asm)
            n_upgrade += 1

            if source_type == "quickmerge" and qm_asm1 and qm_asm2:
                desc = (f"Entire replacement: {backbone_assembler} contig "
                        f"'{bb_name}' replaced by quickmerge of "
                        f"{qm_asm1} x {qm_asm2} contig "
                        f"'{orig_contig}' (class: {repl_class})")
            elif source_type == "quickmerge":
                desc = (f"Entire replacement: {backbone_assembler} contig "
                        f"'{bb_name}' replaced by quickmerge contig "
                        f"'{orig_contig}' (class: {repl_class})")
            else:
                desc = (f"Entire replacement: {backbone_assembler} contig "
                        f"'{bb_name}' replaced by {source_asm} contig "
                        f"'{orig_contig}' (class: {repl_class})")
        elif lookup_name in protected_ids and lookup_name not in donor_names:
            role = "novel_t2t"
            source_asm = pool_asm_map.get(lookup_name, orig_asm)
            n_novel += 1

            if source_type == "quickmerge" and qm_asm1 and qm_asm2:
                desc = (f"Novel T2T addition from quickmerge of "
                        f"{qm_asm1} x {qm_asm2}: '{orig_contig}' "
                        f"covers a chromosomal region absent from "
                        f"{backbone_assembler} backbone")
            elif source_type == "quickmerge":
                desc = (f"Novel T2T addition from quickmerge: '{orig_contig}' "
                        f"covers a chromosomal region absent from "
                        f"{backbone_assembler} backbone")
            else:
                desc = (f"Novel T2T addition from {source_asm}: '{orig_contig}' "
                        f"covers a chromosomal region absent from "
                        f"{backbone_assembler} backbone")
        else:
            role = "backbone"
            repl_class = ""
            source_asm = backbone_assembler
            n_bb += 1

            # Check if this backbone contig has a known original name
            # (backbone contigs may also have provenance if they were
            # in the pool before being selected as backbone)
            bb_prov = pool_provenance_map.get(lookup_name)
            if bb_prov:
                if len(bb_prov) >= 6:
                    orig_asm, orig_contig, source_type = bb_prov[0], bb_prov[1], bb_prov[2]
                else:
                    orig_asm, orig_contig, source_type = bb_prov[0], bb_prov[1], bb_prov[2]
                source_asm = orig_asm
                desc = (f"Retained from {source_asm} ({backbone_assembler} backbone): "
                        f"original contig '{orig_contig}'")
            else:
                # This is a backbone contig not in the pool — it comes directly
                # from the selected assembler.  Set proper provenance fields.
                # Strip purge_dups _N suffix from backbone contig name
                bb_clean = lookup_name
                m_bb = re.match(r'^(.+)_(\d+)$', bb_clean)
                if m_bb:
                    bb_clean = m_bb.group(1)
                source_type = "assembler"
                source_asm = backbone_assembler
                orig_asm = backbone_assembler
                orig_contig = bb_clean
                desc = (f"Retained from {backbone_assembler} backbone: "
                        f"original contig '{bb_clean}'")

        # GFF3 columns: seqid source type start end score strand phase attributes
        attrs = (f"ID={new_name}"
                 f";original_name={old_name}"
                 f";assembler_contig={orig_contig}"
                 f";source_assembler={source_asm}"
                 f";source_type={source_type}"
                 f";role={role}"
                 f";length={length}"
                 f";description={desc}")
        if role == "upgrade_donor":
            attrs += f";replacement_class={repl_class};replaced_contig={bb_name}"
        if source_type == "quickmerge" and qm_asm1 and qm_asm2:
            attrs += f";qm_assembler1={qm_asm1};qm_assembler2={qm_asm2}"

        lines.append(f"{new_name}\tTACO\tcontig\t1\t{length}\t.\t.\t.\t{attrs}")

        # Emit child region-level GFF records for quickmerge contigs
        # Each region shows which assembler contributed that segment
        if source_type == "quickmerge" and qm_regions:
            for ridx, (rstart, rend, rasm, rcontig) in enumerate(
                    sorted(qm_regions, key=lambda r: r[0]), 1):
                region_id = f"{new_name}_region_{ridx}"
                region_desc = (f"Region {rstart}-{rend} from {rasm} "
                               f"contig '{rcontig}'")
                region_attrs = (
                    f"ID={region_id}"
                    f";Parent={new_name}"
                    f";source_assembler={rasm}"
                    f";assembler_contig={rcontig}"
                    f";description={region_desc}")
                lines.append(
                    f"{new_name}\tTACO\tregion\t{rstart}\t{rend}"
                    f"\t.\t.\t.\t{region_attrs}")

    # Write ##sequence-region pragmas for GFF3 compliance
    seq_regions = []
    for new_name, seq in _read_fasta_records(final_fa):
        seq_regions.append(f"##sequence-region {new_name} 1 {len(seq)}")

    header_end = 0
    for i, line in enumerate(lines):
        if not line.startswith("#"):
            header_end = i
            break
    lines = lines[:header_end] + seq_regions + lines[header_end:]

    with open(gff_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    if runner:
        runner.log(f"Wrote provenance GFF: {gff_path} "
                   f"({n_bb} backbone, {n_upgrade} upgrades, {n_novel} novel T2T)")


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
    """Write plausible upgrade/rescue candidates to TSV.

    Handles both pool T2T upgrade candidates (which have 'source' and
    'replacement_class' but no alignment stats) and single-end rescue
    candidates (which have full alignment fields like 'ident', 'cov_backbone').
    """
    cols = ["candidate_rank", "source", "backbone", "donor",
            "replacement_class", "ident", "aligned_bp",
            "cov_backbone", "cov_donor", "ext", "len_gain",
            "touches_left", "touches_right", "structural_score"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for c in candidates:
            ident = c.get("ident")
            cov_bb = c.get("cov_backbone")
            cov_d = c.get("cov_donor")
            score = c.get("structural_score")
            w.writerow([
                c.get("candidate_rank", ""),
                c.get("source", "rescue"),
                c["backbone"], c["donor"],
                c.get("replacement_class", ""),
                f"{ident:.4f}" if ident is not None else "",
                c.get("aligned_bp", ""),
                f"{cov_bb:.4f}" if cov_bb is not None else "",
                f"{cov_d:.4f}" if cov_d is not None else "",
                c.get("ext", ""), c.get("len_gain", ""),
                c.get("touches_left", ""), c.get("touches_right", ""),
                f"{score:.4f}" if score is not None else "",
            ])


def _run_busco_trial(trial_fa, lineage, threads, trial_label, out_dir):
    """Run BUSCO on a trial assembly and return parsed metrics dict or None."""
    busco_bin = shutil.which("busco")
    if not busco_bin:
        return None

    trial_out = os.path.join(out_dir, trial_label)
    if os.path.isdir(trial_out):
        shutil.rmtree(trial_out)

    # Try --offline first, fall back to normal mode if lineage not cached
    cmd = (f"busco -o {trial_label} -i {trial_fa} -l {lineage} "
           f"-m genome -c {threads} --out_path {out_dir} --offline")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        # Retry without --offline (lineage may need download)
        if os.path.isdir(trial_out):
            shutil.rmtree(trial_out)
        cmd = (f"busco -o {trial_label} -i {trial_fa} -l {lineage} "
               f"-m genome -c {threads} --out_path {out_dir}")
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

    # Last resort: parse from BUSCO stdout (contains summary line)
    for output in [result.stdout or "", result.stderr or ""]:
        m = re.search(r"C:(\d+(?:\.\d+)?)%.*?S:(\d+(?:\.\d+)?)%.*?"
                      r"D:(\d+(?:\.\d+)?)%.*?F:(\d+(?:\.\d+)?)%.*?"
                      r"M:(\d+(?:\.\d+)?)%", output, re.S)
        if m:
            return {"C": float(m.group(1)), "S": float(m.group(2)),
                    "D": float(m.group(3)), "F": float(m.group(4)),
                    "M": float(m.group(5))}
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
    # Taxon-aware purge_dups strategy:
    #   - Fungi/haploid: single-round, but adjust calcuts low-coverage cutoff
    #     downward.  In haploid genomes, true haplotigs have ~0.5× coverage
    #     but purge_dups' auto calcuts often sets the boundary too low,
    #     missing them.  Override with PURGE_DUPS_CALCUTS env var if needed.
    #   - Plant/polyploid: conservative single-round to avoid removing
    #     homeologous sequences.  Two-round purging is too aggressive.
    #   - Vertebrate/diploid: two-round purging for thorough haplotig removal.
    #   - Other: single-round, standard calcuts.
    base_cov = os.path.join(pd_dir, "PB.base.cov")
    dups_bed = os.path.join(pd_dir, "dups.bed")

    # Log coverage cutoffs for debugging
    if os.path.isfile(cutoff_file):
        with open(cutoff_file) as cf:
            cutoffs_str = cf.read().strip()
        runner.log(f"purge_dups coverage cutoffs: {cutoffs_str}")

    # Allow user override of calcuts thresholds
    custom_calcuts = os.environ.get("PURGE_DUPS_CALCUTS", "")

    purge_flag = ""
    if taxon in ("vertebrate", "animal"):
        purge_flag = "-2"  # two-round for diploid genomes
        runner.log_info("purge_dups strategy: vertebrate/diploid (two-round)")
    elif taxon == "plant":
        # Single-round, conservative — avoid collapsing homeologs
        runner.log_info("purge_dups strategy: plant (single-round, conservative)")
    elif taxon == "fungal":
        # Fungi are typically haploid — duplicates should have ~0.5× coverage.
        # Use -2 for more aggressive duplicate detection in haploid genomes
        # where coverage differences between primary and duplicate are subtle.
        purge_flag = "-2"
        runner.log_info("purge_dups strategy: fungal/haploid (two-round, aggressive)")
    else:
        runner.log_info("purge_dups strategy: default (single-round)")

    if custom_calcuts:
        # User provided custom calcuts: e.g. "5 10 15 20 25 30"
        with open(cutoff_file, "w") as cf:
            cf.write(custom_calcuts + "\n")
        runner.log_info(f"Using custom calcuts override: {custom_calcuts}")

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
      - HiFi: NextPolish2 by default (k-mer-based, safe for HiFi assemblies).
        Requires yak for k-mer database construction.  Falls back to skip if
        NextPolish2 or yak is not installed.
      - ONT:  Medaka (preferred, neural-network-based) → Racon (fallback)
      - CLR:  Racon (standard for CLR polishing)

    Returns True if polishing ran successfully.
    """
    platform = runner.platform or "pacbio-hifi"

    if platform == "pacbio-hifi":
        # NextPolish2 uses yak k-mer databases (k=21 and k=31) to polish
        # HiFi assemblies.  This is the recommended approach: k-mer-based
        # polishing is safe for HiFi data and corrects residual errors
        # without the risk of introducing new ones from read-alignment.
        np2_bin = shutil.which("nextPolish2") or shutil.which("nextpolish2")
        yak_bin = shutil.which("yak")

        if np2_bin and yak_bin:
            runner.log_info("Polishing HiFi assembly with NextPolish2")
            polish_dir = "assemblies/polish_work"
            os.makedirs(polish_dir, exist_ok=True)

            k21_yak = os.path.join(polish_dir, "k21.yak")
            k31_yak = os.path.join(polish_dir, "k31.yak")

            # Build yak k-mer databases from HiFi reads
            runner.log_info("Building yak k-mer databases (k=21, k=31)...")
            for ksize, yak_out in [(21, k21_yak), (31, k31_yak)]:
                if os.path.isfile(yak_out) and os.path.getsize(yak_out) > 1024:
                    runner.log(f"  Reusing existing {yak_out}")
                    continue
                elif os.path.isfile(yak_out):
                    runner.log(f"  Removing suspect yak file {yak_out} "
                               f"({os.path.getsize(yak_out)} bytes)")
                    os.remove(yak_out)
                # yak count requires reads passed TWICE for two-pass counting
                # (first pass: count k-mers, second pass: build database)
                cmd = (f"{yak_bin} count -k {ksize} -b 37 "
                       f"-t {runner.threads} -o {yak_out} "
                       f"{runner.fastq} {runner.fastq}")
                result = subprocess.run(cmd, shell=True, capture_output=True,
                                        text=True)
                if result.returncode != 0:
                    runner.log_warn(f"yak count k={ksize} failed: "
                                    f"{(result.stderr or '').strip()[:200]}")
                    runner.log_info("Falling back: skipping HiFi polishing")
                    shutil.copy(input_fa, output_fa)
                    return False
                runner.log(f"  Built {yak_out}")

            # Run NextPolish2
            # v0.2+ requires: nextPolish2 -t N reads.sorted.bam genome.fa k21.yak k31.yak
            # Step 1: Map HiFi reads to assembly with minimap2
            # Step 2: Sort BAM with samtools
            # Step 3: Run NextPolish2 with BAM + genome + yak databases
            abs_input = os.path.abspath(input_fa)
            abs_k21 = os.path.abspath(k21_yak)
            abs_k31 = os.path.abspath(k31_yak)

            ver_result = subprocess.run(f"{np2_bin} --version",
                                        shell=True, capture_output=True, text=True)
            np2_ver = (ver_result.stdout or ver_result.stderr or "").strip()
            runner.log_info(f"NextPolish2 version: {np2_ver or 'unknown'}")

            samtools_bin = shutil.which("samtools")
            if not samtools_bin:
                runner.log_warn("samtools not found; NextPolish2 requires sorted BAM. "
                                "Install samtools or skip with --no-polish")
                shutil.copy(input_fa, output_fa)
                return False

            # Map HiFi reads to assembly
            bam_unsorted = os.path.join(polish_dir, "hifi_vs_asm.bam")
            bam_sorted = os.path.join(polish_dir, "hifi_vs_asm.sorted.bam")
            abs_bam = os.path.abspath(bam_sorted)

            runner.log_info("Mapping HiFi reads to assembly for NextPolish2...")
            cmd = (f"minimap2 -ax map-hifi -t {runner.threads} "
                   f"{abs_input} {runner.fastq} "
                   f"| {samtools_bin} sort -@ {runner.threads} "
                   f"-o {bam_sorted}")
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if not os.path.isfile(bam_sorted) or os.path.getsize(bam_sorted) == 0:
                runner.log_warn("HiFi read mapping failed; skipping NextPolish2")
                shutil.copy(input_fa, output_fa)
                return False

            # Index the sorted BAM
            cmd = f"{samtools_bin} index {bam_sorted}"
            subprocess.run(cmd, shell=True, capture_output=True, text=True)
            runner.log("  Mapped and sorted HiFi reads")

            # Run NextPolish2: nextPolish2 -t N reads.bam genome.fa k21.yak k31.yak
            cmd = (f"{np2_bin} -t {runner.threads} "
                   f"{abs_bam} {abs_input} {abs_k21} {abs_k31}")
            runner.log_info(f"Running: {cmd}")
            result = subprocess.run(cmd, shell=True, capture_output=True,
                                    text=True)

            if result.returncode == 0 and result.stdout.strip():
                with open(output_fa, "w") as f:
                    f.write(result.stdout)
                if os.path.getsize(output_fa) > 0:
                    runner.log("NextPolish2 polishing complete")
                    return True

            # If NextPolish2 failed, check if yak databases are stale/corrupt
            err_msg = (result.stderr or "").strip()[:300]
            yak_stale = ("not a valid" in err_msg or "InvalidData" in err_msg
                         or "panicked" in err_msg)

            if yak_stale:
                runner.log_warn(
                    f"NextPolish2 failed (likely stale/incompatible yak databases). "
                    f"Removing cached .yak files — they will be rebuilt on next run. "
                    f"Error: {err_msg}")
                for yf in [k21_yak, k31_yak]:
                    if os.path.isfile(yf):
                        os.remove(yf)
                        runner.log(f"  Removed stale {yf}")
            else:
                runner.log_warn(f"NextPolish2 failed or produced no output"
                                f"{': ' + err_msg if err_msg else ''}")

            shutil.copy(input_fa, output_fa)
            return False
        else:
            missing = []
            if not np2_bin:
                missing.append("nextPolish2")
            if not yak_bin:
                missing.append("yak")
            runner.log_warn(f"HiFi polishing requires {' and '.join(missing)}; "
                            f"install via conda (nextpolish2, yak) or skip with "
                            f"--no-polish")
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

    v1.2.0 workflow — two-tier confidence model:

    TIER 1 (immutable): protected strict-T2T contigs from the merged pool.
      These are NEVER replaced or displaced by rescue candidates unless the
      user explicitly passes --allow-t2t-replace.

    TIER 2 (editable): backbone contigs from the auto-selected assembler.
      These can be deduped, replaced, rescued, purged, and polished.

    Sub-steps:
      12A  Merqury QV scoring (optional, auto-detected)
      12B  Auto-select backbone assembler
      12C  Prepare cleaned backbone + chimera safety check
      12D  T2T-first foundation building:
           12D1  Strict dedup — remove backbone contigs 95%/95% redundant to Tier 1
           12D2  Fragment removal — remove backbone fragments 50%/90% to Tier 1
           12D3  Classify backbone contigs by telomere status
           12D4  Taxon-aware non-telomeric dedup against Tier 1
                 (fungi: aggressive 70%/85%; plant/vertebrate: conservative 85%/92%)
           12D5  Taxon-aware non-telomeric self-dedup
                 (fungi: 80%/90%; plant/vertebrate: 90%/95%)
      12E  Telomere rescue with donor verification + replacement class tagging:
           only accept donors with verified telomere signal;
           never replace Tier 1 contigs (unless --allow-t2t-replace)
      12F  BUSCO trial validation:
           taxon-aware thresholds for C-drop, M-rise, AND D-rise;
           reject if telomere evidence weakens or size drops suspiciously
      12G  Final combine (Tier 1 foundation + Tier 2 rescued backbone)
      12H  purge_dups default cleanup (taxon-aware)
      12I  Platform-aware polishing (Medaka/NextPolish2/Racon)
      12J  Genome-size-aware pruning (never prunes telomere-bearing contigs)
    """
    runner.log("Step 12 - Telomere-aware backbone refinement (v1.2.0)")
    os.makedirs("merqury", exist_ok=True)
    os.makedirs("assemblies", exist_ok=True)

    # ---- 12A. Merqury QV scoring (optional, auto-detected if installed) ----
    if runner.merqury_enable:
        db_src = runner.merqury_db or "(auto-detected)"
        runner.log_info(f"Merqury enabled (db: {db_src})")
        _run_merqury_preselection(runner)
    else:
        runner.log_info("Merqury scoring not available (install merqury + provide .meryl db, "
                        "or use --merqury/--merqury-db to enable)")

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
    # Chimera detection uses two complementary strategies:
    #   1. Size gate: contigs > 1.5× the largest individual assembler contig
    #      are likely chimeric (spurious joins across chromosomes).
    #   2. Cross-assembly mapping: each protected contig is aligned against
    #      all other assembler outputs.  A non-chimeric contig should be
    #      largely covered (≥60%) by at least one other assembler's contigs.
    #      A chimera that fused two separate chromosomes won't be fully
    #      covered by any single assembler's individual contig — it will
    #      only show partial hits.
    backbone_fa = "assemblies/backbone.clean.fa"
    _clean_backbone_headers(asm_fa, backbone_fa)

    if protected_fasta and os.path.isfile(protected_fasta) and os.path.getsize(protected_fasta) > 0:
        # ---- Strategy 1: Size gate ----
        max_individual_len = 0
        for telo_fa in glob.glob("assemblies/*.telo.fasta"):
            for _, seq in _read_fasta_records(telo_fa):
                max_individual_len = max(max_individual_len, len(seq))
        if max_individual_len == 0:
            for _, seq in _read_fasta_records(backbone_fa):
                max_individual_len = max(max_individual_len, len(seq))

        size_suspects = set()
        if max_individual_len > 0:
            chimera_threshold = int(max_individual_len * 1.5)
            protected_recs = list(_read_fasta_records(protected_fasta))
            n_before = len(protected_recs)
            for n, s in protected_recs:
                if len(s) > chimera_threshold:
                    size_suspects.add(n)
        else:
            protected_recs = list(_read_fasta_records(protected_fasta))
            n_before = len(protected_recs)
            chimera_threshold = 0

        # ---- Strategy 2: Cross-assembly mapping consistency ----
        # Align protected contigs against each assembler's output.  For each
        # protected contig, track the BEST single-contig coverage from any
        # assembler.  A legitimate T2T contig should be well-covered by at
        # least one other assembler's contig.  Chimeras will only show
        # partial coverage (two halves matching different contigs).
        mapping_suspects = set()
        if shutil.which("minimap2") and len(protected_recs) > 0:
            # Collect all assembler result FASTAs (excluding the selected backbone)
            asm_fastas = []
            for asm_name in ["canu", "flye", "nextDenovo", "peregrine", "ipa", "hifiasm"]:
                for suffix in [".result.fasta", ".telo.fasta"]:
                    af = f"assemblies/{asm_name}{suffix}"
                    if os.path.isfile(af) and os.path.getsize(af) > 0:
                        asm_fastas.append(af)
                        break

            if len(asm_fastas) >= 2:
                # For each protected contig, find best single-contig coverage
                # across all assembler outputs
                prot_best_cov = {}  # contig_name -> best coverage fraction
                for af in asm_fastas:
                    chk_paf = "assemblies/_chimera_check.paf"
                    cmd = (f"minimap2 -x asm20 -t {runner.threads} "
                           f"{af} {protected_fasta}")
                    result = subprocess.run(cmd, shell=True, capture_output=True,
                                            text=True)
                    with open(chk_paf, "w") as f:
                        f.write(result.stdout)

                    hits = _parse_paf_best_hits(chk_paf)
                    for qname, (cov, ident) in hits.items():
                        if ident >= 0.90:  # only count high-identity hits
                            prev = prot_best_cov.get(qname, 0.0)
                            if cov > prev:
                                prot_best_cov[qname] = cov

                # Protected contigs not well-covered by any assembler are suspect
                min_cross_cov = float(os.environ.get("CHIMERA_MIN_CROSS_COV", "0.60"))
                for n, _ in protected_recs:
                    best_cov = prot_best_cov.get(n, 0.0)
                    if best_cov < min_cross_cov:
                        mapping_suspects.add(n)
                        runner.log(f"  Chimera mapping suspect: {n} "
                                   f"(best cross-assembly cov={best_cov:.2f})")

                if os.path.isfile("assemblies/_chimera_check.paf"):
                    os.remove("assemblies/_chimera_check.paf")

        # ---- Combine both strategies ----
        # Remove contigs flagged by EITHER strategy (conservative for safety)
        chimera_remove = size_suspects | mapping_suspects
        if chimera_remove:
            filtered_recs = [(n, s) for n, s in protected_recs
                             if n not in chimera_remove]
            n_removed = n_before - len(filtered_recs)
            reasons = []
            if size_suspects:
                reasons.append(f"{len(size_suspects)} by size >{chimera_threshold:,} bp")
            if mapping_suspects:
                reasons.append(f"{len(mapping_suspects)} by cross-assembly mapping")
            runner.log_warn(f"Chimera safety: removed {n_removed} protected contigs "
                            f"({', '.join(reasons)})")
            _write_fasta(filtered_recs, protected_fasta)
        else:
            runner.log(f"Chimera safety: all {n_before} protected contigs pass "
                       f"(size threshold {chimera_threshold:,} bp, "
                       f"cross-assembly mapping OK)")

    # ====================================================================
    # BACKBONE-FIRST STRATEGY (v1.2.0)
    #
    # The backbone assembler was selected for the highest composite quality
    # score (BUSCO S, T2T count, N50, contiguity).  Its contigs carry the
    # gene content that earned it the top score — removing them to make room
    # for pool T2T contigs from other assemblers risks losing BUSCO genes
    # in non-overlapping regions.
    #
    # Strategy:
    #   1. Classify backbone contigs by telomere status FIRST.
    #   2. Backbone T2T contigs = Tier 1 (immutable, keep always).
    #   3. Backbone non-T2T contigs = Tier 2 (upgradeable).
    #   4. Pool T2T contigs that are NOT redundant to backbone Tier 1
    #      become "upgrade donors" — they can replace Tier 2 backbone
    #      contigs to add telomere completeness without losing BUSCO genes.
    #   5. Each upgrade is validated via BUSCO trial.
    #   6. Final = backbone Tier 1 + upgraded Tier 2 + remaining Tier 2.
    # ====================================================================

    taxon = getattr(runner, 'taxon', 'other')
    allow_t2t_replace = getattr(runner, 'allow_t2t_replace', False)

    # ---- 12D. Classify backbone contigs by telomere status ----
    # This MUST happen before any dedup, so we know which backbone contigs
    # are T2T (Tier 1) and which are non-T2T (Tier 2).
    runner.log_info("Classifying backbone contigs by telomere status")
    bb_t2t_ids, bb_telo_ids = _classify_contigs_telomere_status(
        backbone_fa, runner)
    bb_total = sum(1 for _ in _read_fasta_records(backbone_fa))
    bb_tier2_count = bb_total - len(bb_t2t_ids)
    runner.log(f"Backbone telomere census: {len(bb_t2t_ids)} T2T (Tier 1), "
               f"{len(bb_telo_ids)} any-telomere, {bb_total} total")
    runner.log_info(f"Tier 1 (backbone T2T, immutable): {len(bb_t2t_ids)} contigs"
                    f"{' (--allow-t2t-replace ON)' if allow_t2t_replace else ''}")
    runner.log_info(f"Tier 2 (backbone non-T2T, upgradeable): {bb_tier2_count} contigs")

    # Copy pool T2T for reference (chimera-checked version)
    if protected_fasta and os.path.isfile(protected_fasta) and \
       os.path.getsize(protected_fasta) > 0:
        shutil.copy(protected_fasta, "assemblies/protected.telomere.fa")
    else:
        with open("assemblies/protected.telomere.fa", "w") as f:
            pass

    pool_t2t_fa = "assemblies/protected.telomere.fa"

    # ---- 12D2. Find novel pool T2T contigs (not redundant to backbone) ----
    # Pool T2T contigs that cover the same chromosomes as backbone Tier 1
    # are redundant — the backbone version is preferred (preserves BUSCO).
    # Only pool T2T contigs covering regions where backbone LACKS T2T are
    # useful as upgrade donors.
    novel_pool_t2t_fa = "assemblies/novel_pool_t2t.fa"
    novel_pool_ids = set()

    if os.path.isfile(pool_t2t_fa) and os.path.getsize(pool_t2t_fa) > 0 \
       and shutil.which("minimap2"):
        runner.log_version("minimap2", "minimap2")

        # Align pool T2T against backbone (full backbone, not just Tier 1)
        pool_vs_bb_paf = "assemblies/pool_t2t_vs_backbone.paf"
        cmd = (f"minimap2 -x asm20 -t {runner.threads} "
               f"{backbone_fa} {pool_t2t_fa}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        with open(pool_vs_bb_paf, "w") as f:
            f.write(result.stdout)

        pool_hits = _parse_paf_best_hits(pool_vs_bb_paf)

        # A pool T2T contig is "redundant" if it aligns well (≥80% cov,
        # ≥90% identity) to a backbone Tier 1 (T2T) contig.
        # Pool T2T contigs that align to Tier 2 backbone contigs are
        # "upgrade donors" — they can replace the non-T2T version.
        # Pool T2T contigs that don't align to anything are novel additions.
        redundant_to_tier1 = set()
        upgrade_donors = {}  # pool_contig -> best backbone Tier 2 target

        # Parse PAF to find which backbone contig each pool contig maps to.
        # Track BOTH query coverage (pool) and target coverage (backbone).
        # For upgrade decisions, the backbone coverage matters — we must not
        # replace a large backbone contig with a smaller pool contig that only
        # covers part of it.
        pool_to_target = {}  # pool_name -> (backbone_name, qcov, tcov, ident)
        if os.path.isfile(pool_vs_bb_paf):
            with open(pool_vs_bb_paf) as f:
                for ln in f:
                    if not ln.strip() or ln.startswith("#"):
                        continue
                    p = ln.rstrip().split("\t")
                    if len(p) < 12:
                        continue
                    qname, qlen = p[0], int(p[1])
                    qs, qe = int(p[2]), int(p[3])
                    tname, tlen = p[5], int(p[6])
                    ts, te = int(p[7]), int(p[8])
                    matches, alnlen = int(p[9]), int(p[10])
                    if qlen <= 0 or alnlen <= 0:
                        continue
                    ident = matches / max(1, alnlen)
                    qcov = (qe - qs) / max(1, qlen)   # how much of pool contig aligns
                    tcov = (te - ts) / max(1, tlen)    # how much of backbone contig is covered
                    cur = pool_to_target.get(qname, ("", 0.0, 0.0, 0.0))
                    if qcov > cur[1] or (abs(qcov - cur[1]) < 1e-12 and ident > cur[3]):
                        pool_to_target[qname] = (tname, qcov, tcov, ident)

        pool_recs = list(_read_fasta_records(pool_t2t_fa))
        novel_recs = []
        for pname, pseq in pool_recs:
            if pname in pool_to_target:
                target, qcov, tcov, ident = pool_to_target[pname]
                if target in bb_t2t_ids and qcov >= 0.80 and ident >= 0.90:
                    # Pool T2T covers same chromosome as backbone T2T → redundant
                    redundant_to_tier1.add(pname)
                    runner.log(f"  Pool T2T '{pname}' redundant to backbone "
                               f"T2T '{target}' (qcov={qcov:.2f}, tcov={tcov:.2f}, "
                               f"id={ident:.2f})")
                elif target not in bb_t2t_ids and qcov >= 0.60 and \
                     tcov >= 0.80 and ident >= 0.85:
                    # Pool T2T covers a Tier 2 backbone contig:
                    # - qcov >= 0.60: pool contig aligns substantially
                    # - tcov >= 0.80: pool contig covers most of backbone contig
                    #   (prevents replacing large backbone with smaller pool contig)
                    # - ident >= 0.85: alignment quality
                    upgrade_donors[pname] = target
                    novel_recs.append((pname, pseq))
                    novel_pool_ids.add(pname)
                    runner.log(f"  Pool T2T '{pname}' can upgrade backbone "
                               f"Tier 2 '{target}' (qcov={qcov:.2f}, tcov={tcov:.2f}, "
                               f"id={ident:.2f})")
                else:
                    # Pool T2T has low coverage against backbone — may cover a
                    # new region not in backbone at all → add as novel
                    novel_recs.append((pname, pseq))
                    novel_pool_ids.add(pname)
                    if tcov < 0.80 and target not in bb_t2t_ids:
                        runner.log(f"  Pool T2T '{pname}' rejected as upgrade for "
                                   f"'{target}': only covers {tcov:.0%} of backbone contig "
                                   f"(need ≥80%); D-aware filter will evaluate")
                    else:
                        runner.log(f"  Pool T2T '{pname}' covers novel region "
                                   f"(best bb hit: {target}, qcov={qcov:.2f}, "
                                   f"tcov={tcov:.2f})")
            else:
                # No alignment to backbone at all → novel region
                novel_recs.append((pname, pseq))
                novel_pool_ids.add(pname)
                runner.log(f"  Pool T2T '{pname}' no backbone alignment → novel")

        _write_fasta(novel_recs, novel_pool_t2t_fa)
        runner.log(f"Pool T2T analysis: {len(pool_recs)} total, "
                   f"{len(redundant_to_tier1)} redundant to backbone T2T, "
                   f"{len(upgrade_donors)} upgrade donors, "
                   f"{len(novel_pool_ids) - len(upgrade_donors)} novel additions")
        # ---- 12D2b. Report un-upgraded Tier 2 backbone contigs ----
        # For Tier 2 contigs that no pool T2T covers ≥80%, log a diagnostic
        # to help identify potential chimeric backbone contigs.
        upgraded_bb = set(upgrade_donors.values())
        unupgraded_tier2 = [n for n in backbone_seqs
                            if n not in bb_t2t_ids and n not in upgraded_bb]
        if unupgraded_tier2:
            runner.log_info(f"Tier 2 backbone contigs without T2T upgrade: "
                            f"{len(unupgraded_tier2)}")
            for bbname in unupgraded_tier2:
                bblen = len(backbone_seqs[bbname])
                # Check if any pool T2T partially covers this contig
                partial_hits = [(pn, qc, tc, idt)
                                for pn, (tgt, qc, tc, idt) in pool_to_target.items()
                                if tgt == bbname and qc >= 0.50]
                if partial_hits:
                    best = max(partial_hits, key=lambda x: x[2])
                    runner.log_info(
                        f"  '{bbname}' ({bblen:,} bp): best partial T2T hit "
                        f"'{best[0]}' covers {best[2]:.0%} of backbone "
                        f"(qcov={best[1]:.0%}, id={best[3]:.1%}). "
                        f"Backbone may be chimeric or contain extra sequence. "
                        f"No replacement available — keeping backbone as-is")
                else:
                    runner.log_info(
                        f"  '{bbname}' ({bblen:,} bp): no T2T contig in pool "
                        f"covers this region — no telomere improvement possible")
    else:
        if not shutil.which("minimap2"):
            runner.log_warn("minimap2 not found; cannot compare pool T2T to backbone")
        with open(novel_pool_t2t_fa, "w") as f:
            pass

    # ---- 12D3. Backbone internal redundancy ----
    # Backbone self-dedup is DISABLED by default.  Even conservative thresholds
    # (95%/95%) can remove contigs whose 5% unique region contains BUSCO genes,
    # causing significant completeness loss (e.g. C drops from 97.4% to 92.9%
    # when 2 "redundant" contigs carry ~77 unique genes).  purge_dups at Step
    # 12H handles haplotig removal more safely using read-coverage evidence.
    #
    # Enable with SELFDEDUP_ENABLE=1 if you have known internal duplications
    # that purge_dups does not catch.
    shutil.copy(backbone_fa, "assemblies/backbone.working.fa")
    working_backbone_fa = "assemblies/backbone.working.fa"

    if os.environ.get("SELFDEDUP_ENABLE", "0") == "1":
        if taxon == "fungal":
            default_sd_cov, default_sd_id = "0.95", "0.95"
        elif taxon in ("plant", "vertebrate", "animal"):
            default_sd_cov, default_sd_id = "0.95", "0.95"
        else:
            default_sd_cov, default_sd_id = "0.95", "0.95"

        n_self_dedup = _self_dedup_non_telomeric(
            working_backbone_fa, bb_telo_ids,
            "assemblies/backbone.selfdedup.fa",
            runner.threads,
            cov_thr=float(os.environ.get("SELFDEDUP_COV", default_sd_cov)),
            id_thr=float(os.environ.get("SELFDEDUP_ID", default_sd_id)),
        )
        if n_self_dedup > 0:
            runner.log(f"Tier 2 self-dedup (taxon={taxon}): removed "
                       f"{n_self_dedup} redundant non-telomeric contigs")
            shutil.copy("assemblies/backbone.selfdedup.fa", working_backbone_fa)
        else:
            runner.log(f"Tier 2 self-dedup (taxon={taxon}): "
                       f"no redundant contigs found")
    else:
        runner.log(f"Backbone self-dedup: skipped (preserving all backbone contigs "
                   f"for BUSCO completeness; purge_dups handles haplotig removal)")

    # ---- 12E. Telomere upgrade: replace Tier 2 with pool T2T donors ----
    # Two types of upgrades:
    #   1. Direct upgrade: pool T2T contig replaces a specific Tier 2 contig
    #      (aligned as covering the same chromosome region).
    #   2. Novel addition: pool T2T contig covers a region not in backbone.
    #
    # Also check single-end telomere pool for additional rescue candidates.
    #
    # Replacement classes:
    #   upgrade_tier2_to_t2t     — Tier 2 backbone replaced by pool T2T
    #   fill_missing_end         — backbone non-telo → donor adds telomere
    #   replace_single_with_better — backbone single-telo → donor is T2T
    #   add_novel_t2t            — pool T2T covers region not in backbone

    # Collect all donor sources: pool T2T upgrades + single-end telomere pool
    single_tel_src = ""
    for candidate_file in ["single_tel_best_clean.fasta", "single_tel_best.fasta",
                           "telomere_supported_best.fasta", "single_tel_clean.fasta",
                           "single_tel.fasta"]:
        if os.path.isfile(candidate_file) and os.path.getsize(candidate_file) > 0:
            single_tel_src = candidate_file
            break

    # Pre-initialize rescue variables
    candidates = []
    donor_seqs = {}
    donor_telo_verified = set()
    donor_t2t_verified = set()
    tier1_ids = bb_t2t_ids.copy()  # Tier 1 = backbone T2T contigs

    # Read current working backbone
    backbone_seqs = dict(_read_fasta_records(working_backbone_fa))
    backbone_nodup_fa = working_backbone_fa  # alias for downstream compat

    # ---- 12E1. Build upgrade candidates from pool T2T donors ----
    if upgrade_donors and shutil.which("minimap2"):
        pool_donor_seqs = dict(_read_fasta_records(novel_pool_t2t_fa))

        for pool_name, bb_target in upgrade_donors.items():
            if bb_target in backbone_seqs and pool_name in pool_donor_seqs:
                # Determine replacement class
                if bb_target in bb_t2t_ids:
                    repl_class = "replace_protected_t2t"
                    if not allow_t2t_replace:
                        continue  # skip Tier 1 targets
                elif bb_target in bb_telo_ids:
                    repl_class = "replace_single_with_better"
                else:
                    repl_class = "upgrade_tier2_to_t2t"

                candidates.append({
                    "donor": pool_name,
                    "backbone": bb_target,
                    "replacement_class": repl_class,
                    "source": "pool_t2t",
                })

        runner.log(f"Pool T2T upgrade candidates: {len(candidates)} "
                   f"(targeting Tier 2 backbone contigs)")

        # Add donor sequences
        donor_seqs.update(pool_donor_seqs)
        donor_t2t_verified = novel_pool_ids.copy()  # pool T2T are already T2T-verified
        donor_telo_verified = novel_pool_ids.copy()

    # ---- 12E2. Also screen single-end telomere pool for rescue ----
    if single_tel_src and shutil.which("minimap2"):
        runner.log_info(f"Screening additional rescue candidates from {single_tel_src}")

        rescue_paf = "assemblies/single_tel_vs_backbone.paf"
        cmd = (f"minimap2 -x asm20 -t {runner.threads} "
               f"{working_backbone_fa} {single_tel_src}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        with open(rescue_paf, "w") as f:
            f.write(result.stdout)

        stel_donor_seqs = dict(_read_fasta_records(single_tel_src))
        all_hits = _parse_paf_rescue_hits(rescue_paf, stel_donor_seqs, backbone_seqs)
        runner.log(f"Single-end telomere rescue: {len(all_hits)} raw alignment hits")

        _write_rescue_debug_tsv(all_hits, "assemblies/single_tel.replaced.debug.tsv")

        # Structural screening
        rejected, stel_candidates = _screen_rescue_candidates(
            all_hits,
            min_ident=float(os.environ.get("RESCUE_MIN_IDENT", "0.85")),
            min_aligned_bp=int(os.environ.get("RESCUE_MIN_ALN_BP", "8000")),
            min_cov_backbone=float(os.environ.get("RESCUE_MIN_COV_BB", "0.60")),
            min_cov_donor=float(os.environ.get("RESCUE_MIN_COV_DONOR", "0.50")),
            min_ext=int(os.environ.get("RESCUE_MIN_EXT", "1000")),
        )
        runner.log(f"Structural screening: {len(rejected)} rejected, "
                   f"{len(stel_candidates)} plausible candidates")

        # Tier 1 protection: reject candidates targeting backbone T2T
        if not allow_t2t_replace:
            tier1_filtered = []
            for c in stel_candidates:
                if c["backbone"] in tier1_ids:
                    c["reject_reason"] = "target_is_tier1_backbone_t2t"
                    rejected.append(c)
                else:
                    tier1_filtered.append(c)
            n_tier1_reject = len(stel_candidates) - len(tier1_filtered)
            if n_tier1_reject > 0:
                runner.log(f"Tier 1 protection: blocked {n_tier1_reject} candidates "
                           f"targeting backbone T2T contigs")
            stel_candidates = tier1_filtered

        # Donor telomere verification for single-end candidates
        if stel_candidates:
            donor_check_fa = "assemblies/rescue_donors_check.fa"
            unique_donors = {c["donor"] for c in stel_candidates}
            _write_fasta([(n, stel_donor_seqs[n]) for n in unique_donors
                          if n in stel_donor_seqs], donor_check_fa)
            d_t2t_set, d_telo_set = _classify_contigs_telomere_status(
                donor_check_fa, runner)
            donor_telo_verified |= d_telo_set
            donor_t2t_verified |= d_t2t_set

            telo_verified = []
            for c in stel_candidates:
                if c["donor"] in donor_telo_verified:
                    telo_verified.append(c)
                else:
                    c["reject_reason"] = "donor_no_telomere_signal"
                    rejected.append(c)

            n_donor_reject = len(stel_candidates) - len(telo_verified)
            if n_donor_reject > 0:
                runner.log(f"Donor telomere verification: rejected {n_donor_reject} "
                           f"donors lacking telomere signal")
            stel_candidates = telo_verified

        # Assign replacement classes for single-end rescue
        for c in stel_candidates:
            bb_name = c["backbone"]
            d_name = c["donor"]
            if bb_name in tier1_ids:
                c["replacement_class"] = "replace_protected_t2t"
            elif bb_name not in bb_telo_ids:
                if d_name in donor_t2t_verified:
                    c["replacement_class"] = "fill_missing_end"
                else:
                    c["replacement_class"] = "replace_non_telo_backbone"
            elif bb_name in bb_telo_ids and bb_name not in bb_t2t_ids:
                c["replacement_class"] = "replace_single_with_better"
            else:
                c["replacement_class"] = "replace_non_telo_backbone"
            c["source"] = "single_tel"

        # Merge single-end candidates with pool T2T candidates
        # (pool T2T upgrades take priority — they're T2T-verified by definition)
        already_targeted = {c["backbone"] for c in candidates}
        for c in stel_candidates:
            if c["backbone"] not in already_targeted:
                candidates.append(c)
                already_targeted.add(c["backbone"])

        donor_seqs.update(stel_donor_seqs)

        # Write rejection summary
        with open("assemblies/rescue_rejection_summary.txt", "w") as f:
            f.write(f"Total hits evaluated: {len(all_hits)}\n")
            f.write(f"Obvious rejections: {len(rejected)}\n")
            f.write(f"Plausible candidates (single-end): {len(stel_candidates)}\n")
            f.write(f"Tier 1 (backbone T2T, immutable): {len(tier1_ids)}\n")
            f.write(f"Allow T2T replace: {allow_t2t_replace}\n\n")
            for r in rejected:
                f.write(f"REJECT {r.get('donor','')} → {r.get('backbone','')}: "
                        f"{r.get('reject_reason','')}\n")
    else:
        if not single_tel_src:
            runner.log_info("No single-end telomeric contigs for additional rescue")
        for p in ["assemblies/single_tel.replaced.debug.tsv",
                   "assemblies/rescue_rejection_summary.txt"]:
            with open(p, "w") as f:
                pass

    # Rank all candidates
    for i, c in enumerate(candidates):
        c["candidate_rank"] = i + 1
    runner.log(f"Total upgrade/rescue candidates: {len(candidates)}")

    _write_candidates_tsv(candidates, "assemblies/single_tel.candidates.tsv")

    # ---- 12F. BUSCO trial validation (enhanced) ----
    # Taxon-aware max rescue limits:
    if taxon == "fungal":
        default_max = "20"
    elif taxon == "plant":
        default_max = "8"
    elif taxon in ("vertebrate", "animal"):
        default_max = "10"
    else:
        default_max = "15"
    max_accepted = int(os.environ.get("STEP12_MAX_ACCEPTED", default_max))
    min_bp_ratio = float(os.environ.get("STEP12_MIN_BP_RATIO", "0.90"))

    # Taxon-aware BUSCO thresholds (C-drop, M-rise, D-rise):
    default_c_drop, default_m_rise, default_d_rise = "2.0", "0.3", "3.0"
    if taxon == "plant":
        default_c_drop, default_m_rise, default_d_rise = "4.0", "1.0", "6.0"
    elif taxon in ("vertebrate", "animal"):
        default_c_drop, default_m_rise, default_d_rise = "3.0", "0.5", "4.0"
    elif taxon == "fungal":
        default_d_rise = "2.0"  # fungi: duplication almost always artefactual

    max_busco_c_drop = float(os.environ.get("STEP12_MAX_BUSCO_C_DROP", default_c_drop))
    max_busco_m_rise = float(os.environ.get("STEP12_MAX_BUSCO_M_RISE", default_m_rise))
    max_busco_d_rise = float(os.environ.get("STEP12_MAX_BUSCO_D_RISE", default_d_rise))

    lineage = runner.busco_lineage or "ascomycota_odb10"
    busco_available = shutil.which("busco") is not None

    trial_dir = "assemblies/rescue_trials"
    os.makedirs(trial_dir, exist_ok=True)

    # Baseline = current backbone (preserves original BUSCO quality)
    baseline_recs = list(_read_fasta_records(working_backbone_fa))
    baseline_bp = sum(len(s) for _, s in baseline_recs)
    baseline_contigs = len(baseline_recs)
    baseline_busco = None

    runner.log_info(f"BUSCO trial validation: {len(candidates)} candidates, "
                    f"busco_available={busco_available}")
    if busco_available and candidates:
        runner.log_info("Computing baseline BUSCO for trial validation")
        baseline_busco = _run_busco_trial(working_backbone_fa, lineage,
                                           runner.threads, "baseline", trial_dir)
        if baseline_busco:
            runner.log(f"Baseline BUSCO (backbone): C={baseline_busco['C']:.1f}% "
                       f"D={baseline_busco.get('D', 0):.1f}% "
                       f"M={baseline_busco['M']:.1f}%")
        else:
            runner.log_warn("Baseline BUSCO parse failed; using structural checks only")

    accepted_count = 0
    trial_results = []
    current_backbone = dict(backbone_seqs)

    if not candidates:
        runner.log_info("No upgrade/rescue candidates; backbone preserved as-is")

    for cand in candidates:
        if accepted_count >= max_accepted:
            runner.log_info(f"Reached max accepted upgrades ({max_accepted}) "
                            f"for taxon={taxon}")
            break

        backbone_name = cand["backbone"]
        donor_name = cand["donor"]
        repl_class = cand.get("replacement_class", "unknown")

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

        # Size check — reject suspicious size drops
        if baseline_bp > 0 and trial_bp < baseline_bp * min_bp_ratio:
            reason = "bp_drop_suspicious"
            trial_results.append({
                "candidate_rank": cand["candidate_rank"],
                "backbone": backbone_name, "donor": donor_name,
                "replacement_class": repl_class,
                "accepted": False, "reason": reason,
                "trial_bp": trial_bp, "trial_contigs": trial_contigs,
                "baseline_busco_c": "", "trial_busco_c": "",
                "baseline_busco_d": "", "trial_busco_d": "",
                "baseline_busco_m": "", "trial_busco_m": "",
                "busco_c_delta": "", "busco_d_delta": "", "busco_m_delta": "",
            })
            runner.log(f"  Candidate {cand['candidate_rank']} "
                       f"({donor_name} → {backbone_name}): "
                       f"REJECTED (bp_drop_suspicious, {repl_class})")
            continue

        # BUSCO trial
        trial_busco = None
        c_delta = d_delta = m_delta = 0.0
        if busco_available and baseline_busco:
            trial_label = f"trial_{cand['candidate_rank']}"
            trial_busco = _run_busco_trial(trial_fa, lineage,
                                            runner.threads, trial_label, trial_dir)

        accepted = True
        reason = "accepted"

        if trial_busco and baseline_busco:
            c_delta = trial_busco["C"] - baseline_busco["C"]
            m_delta = trial_busco["M"] - baseline_busco["M"]
            d_delta = trial_busco.get("D", 0) - baseline_busco.get("D", 0)

            if c_delta < -max_busco_c_drop:
                accepted = False
                reason = f"busco_c_drop ({c_delta:.1f}% < -{max_busco_c_drop}%)"
            elif m_delta > max_busco_m_rise:
                accepted = False
                reason = f"busco_m_rise ({m_delta:.1f}% > +{max_busco_m_rise}%)"
            elif d_delta > max_busco_d_rise:
                accepted = False
                reason = f"busco_d_rise ({d_delta:.1f}% > +{max_busco_d_rise}%)"
        elif busco_available and baseline_busco and trial_busco is None:
            accepted = False
            reason = "trial_busco_failed"

        # Telomere evidence safety for single-end replacements
        if accepted and repl_class == "replace_single_with_better":
            if donor_name not in donor_t2t_verified and backbone_name in bb_telo_ids:
                if cand.get("structural_score", 0) < 7.0:
                    accepted = False
                    reason = "donor_telo_not_stronger_than_backbone"

        trial_results.append({
            "candidate_rank": cand["candidate_rank"],
            "backbone": backbone_name, "donor": donor_name,
            "replacement_class": repl_class,
            "accepted": accepted, "reason": reason,
            "trial_bp": trial_bp, "trial_contigs": trial_contigs,
            "baseline_busco_c": f"{baseline_busco['C']:.1f}" if baseline_busco else "",
            "trial_busco_c": f"{trial_busco['C']:.1f}" if trial_busco else "",
            "baseline_busco_d": f"{baseline_busco.get('D', 0):.1f}" if baseline_busco else "",
            "trial_busco_d": f"{trial_busco.get('D', 0):.1f}" if trial_busco else "",
            "baseline_busco_m": f"{baseline_busco['M']:.1f}" if baseline_busco else "",
            "trial_busco_m": f"{trial_busco['M']:.1f}" if trial_busco else "",
            "busco_c_delta": f"{c_delta:.1f}" if (trial_busco and baseline_busco) else "",
            "busco_d_delta": f"{d_delta:.1f}" if (trial_busco and baseline_busco) else "",
            "busco_m_delta": f"{m_delta:.1f}" if (trial_busco and baseline_busco) else "",
        })

        if accepted:
            accepted_count += 1
            del current_backbone[backbone_name]
            current_backbone[donor_name] = donor_seqs[donor_name]
            if trial_busco:
                baseline_busco = trial_busco
            baseline_bp = trial_bp
            runner.log(f"  Candidate {cand['candidate_rank']} "
                       f"({donor_name} → {backbone_name}): "
                       f"ACCEPTED ({repl_class})")
        else:
            runner.log(f"  Candidate {cand['candidate_rank']} "
                       f"({donor_name} → {backbone_name}): "
                       f"REJECTED ({reason}, {repl_class})")

    runner.log(f"Telomere upgrades: {accepted_count} accepted "
               f"(max={max_accepted}, taxon={taxon})")

    # Write trial summary
    trial_cols = ["candidate_rank", "backbone", "donor",
                  "replacement_class", "accepted", "reason",
                  "trial_bp", "trial_contigs",
                  "baseline_busco_c", "trial_busco_c", "busco_c_delta",
                  "baseline_busco_d", "trial_busco_d", "busco_d_delta",
                  "baseline_busco_m", "trial_busco_m", "busco_m_delta"]
    with open("assemblies/rescue_trial_summary.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(trial_cols)
        for tr in trial_results:
            w.writerow([tr.get(c, "") for c in trial_cols])

    with open("assemblies/single_tel.replaced.ids", "w") as f:
        for tr in trial_results:
            if tr["accepted"]:
                f.write(f"{tr['backbone']}\t{tr['donor']}\t"
                        f"{tr['replacement_class']}\n")

    # Build replaced_map early so novel filter can add to it.
    # Format: {donor_name: (backbone_name, replacement_class)}
    replaced_map = {}
    for tr in trial_results:
        if tr["accepted"]:
            replaced_map[tr["donor"]] = (tr["backbone"], tr["replacement_class"])

    # ---- 12F2. Add novel pool T2T contigs (D-aware duplicate filter) ----
    # Before adding a "novel" T2T contig, check if it overlaps existing backbone.
    # Three possible outcomes for each candidate:
    #   (a) Overlaps a Tier 2 (non-T2T) backbone contig → UPGRADE: replace the
    #       backbone contig with this T2T version (it's better — has telomeres).
    #   (b) Overlaps a Tier 1 (already T2T) backbone contig → REJECT: the
    #       backbone already has a T2T version, this is a pure duplicate.
    #   (c) No significant overlap → ADD as genuinely novel region.
    # Optional BUSCO D check rejects additions that increase D excessively.
    novel_additions = 0
    novel_upgrades = 0
    novel_rejected = 0
    novel_dup_cov = float(os.environ.get("NOVEL_DUP_COV", "0.80"))
    novel_dup_id = float(os.environ.get("NOVEL_DUP_ID", "0.90"))
    max_novel_d_rise = float(os.environ.get("NOVEL_MAX_D_RISE",
                                             default_d_rise))

    if os.path.isfile(novel_pool_t2t_fa) and os.path.getsize(novel_pool_t2t_fa) > 0:
        replaced_donors = {tr["donor"] for tr in trial_results if tr["accepted"]}
        novel_candidates = []
        for pname, pseq in _read_fasta_records(novel_pool_t2t_fa):
            if pname not in replaced_donors and pname not in current_backbone \
               and pname not in upgrade_donors:
                novel_candidates.append((pname, pseq))

        if novel_candidates and shutil.which("minimap2"):
            # Write current backbone for alignment
            bb_check_fa = "assemblies/backbone_for_novel_check.fa"
            _write_fasta(list(current_backbone.items()), bb_check_fa)

            for pname, pseq in novel_candidates:
                # Align novel T2T contig against current backbone
                novel_single_fa = "assemblies/novel_check_single.fa"
                _write_fasta([(pname, pseq)], novel_single_fa)

                cmd = (f"minimap2 -x asm20 -t {runner.threads} "
                       f"{bb_check_fa} {novel_single_fa}")
                result = subprocess.run(cmd, shell=True, capture_output=True,
                                        text=True)

                # Parse: which backbone contig does the novel contig align to?
                best_qcov = 0.0
                best_tcov = 0.0
                best_ident = 0.0
                best_target = ""
                for ln in (result.stdout or "").strip().split("\n"):
                    if not ln.strip():
                        continue
                    p = ln.split("\t")
                    if len(p) < 12:
                        continue
                    qlen = int(p[1])
                    qs, qe = int(p[2]), int(p[3])
                    tname, tlen = p[5], int(p[6])
                    ts, te = int(p[7]), int(p[8])
                    matches, alnlen = int(p[9]), int(p[10])
                    if qlen <= 0 or alnlen <= 0:
                        continue
                    qcov = (qe - qs) / max(1, qlen)
                    tcov = (te - ts) / max(1, tlen)
                    ident = matches / max(1, alnlen)
                    if qcov > best_qcov or (abs(qcov - best_qcov) < 1e-12
                                            and ident > best_ident):
                        best_qcov = qcov
                        best_tcov = tcov
                        best_ident = ident
                        best_target = tname

                is_dup = best_qcov >= novel_dup_cov and best_ident >= novel_dup_id

                if is_dup and best_target in current_backbone:
                    # This novel T2T overlaps an existing backbone contig.
                    # Decision depends on:
                    #   1. Whether backbone is Tier 1 (T2T) or Tier 2 (non-T2T)
                    #   2. How much of the BACKBONE the novel contig covers (tcov)
                    target_is_t2t = best_target in tier1_ids
                    novel_len = len(pseq)
                    bb_len = len(current_backbone[best_target])

                    # Minimum backbone coverage required for replacement
                    min_upgrade_tcov = float(os.environ.get(
                        "NOVEL_UPGRADE_TCOV", "0.80"))

                    if target_is_t2t:
                        # Backbone already T2T — pure duplicate, reject
                        runner.log(
                            f"  Novel T2T '{pname}' ({novel_len:,} bp) rejected: "
                            f"duplicates Tier 1 T2T '{best_target}' "
                            f"({bb_len:,} bp) at qcov={best_qcov:.0%} / "
                            f"id={best_ident:.1%}")
                        novel_rejected += 1
                        continue
                    elif best_tcov >= min_upgrade_tcov:
                        # Novel T2T covers ≥80% of Tier 2 backbone →
                        # safe to upgrade (replacement preserves most content)
                        runner.log(
                            f"  Novel T2T '{pname}' ({novel_len:,} bp) upgrades "
                            f"Tier 2 '{best_target}' ({bb_len:,} bp): "
                            f"T2T replaces non-T2T "
                            f"(qcov={best_qcov:.0%}, tcov={best_tcov:.0%}, "
                            f"id={best_ident:.1%})")
                        del current_backbone[best_target]
                        current_backbone[pname] = pseq
                        novel_upgrades += 1
                        replaced_map[pname] = (best_target,
                                               "novel_t2t_upgrades_tier2")
                        _write_fasta(list(current_backbone.items()), bb_check_fa)
                        continue
                    else:
                        # Novel T2T only partially covers Tier 2 backbone
                        # (e.g., 61% of a 1.7M contig).  Cannot upgrade
                        # (would lose backbone content) and cannot add as
                        # novel (would increase BUSCO D).  Reject entirely.
                        runner.log(
                            f"  Novel T2T '{pname}' ({novel_len:,} bp) rejected: "
                            f"partially covers Tier 2 '{best_target}' "
                            f"({bb_len:,} bp, tcov={best_tcov:.0%}). "
                            f"Too small to replace backbone "
                            f"(need ≥{min_upgrade_tcov:.0%} tcov); "
                            f"adding as novel would increase BUSCO D")
                        novel_rejected += 1
                        continue

                elif is_dup:
                    # Duplicate of something no longer in current_backbone
                    # (already replaced) — reject
                    runner.log(
                        f"  Novel T2T '{pname}' rejected: "
                        f"{best_qcov:.0%} aligns to '{best_target}' "
                        f"at {best_ident:.1%} id (duplicate)")
                    novel_rejected += 1
                    continue

                # Not a duplicate — check BUSCO D if available
                if busco_available and baseline_busco:
                    trial_recs = list(current_backbone.items()) + [(pname, pseq)]
                    trial_fa = os.path.join(trial_dir,
                                            f"novel_trial_{pname}.fa")
                    _write_fasta(trial_recs, trial_fa)
                    trial_busco = _run_busco_trial(
                        trial_fa, lineage, runner.threads,
                        f"novel_{pname}", trial_dir)
                    if trial_busco and baseline_busco:
                        d_delta = trial_busco.get("D", 0) - \
                                  baseline_busco.get("D", 0)
                        if d_delta > max_novel_d_rise:
                            runner.log(
                                f"  Novel T2T '{pname}' rejected: "
                                f"BUSCO D rises by {d_delta:+.1f}% "
                                f"(threshold: +{max_novel_d_rise}%)")
                            novel_rejected += 1
                            continue
                        runner.log(
                            f"  Novel T2T '{pname}' BUSCO D check: "
                            f"D delta {d_delta:+.1f}% (OK)")

                # Passed all filters — add as genuinely novel
                current_backbone[pname] = pseq
                novel_additions += 1
                runner.log(f"  Added novel pool T2T '{pname}' "
                           f"({len(pseq):,} bp, genuinely novel region)")
                # Update backbone FASTA for next candidate
                _write_fasta(list(current_backbone.items()), bb_check_fa)

            # Cleanup temp files
            for tmp in [bb_check_fa, "assemblies/novel_check_single.fa"]:
                if os.path.isfile(tmp):
                    os.remove(tmp)
        elif novel_candidates:
            for pname, pseq in novel_candidates:
                current_backbone[pname] = pseq
                novel_additions += 1
                runner.log(f"  Added novel pool T2T '{pname}' "
                           f"({len(pseq):,} bp) [minimap2 unavailable]")

    if novel_additions > 0 or novel_upgrades > 0 or novel_rejected > 0:
        runner.log(f"Novel pool T2T: {novel_additions} added, "
                   f"{novel_upgrades} upgraded Tier 2 backbone, "
                   f"{novel_rejected} rejected as Tier 1 duplicates")

    # Write the upgraded backbone
    _write_fasta(list(current_backbone.items()),
                 "assemblies/backbone.telomere_rescued.fa")

    # ---- 12G. Final assembly ----
    # In backbone-first mode, the assembly IS the backbone (with upgrades).
    # No separate pool + backbone merge needed.
    raw_out = "assemblies/final_merge.raw.fasta"
    shutil.copy("assemblies/backbone.telomere_rescued.fa", raw_out)

    # ---- 12G2. Post-upgrade dedup (safety net) ----
    # Only check if novel additions are redundant to existing backbone contigs.
    # Protect ALL backbone contigs (including non-telomeric Tier 2) to avoid
    # BUSCO completeness loss.  purge_dups at 12H handles true haplotig removal.
    if novel_additions > 0 and shutil.which("minimap2"):
        # Protect everything that was in the backbone (all original names)
        all_backbone_names = set(backbone_seqs.keys())
        # Also protect upgrade donors that replaced backbone contigs
        for tr in trial_results:
            if tr["accepted"]:
                all_backbone_names.add(tr["donor"])
                all_backbone_names.discard(tr["backbone"])  # replaced, not in assembly

        n_post_dedup = _self_dedup_non_telomeric(
            raw_out, all_backbone_names,
            "assemblies/final_merge.deduped.fa",
            runner.threads,
            cov_thr=0.95, id_thr=0.95,
        )
        if n_post_dedup > 0:
            shutil.copy("assemblies/final_merge.deduped.fa", raw_out)
            runner.log(f"Post-upgrade dedup: removed {n_post_dedup} redundant contigs")
        else:
            runner.log("Post-upgrade dedup: no redundancy found")

    # Keep prot variable pointing to pool T2T for downstream provenance
    prot = "assemblies/protected.telomere.fa"

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

    n_sorted, name_map = _fasta_sort_minlen_with_map(
        raw_out, f"merged_{assembler}_sort.fa", prefix="contig", minlen=500)
    runner.log(f"Sorted final merged assembly: {n_sorted} contigs >= 500 bp")

    shutil.copy(f"merged_{assembler}_sort.fa", "assemblies/final.merged.fasta")
    runner.log("Wrote assemblies/final.merged.fasta")

    # ---- Provenance GFF ----
    # Build the data structures needed for provenance annotation:
    # 1. protected_ids: contig names from the T2T protected pool
    prot_ids = set()
    if os.path.isfile(prot) and os.path.getsize(prot) > 0:
        for pname, _ in _read_fasta_records(prot):
            prot_ids.add(pname)

    # 2. replaced_map: {donor_name: (backbone_name, replacement_class)}
    # Already initialized before 12F2 with trial results + novel upgrades.
    # Merge any entries from the .ids file (in case of restart from mid-step).
    repl_ids_path = "assemblies/single_tel.replaced.ids"
    if os.path.isfile(repl_ids_path):
        with open(repl_ids_path) as f:
            for ln in f:
                parts = ln.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    replaced_map.setdefault(parts[1], (parts[0], parts[2]))
                elif len(parts) == 2:
                    replaced_map.setdefault(parts[1], (parts[0], "unknown"))

    # 3. pool_provenance_map: read provenance TSV from Step 10
    #    Maps pool contig names → (asm, orig, source_type, qm_asm1, qm_asm2, qm_regions)
    pool_provenance_map = {}
    prov_tsv = "pool_contig_provenance.tsv"
    # Also check final_results/ (cleanup step may have moved it)
    if not os.path.isfile(prov_tsv) and \
       os.path.isfile("final_results/pool_contig_provenance.tsv"):
        prov_tsv = "final_results/pool_contig_provenance.tsv"
    if os.path.isfile(prov_tsv):
        with open(prov_tsv) as f:
            header = f.readline()
            for ln in f:
                parts = ln.rstrip("\n").split("\t")
                if len(parts) >= 7 and parts[3] == "quickmerge":
                    # Parse region string: "start-end:asm:contig;..."
                    regions = []
                    if parts[6]:
                        for seg in parts[6].split(";"):
                            try:
                                rng, asm_r, contig_r = seg.split(":", 2)
                                rs, re = rng.split("-")
                                regions.append((int(rs), int(re), asm_r, contig_r))
                            except (ValueError, IndexError):
                                pass
                    pool_provenance_map[parts[0]] = (
                        parts[1], parts[2], parts[3],
                        parts[4], parts[5], regions
                    )
                elif len(parts) >= 4:
                    pool_provenance_map[parts[0]] = (parts[1], parts[2], parts[3])
                elif len(parts) >= 3:
                    pool_provenance_map[parts[0]] = (parts[1], parts[2], "unknown")
        runner.log(f"Loaded pool provenance: {len(pool_provenance_map)} entries")
    else:
        runner.log_warn("pool_contig_provenance.tsv not found; "
                        "GFF provenance will use backbone assembler for all contigs")

    # Build pool_asm_map for compatibility
    pool_asm_map = {}
    for pn, prov in pool_provenance_map.items():
        pool_asm_map[pn] = prov[0]

    _write_provenance_gff(
        "assemblies/final.merged.fasta",
        "assemblies/final.merged.provenance.gff3",
        name_map, prot_ids, replaced_map, backbone_assembler=assembler,
        pool_asm_map=pool_asm_map,
        pool_provenance_map=pool_provenance_map,
        runner=runner)

    # Also copy GFF to final_results when cleanup runs
    runner.log("Wrote assemblies/final.merged.provenance.gff3")


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

    # Always re-run BUSCO on final assembly — stale cached results from
    # previous runs cause incorrect metrics in the comparison table.
    # Remove old BUSCO output so we get fresh results.
    for old_dir in ["busco/final", "busco/run_final"]:
        if os.path.isdir(old_dir):
            shutil.rmtree(old_dir)
            runner.log(f"  Removed stale BUSCO cache: {old_dir}")

    if True:
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
    # Use shutil.copy + os.remove instead of shutil.move to avoid silent
    # failures when destination already exists from a previous run.
    for src in [
        "assemblies/final.merged.fasta",
        "assemblies/final.merged.provenance.gff3",
        "pool_contig_provenance.tsv",
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
            dst = os.path.join("final_results", os.path.basename(src))
            try:
                if os.path.isfile(dst):
                    os.remove(dst)
                shutil.copy2(src, dst)
                os.remove(src)
            except Exception as e:
                runner.log_warn(f"Could not move {src} to final_results/: {e}")

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
