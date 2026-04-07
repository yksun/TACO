# ------------------------------------------------------------
# v0.5.6
# - MAJOR CHANGE (Steps 9, 10, 14): Replaced seqtk telo exact-motif telomere detection
#   with hybrid scoring system supporting de novo k-mer discovery, built-in motif families
#   (canonical TTAGGG, budding yeast TG1-3, Candida), and per-end composite scoring.
#   Three-tier classification: strict_t2t (both ends ≥0.25), single_tel_strong (one end ≥0.25),
#   telomere_supported (one end ≥0.08).
#
# - FIX (Step 12): Updated backbone scoring to read new telomere column names
#   ("Telomere strict T2T contigs", "Telomere single-end strong contigs") instead of
#   old names ("Telomere double-end contigs", "Telomere single-end contigs").
#
# - FIX (Step 10): Fixed awk variable scope bug where LIST was passed as a positional
#   argument instead of -v, making it unavailable in the BEGIN block and causing empty
#   telomere pool extraction.
#
# - FIX (Step 16): Updated final report to use score-based telomere classification from
#   .telomere_end_scores.tsv instead of position-based heuristics, with consistent column names.
#
# - CHANGE (Step 12): Backbone scoring now uses BUSCO S% (single-copy) instead of C%
#   (completeness) to avoid rewarding highly duplicated assemblies as backbone candidates.
#   Reduced T2T/single-end telomere weights (pool can rescue), increased contigs penalty.
#   Formula: BUSCO_S*1000 + T2T*300 + single*150 + MerquryComp*200 + MerquryQV*20
#   - contigs*30 + log10(N50)*150.
#
# - FIX (Step 10): Pool classification .list files are now written directly by Python,
#   eliminating awk TSV extraction as a potential failure point.
#
# - MAJOR CHANGE (Step 10): Reworked telomere-pool construction from seqtk telo
#   coordinate output. Contigs are now separated into biologically distinct classes:
#   strict T2T contigs, best single-end telomeric contigs, and optimized telomere-supported contigs.
#
# - MAJOR CHANGE (Step 10): Preserved the strict meaning of t2t.fasta.
#   It contains only true double-end telomere-to-telomere contigs.
#   Best single-end telomeric representatives are written to single_tel_best.fasta,
#   and the optimized combined pool is written to telomere_supported_best.fasta.
#
# - MAJOR CHANGE (Step 10): Added redundancy reduction for single-end telomeric contigs
#   using all-vs-all minimap2 clustering and longest-representative selection.
#
# - MAJOR CHANGE (Step 10): Updated protected telomere-pool priority to:
#   strict T2T > best single-end telomeric representatives > optimized telomere-supported pool.
#
# - MAJOR CHANGE (Step 12): Replaced older final merge behavior with cleaned backbone
#   refinement logic. The selected assembler output is used as the backbone assembly,
#   while redundant backbone contigs can be replaced or rescued using the optimized telomere pool.
#
# - MAJOR CHANGE (Step 12): Added support for optional automatic backbone selection modes:
#   smart scoring or legacy N50-only selection.
#
# - CHANGE (Step 12): Increased telomere contribution in smart backbone scoring to better
#   reflect telomere-end retention during final assembly refinement.
#
# - CHANGE (Step 12): Added optional Merqury integration for assembler ranking and final
#   reporting when --merqury or --merqury-db is supplied.
#
# - FIX (Step 12): Improved fungal single-end telomere rescue by relaxing terminal overhang,
#   alignment identity, and coverage thresholds.
#
# - CHANGE (Step 12): Telomere rescue now prioritizes single_tel_best_clean.fasta as the
#   preferred rescue pool before falling back to broader telomeric sets.
#
# - CHANGE (Step 14): Updated final telomere analysis to report strict T2T, single-end,
#   total telomere-supported contigs, protected telomere mode, and rescue counts.
#
# - CHANGE (Step 16): Updated final comparison reporting to include Merqury metrics,
#   telomere-pool statistics, rescue counts, selection score, selected assembler,
#   auto-selection mode, and score formula in final_results/final_result.csv.
#
# - CHANGE (Step 18): Added assembly-only mode to run assembler benchmarking/comparison
#   without final telomere-aware refinement, while still generating assemblies/assembly_info.csv
#   and optional Merqury comparison output.
#
# - FIX: Removed duplicated old parser/selection code, aligned usage text with current logic,
#   and cleaned older command fragments left from pre-0.5.x versions.
#
# - DESIGN UPDATE: TACO now follows a telomere-pool optimization and backbone-refinement
#   strategy rather than repeated structural merging, improving chromosome-end retention
#   while reducing non-telomeric redundancy.
# ------------------------------------------------------------
# v0.5.5
# - MAJOR CHANGE (Steps 9, 10, 14): Replaced seqtk telo exact-motif telomere detection
#   with hybrid scoring system supporting de novo k-mer discovery, built-in motif families
#   (canonical TTAGGG, budding yeast TG1-3, Candida), and per-end composite scoring.
#   Three-tier classification: strict_t2t (both ends ≥0.25), single_tel_strong (one end ≥0.25),
#   telomere_supported (one end ≥0.08).
#
# - FIX (Step 12): Updated backbone scoring to read new telomere column names
#   ("Telomere strict T2T contigs", "Telomere single-end strong contigs") instead of
#   old names ("Telomere double-end contigs", "Telomere single-end contigs").
#
# - FIX (Step 10): Fixed awk variable scope bug where LIST was passed as a positional
#   argument instead of -v, making it unavailable in the BEGIN block and causing empty
#   telomere pool extraction.
#
# - FIX (Step 16): Updated final report to use score-based telomere classification from
#   .telomere_end_scores.tsv instead of position-based heuristics, with consistent column names.
#
# - CHANGE (Step 12): Backbone scoring now uses BUSCO S% (single-copy) instead of C%
#   (completeness) to avoid rewarding highly duplicated assemblies as backbone candidates.
#   Reduced T2T/single-end telomere weights (pool can rescue), increased contigs penalty.
#   Formula: BUSCO_S*1000 + T2T*300 + single*150 + MerquryComp*200 + MerquryQV*20
#   - contigs*30 + log10(N50)*150.
#
# - FIX (Step 10): Pool classification .list files are now written directly by Python,
#   eliminating awk TSV extraction as a potential failure point.
#
# - MAJOR CHANGE (Step 10): Reworked telomere-pool construction from seqtk telo
#   coordinate output. Contigs are now separated into biologically distinct classes:
#   strict T2T contigs, best single-end telomeric contigs, and optimized telomere-supported contigs.
#
# - MAJOR CHANGE (Step 10): Preserved the strict meaning of t2t.fasta.
#   It contains only true double-end telomere-to-telomere contigs.
#   Best single-end telomeric representatives are written to single_tel_best.fasta,
#   and the optimized combined pool is written to telomere_supported_best.fasta.
#
# - MAJOR CHANGE (Step 10): Added redundancy reduction for single-end telomeric contigs
#   using all-vs-all minimap2 clustering and longest-representative selection.
#
# - MAJOR CHANGE (Step 10): Updated protected telomere-pool priority to:
#   strict T2T > best single-end telomeric representatives > optimized telomere-supported pool.
#
# - MAJOR CHANGE (Step 12): Replaced older final merge behavior with cleaned backbone
#   refinement logic. The selected assembler output is used as the backbone assembly,
#   while redundant backbone contigs can be replaced or rescued using the optimized telomere pool.
#
# - MAJOR CHANGE (Step 12): Added support for optional automatic backbone selection modes:
#   smart scoring or legacy N50-only selection.
#
# - CHANGE (Step 12): Increased telomere contribution in smart backbone scoring to better
#   reflect telomere-end retention during final assembly refinement.
#
# - CHANGE (Step 12): Added optional Merqury integration for assembler ranking and final
#   reporting when --merqury or --merqury-db is supplied.
#
# - FIX (Step 12): Improved fungal single-end telomere rescue by relaxing terminal overhang,
#   alignment identity, and coverage thresholds.
#
# - CHANGE (Step 12): Telomere rescue now prioritizes single_tel_best_clean.fasta as the
#   preferred rescue pool before falling back to broader telomeric sets.
#
# - CHANGE (Step 14): Updated final telomere analysis to report strict T2T, single-end,
#   total telomere-supported contigs, protected telomere mode, and rescue counts.
#
# - CHANGE (Step 16): Updated final comparison reporting to include Merqury metrics,
#   telomere-pool statistics, rescue counts, selection score, selected assembler,
#   auto-selection mode, and score formula in final_results/final_result.csv.
#
# - CHANGE (Step 18): Added assembly-only mode to run assembler benchmarking/comparison
#   without final telomere-aware refinement, while still generating assemblies/assembly_info.csv
#   and optional Merqury comparison output.
#
# - FIX: Removed duplicated old parser/selection code, aligned usage text with current logic,
#   and cleaned older command fragments left from pre-0.5.x versions.
#
# - DESIGN UPDATE: TACO now follows a telomere-pool optimization and backbone-refinement
#   strategy rather than repeated structural merging, improving chromosome-end retention
#   while reducing non-telomeric redundancy.
# ------------------------------------------------------------
# v0.5.4
# - MAJOR CHANGE (Step 10): Reworked telomere-pool construction from seqtk telo
#   coordinate output. Contigs are now separated into biologically distinct classes:
#   strict T2T contigs, best single-end telomeric contigs, and optimized telomere-supported contigs.
#
# - MAJOR CHANGE (Step 10): Preserved the strict meaning of t2t.fasta.
#   It contains only true double-end telomere-to-telomere contigs.
#   Best single-end telomeric representatives are written to single_tel_best.fasta,
#   and the optimized combined pool is written to telomere_supported_best.fasta.
#
# - MAJOR CHANGE (Step 10): Added redundancy reduction for single-end telomeric contigs
#   using all-vs-all minimap2 clustering and longest-representative selection.
#
# - MAJOR CHANGE (Step 10): Updated protected telomere-pool priority to:
#   strict T2T > best single-end telomeric representatives > optimized telomere-supported pool.
#
# - MAJOR CHANGE (Step 12): Replaced older final merge behavior with cleaned backbone
#   refinement logic. The selected assembler output is used as the backbone assembly,
#   while redundant backbone contigs can be replaced or rescued using the optimized telomere pool.
#
# - MAJOR CHANGE (Step 12): Added support for optional automatic backbone selection modes:
#   smart scoring or legacy N50-only selection.
#
# - CHANGE (Step 12): Increased telomere contribution in smart backbone scoring to better
#   reflect telomere-end retention during final assembly refinement.
#
# - CHANGE (Step 12): Added optional Merqury integration for assembler ranking and final
#   reporting when --merqury or --merqury-db is supplied.
#
# - FIX (Step 12): Improved fungal single-end telomere rescue by relaxing terminal overhang,
#   alignment identity, and coverage thresholds.
#
# - CHANGE (Step 12): Telomere rescue now prioritizes single_tel_best_clean.fasta as the
#   preferred rescue pool before falling back to broader telomeric sets.
#
# - CHANGE (Step 14): Updated final telomere analysis to report strict T2T, single-end,
#   total telomere-supported contigs, protected telomere mode, and rescue counts.
#
# - CHANGE (Step 16): Updated final comparison reporting to include Merqury metrics,
#   telomere-pool statistics, rescue counts, selection score, selected assembler,
#   auto-selection mode, and score formula in final_results/final_result.csv.
#
# - CHANGE (Step 18): Added assembly-only mode to run assembler benchmarking/comparison
#   without final telomere-aware refinement, while still generating assemblies/assembly_info.csv
#   and optional Merqury comparison output.
#
# - FIX: Removed duplicated old parser/selection code, aligned usage text with current logic,
#   and cleaned older command fragments left from pre-0.5.x versions.
#
# - DESIGN UPDATE: TACO now follows a telomere-pool optimization and backbone-refinement
#   strategy rather than repeated structural merging, improving chromosome-end retention
#   while reducing non-telomeric redundancy.
# ------------------------------------------------------------
# v0.5.3
# - MAJOR CHANGE (Step 10): Reworked telomere-pool construction from seqtk telo
#   coordinate output. Contigs are now separated into biologically distinct classes:
#   strict T2T contigs, best single-end telomeric contigs, and optimized telomere-supported contigs.
#
# - MAJOR CHANGE (Step 10): Preserved the strict meaning of t2t.fasta.
#   It contains only true double-end telomere-to-telomere contigs.
#   Best single-end telomeric representatives are written to single_tel_best.fasta,
#   and the optimized combined pool is written to telomere_supported_best.fasta.
#
# - MAJOR CHANGE (Step 10): Added redundancy reduction for single-end telomeric contigs
#   using all-vs-all minimap2 clustering and longest-representative selection.
#
# - MAJOR CHANGE (Step 10): Updated protected telomere-pool priority to:
#   strict T2T > best single-end telomeric representatives > optimized telomere-supported pool.
#
# - MAJOR CHANGE (Step 12): Replaced older final merge behavior with cleaned backbone
#   refinement logic. The selected assembler output is used as the backbone assembly,
#   while redundant backbone contigs can be replaced or rescued using the optimized telomere pool.
#
# - MAJOR CHANGE (Step 12): Added support for optional automatic backbone selection modes:
#   smart scoring or legacy N50-only selection.
#
# - CHANGE (Step 12): Increased telomere contribution in smart backbone scoring to better
#   reflect telomere-end retention during final assembly refinement.
#
# - CHANGE (Step 12): Added optional Merqury integration for assembler ranking and final
#   reporting when --merqury or --merqury-db is supplied.
#
# - FIX (Step 12): Improved fungal single-end telomere rescue by relaxing terminal overhang,
#   alignment identity, and coverage thresholds.
#
# - CHANGE (Step 12): Telomere rescue now prioritizes single_tel_best_clean.fasta as the
#   preferred rescue pool before falling back to broader telomeric sets.
#
# - CHANGE (Step 14): Updated final telomere analysis to report strict T2T, single-end,
#   total telomere-supported contigs, protected telomere mode, and rescue counts.
#
# - CHANGE (Step 16): Updated final comparison reporting to include Merqury metrics,
#   telomere-pool statistics, rescue counts, selection score, selected assembler,
#   auto-selection mode, and score formula in final_results/final_result.csv.
#
# - CHANGE (Step 18): Added assembly-only mode to run assembler benchmarking/comparison
#   without final telomere-aware refinement, while still generating assemblies/assembly_info.csv
#   and optional Merqury comparison output.
#
# - FIX: Removed duplicated old parser/selection code, aligned usage text with current logic,
#   and cleaned older command fragments left from pre-0.5.x versions.
#
# - DESIGN UPDATE: TACO now follows a telomere-pool optimization and backbone-refinement
#   strategy rather than repeated structural merging, improving chromosome-end retention
#   while reducing non-telomeric redundancy.
# ------------------------------------------------------------
# v0.5.2
# - CHANGE (Step 10): Reworked telomere-pool construction to prioritize strict T2T contigs
#   and select best single-end telomeric representatives from redundant alternatives.
# - CHANGE (Step 10): Added redundancy reduction for single-end telomeric contigs using
#   all-vs-all minimap2 clustering and longest-representative selection.
# - CHANGE (Step 10): Replaced broad telomere-supported fallback with an optimized telomere pool:
#   strict T2T > best single-end telomeric contigs > optimized telomere-supported contigs.
# - CHANGE (Step 10): Added telomere_cluster_summary.tsv and updated telomere support summary
#   reporting for representative telomeric contig selection.
# - CHANGE (Step 12): Updated backbone refinement to use the optimized Step 10 telomere pool
#   rather than an undifferentiated telomere-supported set.
# - CHANGE (Step 12): Increased telomere contribution in smart scoring to better reflect
#   telomere-end retention during backbone selection.
# - FIX (Step 12): Improved fungal single-end telomere rescue by relaxing thresholds for
#   terminal overhang, alignment identity, and coverage.
# - CHANGE (Step 12): Telomere rescue now prioritizes single_tel_best_clean.fasta as the
#   preferred rescue pool before falling back to broader telomeric sets.
# - CHANGE: TACO now more explicitly optimizes for the best telomere-end pool before final
#   backbone refinement, improving chromosome-end retention in fungal assemblies.
# ------------------------------------------------------------
# v0.5.1
# - FIX (Step 12): Unified argument parsing for --choose and --auto-mode to prevent
#   auto-mode from being ignored during final backbone selection.
# - CHANGE (Step 12): Added Merqury pre-selection evaluation for all available assembler
#   outputs when merqury.sh and a .meryl database are present.
# - CHANGE (Step 12): Added Merqury-derived metrics (QV and completeness) to
#   assemblies/assembly.merqury.csv and integrated them into assemblies/assembly_info.csv.
# - CHANGE (Step 12): Updated automatic backbone scoring to include BUSCO completeness,
#   telomere support, Merqury completeness, Merqury QV, contig count, and N50.
# - CHANGE (Step 12): Added selection_debug.tsv and selection_decision.txt for transparent
#   scoring, reproducibility, and reporting of backbone selection.
# - FIX (Step 12): Improved telomere rescue logic for single-end telomeric contigs by using
#   terminal overhang rather than whole-contig length differences.
# - FIX (Step 12): Added post-rescue deduplication against protected telomere contigs to
#   reduce redundancy before final assembly construction.
# - CHANGE (Step 16): Added final Merqury evaluation for assemblies/final.merged.fasta and
#   wrote assemblies/merged.merqury.csv for final assembly comparison.
# - CHANGE (Step 16): Renamed the final comparison table output to
#   final_results/final_result.csv and expanded it to include Merqury metrics,
#   selection score, selected assembler, auto-selection mode, and scoring formula.
# - FIX (Step 17): Updated cleanup to preserve new final outputs, including
#   merged.merqury.csv, selection_debug.tsv, selection_decision.txt, and final_result.csv.
# - CHANGE: TACO now supports telomere-aware backbone refinement with optional
#   k-mer-based assembly quality scoring through Merqury.
# ------------------------------------------------------------
# v0.5.0
# - MAJOR CHANGE (Step 10): Reworked telomere-support classification from seqtk telo
#   coordinate output. Contigs are now separated into three biologically distinct classes:
#   strict T2T contigs, single-end telomeric contigs, and combined telomere-supported contigs.
#
# - MAJOR CHANGE (Step 10): Preserved the strict meaning of t2t.fasta.
#   It now contains only true double-end telomere-to-telomere contigs.
#   Single-end telomeric contigs are written to single_tel.fasta, and the union of both
#   classes is written to telomere_supported.fasta.
#
# - MAJOR CHANGE (Step 10): Added telomere-supported fallback protection for downstream
#   refinement. The pipeline now writes protected_telomere_contigs.fasta and
#   protected_telomere_mode.txt so downstream logic can prioritize strict T2T contigs first
#   and fall back to telomere-supported contigs when no strict T2T contigs are available.
#
# - MAJOR CHANGE (Step 12): Replaced merge_wrapper-based final merge logic with
#   backbone refinement logic. The selected assembler output is now used as the primary
#   backbone assembly, and redundant backbone contigs are replaced by protected
#   telomere-supported contigs rather than merged again.
#
# - MAJOR CHANGE (Step 12): Updated final assembly refinement to use
#   protected_telomere_contigs.fasta instead of relying only on t2t_clean.fasta.
#   Strict T2T contigs are prioritized first; if absent, telomere-supported contigs are used
#   as the protected replacement set.
#
# - CHANGE (Step 12): Improved redundancy filtering between protected telomere contigs
#   and selected-assembly backbone contigs using minimap2, reducing non-telomeric duplication
#   while preserving telomere-supported sequence.
#
# - CHANGE (Step 12): Renamed Step 12 from "Final merge" to
#   "Final assembly refinement with telomere-supported contig replacement"
#   to better reflect the updated pipeline design.
#
# - CHANGE (Step 14): Updated final telomere analysis to use the same coordinate-window
#   logic as Step 10, ensuring consistent reporting of strict T2T and single-end telomeric
#   contigs between intermediate and final outputs.
#
# - CHANGE (Step 16): Updated final telomere metric recomputation to match the new
#   coordinate-window logic, preventing discrepancies between Step 14 results and the final
#   comparison table.
#
# - CHANGE: Added telomere support summary reporting, including counts of strict T2T,
#   single-end telomeric, and total telomere-supported contigs, for clearer downstream
#   interpretation and manuscript reporting.
#
# - CHANGE: Improved robustness when no strict T2T contigs are present, preventing empty
#   FASTA outputs and allowing the pipeline to continue using telomere-supported fallback logic.
#
# - DESIGN UPDATE: The final pipeline now follows a protection/replacement strategy rather
#   than repeated structural merging: telomere-supported contigs are preserved first, and the
#   best selected assembly is used as a nonredundant backbone.
# ------------------------------------------------------------
# v0.30 (2026-03-19)
# - CHANGE: integrated the former external helper t2t_list.sh as an internal function.
# - CHANGE: extract_contig_t2t_V3.sh and extract_contig_T_V3.sh are not called by this version, so no external helper scripts are needed for those steps.
# - FIX (Step 12): Pass CLI args to Python heredocs correctly (no post-heredoc args).
#                  Fixes IndexError and "Permission denied" after merge.
# - FIX: Robust project name parsing for .fastq.gz.
# - Tweak: expose PROTECT_COV/PROTECT_ID env vars for redundancy filter (default 0.95/0.95).
# - FIX (Step 12): Pre-filter 'others' contigs against protected T2T via minimap2
#                  (drop if best hit has cov>=0.95 & id>=0.95), then run redundans
#                  ONLY on the non-redundant remainder. Recombine with protected T2T
#                  (unchanged). Tunable via PROTECT_COV and PROTECT_ID.
# - Tweak: Correct minor typos in --choose prompt.
# - FIX (Step 16): Ensure 'final' telomere metrics are computed reliably by recomputing
#                   from assemblies/final.telo.list during final_result merge, and by
#                   normalizing whitespace/CRs. Prevents 'Telomere double-end contigs'
#                   from showing 0 when list indicates >0.
# - FIX (Step 14): Compute 'Telomere double-end contigs' by aggregating left/right telomeric hits
#                  per contig (s[c]>0 && e[c]>0), matching Step 9 logic. Also strip CRs before parsing.
# - FIX (Step 12): Replaced awk-based assembly_info builder with a Python-based
#   build_assembly_info_v2() that strips CRs, normalizes headers, and merges BUSCO/QUAST/TELO
#   reliably across awk variants. Called at the start of Step 12.
# - FIX (Step 12): Make AWK CR stripping portable and robust:
#     • Use gsub("\\r","",$i) instead of regex literal.
#     • Pipe each CSV through tr -d '\\r' before AWK.
#   This eliminates the 'gsub(/' parse error on some AWK builds.
# - FIX (Step 12): Always rebuild and overwrite assemblies/assembly_info.csv by merging BUSCO/QUAST/TELO CSVs if present;
#                  corrected awk CR-strip and stabilized header/rows; print matrix before prompt/--choose.
# - FIX (Step 7): Removed stray Bash inside Python heredoc; log Python version from Bash before invoking the temp script.
# - FIX (Steps 9 & 14): Telomere double-end contigs counted when both ends have telomeric signal (start==0 and end==len).
# - Version bump from 0.2 → 0.2.1.
# - Step 12 note: t2t_clean.fasta contigs are preserved; redundans runs only on non-T2T contigs.
#   (--choose still supported; interactive prompt retained.)
# ------------------------------------------------------------
# v0.2.7 (2025-11-14)
# - CHANGE (Step 12): Normalize others.fa headers, run redundans with minimap2-based reduction,
#                     and use corrected others.filtered.fa writer; log contig reduction stats.
# ------------------------------------------------------------
# v0.2.6.8 (2025-11-12)
# - CHANGE (Step 12): Auto-select assembler with highest N50 from assemblies/assembly_info.csv
#                     when --choose is not provided; explicit --choose still overrides.
# ------------------------------------------------------------
# v0.2.6.7 (2025-11-12)
# - FIX (Step 12): Guard against missing/empty t2t_clean.fasta. If absent,
#   skip merge_wrapper and pass through the chosen assembler result directly.
#   Prevents FileNotFoundError at merge_wrapper.py when no T2T contigs exist.
# ------------------------------------------------------------
# ------------------------------------------------------------
# v0.2.6.6 (2025-11-12)
# - FIX (Step 10): If one or more assemblies are missing, proceed with the
#   available FASTAs in ./assemblies for the merge. For a single available
#   assembly, skip pairwise merge and copy it directly to allmerged_telo.fasta.
#   Continue to generate t2t.fasta and t2t_clean.fasta from what is available.
#   (Avoids failures when some assemblers are absent.)
# ------------------------------------------------------------
# ------------------------------------------------------------