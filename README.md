<<<<<<< HEAD
# BenFolder summary

This directory is a curated copy of plots, result tables, annotations, QC
outputs, and analysis notes collected recursively from the project. Original
files remain in their source locations.

## Sections

- `plots/` — figures and plot PDFs, grouped by dataset.
- `tables/` — analysis result tables and signal summaries.
- `annotations/` — ChIPseeker, peak, and gene annotation outputs.
- `overlap_analysis/` — the p53KD versus YEATS2KD strict-to-lenient comparison,
  including its plots, BED/TSV/CSV results, and interpretation files.
- `qc/` — MultiQC, FastQC, and other quality-control outputs.
- `metadata/` — sample sheets, manifests, normalization factors, and QC summary
  spreadsheets.
- `reports/` — figure notes, status files, and analysis documentation.

Dataset names are retained as subdirectories where they were identifiable.
Files with identical names from different source folders were renamed with a
source-parent suffix to prevent overwriting.

Current inventory: 265 files, approximately 499 MB. Redundant byte-identical
copies created when the overlap-analysis summary was collected recursively
were removed; the canonical copy is retained under `overlap_analysis/`.
=======
# Project scripts

Scripts are grouped by pipeline stage:

- `alignment/` — alignment programs.
- `endtoend/` — end-to-end workflows.
- `preprocessing/` — FASTQ, BAM, consensus, spike-in, and count preparation.
- `qc/` — trimming, post-alignment QC, shift checks, and signal extraction.
- `peak_annotation/` — MACS peak workflows and ChIPseeker annotation.
- `differential_expression/` — DESeq2 and DiffBind workflows.
- `chipseq_summary/` — ChIP-seq summary plots; `package/` contains package copies.
- `visualization/` — metaplots, comparisons, and plotting utilities.
- `utilities/` — bigWig comparison and local serving helpers.
- `overlap_analysis/` — p53/YEATS2 overlap and GO-enrichment analysis.

Batch scripts change to the project directory before running, so data and result paths remain relative to the project root.
>>>>>>> 41782f5 (all scripts)
