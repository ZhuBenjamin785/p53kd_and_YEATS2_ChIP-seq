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
