# ChIP-seq summary plotting package

This folder contains the annotation, plotting, and comparison programs for:

- `p53KD` H4K16ac versus Scramble
- `YEATS2KD` H4K16ac versus Scramble
- `p53_0hr` p53 ChIP at 0 hr TNF-alpha

The R stage uses ChIPseeker to create detailed peak-annotation tables. The
Python stage uses matplotlib, seaborn, pandas, and pyBigWig to create:

- TSS-centered heatmaps and average profiles
- TSS-to-TES metagene heatmaps
- detailed peak-annotation percentage barplots
- DiffBind MA, volcano, and effect-size plots where differential results exist
- dataset master PDFs containing multiple plots per page

## Environment setup

The YAML files are portable environment specifications:

```bash
conda env create -f environment_chipseeker.yml
conda env create -f environment_pybw.yml
```

The `chipseq_chipseeker` environment runs the R annotation stage. The
`chipseq_pybw` environment runs the Python plotting stage. Conda will solve
platform-appropriate package builds on Linux, macOS, or Windows.

## Required input layout

Run from the project root containing these paths:

```text
gencode.v50.basic.annotation.gtf
bamcompare/
bamcompareY2/
fastqchip_macs3_results/bamcompare/
fastqchip_macs3_results/peaks/p53_0hr/
macs3_results_p53kd/peaks/
macs3_results_yeats2/peaks/
diffbind_results/DiffBind_all_peaks.csv
diffbind_results_yeats2/DiffBind_YEATS2_all_peaks.csv
```

Large BAM and BigWig files are intentionally not included in this package.
Copy the required data into the matching paths before running it.

## Run locally

From the project root:

```bash
conda run -n chipseq_chipseeker Rscript --vanilla chipseq_summary_package/chipseq_summary_plots.r p53KD
conda run -n chipseq_pybw python chipseq_summary_package/chipseq_summary_plots.py p53KD
```

Replace `p53KD` with `YEATS2KD` or `p53_0hr`. The output is written to:

```text
chipseq_summary_plots/<dataset>/
```

The comparison program requires completed p53KD and YEATS2KD outputs:

```bash
conda run -n chipseq_pybw python chipseq_summary_package/chipseq_compare_p53_YEATS2.py
```

## Run on Slurm

The included Slurm scripts assume the original project path and module setup.
Edit the `cd`, module, and Conda paths if the package is moved to another
cluster.

The master PDFs are named:

```text
chipseq_summary_plots/<dataset>/MASTER_<dataset>_figures.pdf
```

Signal notes: the ChIP/Input BigWig tracks are log2(ChIP/Input); the current
fastqchip p53 0-hour tracks were generated from dm6 spike-in-normalized BAMs.
The p53KD and YEATS2KD differential results come from their DiffBind outputs.
