# p53KD versus YEATS2KD overlap analysis

This workflow compares significant H4K16ac DiffBind peaks (`FDR < 0.05`) using
strict genomic overlap, direction-stratified overlap, nearest genes, promoter
neighborhoods, and ±1/±2 kb lenient overlaps.

Submit from this directory so the SLURM log paths resolve here:

```bash
cd /gpfs/projects/b1042/LauberthLab/BenFolder/p53KD_YEATS2KD_overlap_analysis
sbatch run_overlap_analysis.slurm
```

The central SLURM script calls bedtools directly from the project-local
`intervene_env`, then uses the existing `pybw` Conda environment for pandas,
seaborn, matplotlib, and matplotlib-venn. Intervene is attempted for a
supplementary peak Venn; if it cannot import, the Python Venn fallback is used.
GO Biological Process over-representation is calculated with clusterProfiler
and org.Hs.eg.db in the existing `chipseeker` environment; all GO figures are
still generated in Python.

Results are organized into `beds/`, `tables/`, `plots/`, `intervene_plots/`, and
`logs/`. No R step is needed because GENCODE v50 gene TSS annotation and the
existing ChIPseeker peak categories are already available.
