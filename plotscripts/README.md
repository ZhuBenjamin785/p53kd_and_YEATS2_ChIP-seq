# Simple matplotlib plots

Each script reads tables from the parent `biological_consensus_repaired` directory and writes matching PNG and PDF files to `../plots/`.

Run everything with:

```bash
MPLCONFIGDIR=/tmp/matplotlib_consensus_simple \
  /home/nqp9093/.conda/envs/pydeseq2/bin/python \
  shared/biological_consensus_repaired/plotscripts/make_all_plots.py
```

`plot_chip_qc.py` deliberately shows an incomplete-QC notice until the overnight workflow has regenerated all four full-depth ChIP rows. The final postprocessing job reruns all plots automatically.
