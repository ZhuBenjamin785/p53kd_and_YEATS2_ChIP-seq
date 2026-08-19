#!/usr/bin/env python3
"""Run every supplementary plot script; missing optional inputs are reported and skipped."""
from pathlib import Path
import subprocess,sys
HERE=Path(__file__).resolve().parent
SCRIPTS=["plot_rna_pca.py","plot_rna_replicate_correlation.py","plot_rna_ma_volcano.py","plot_rna_top_gene_heatmap.py","plot_chip_qc_supplementary.py","plot_chip_replicate_correlation.py","plot_mapq30_crossmapping.py","plot_corrected_chip_pca_correlation.py","plot_diffbind_primary_vs_mapq30.py","plot_rna_chip_scatter.py","plot_fisher_odds_ratio.py","plot_integrated_category_counts.py","plot_genome_browser_panels.py"]
for s in SCRIPTS:
    print(f"=== {s} ===",flush=True); subprocess.run([sys.executable,str(HERE/s)],check=False)
