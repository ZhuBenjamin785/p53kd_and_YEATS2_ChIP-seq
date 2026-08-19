#!/usr/bin/env python3
"""Run every simple matplotlib visualization script in this directory."""
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
for name in ["plot_fisher_tables.py", "plot_rna_chip_scatter.py",
             "plot_enrichment_summary.py", "plot_chip_qc.py"]:
    subprocess.run([sys.executable, str(HERE / name)], check=True)
print(f"Plots written to {HERE.parent / 'plots'}")
