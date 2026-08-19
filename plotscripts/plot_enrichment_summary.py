#!/usr/bin/env python3
"""Plot compact pathway-count summaries for corrected FEA, GSEA, and ORA."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "plots"
OUT.mkdir(exist_ok=True)
colors = {"GO_BP":"#4C78A8", "KEGG":"#F58518", "Reactome":"#54A24B"}

fea = pd.read_csv(BASE / "fea/directional/tables/fea_run_summary.csv")
gsea = pd.read_csv(BASE / "gsea/all_eligible/tables/gsea_run_summary.csv")
ora = pd.read_csv(BASE / "ora/all_peaks/tables/ora_run_summary.csv")
fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)

# Directional RNA FEA
x = np.arange(2); width = 0.23
for i, db in enumerate(colors):
    vals = [int(fea[(fea.Direction == direction) & (fea.Database == db)].Significant_pathways.iloc[0]) for direction in ["RNA_up", "RNA_down"]]
    axes[0].bar(x + (i-1)*width, vals, width, label=db, color=colors[db])
axes[0].set_xticks(x, ["RNA up", "RNA down"]); axes[0].set_ylabel("Significant pathways (FDR < 0.05)")
axes[0].set_title("Directional RNA FEA"); axes[0].legend(frameon=False, fontsize=8)

# GSEA positive/negative directions
x = np.arange(len(gsea))
axes[1].bar(x, gsea.Positive_NES_significant, color="#D55E00", label="Positive NES")
axes[1].bar(x, gsea.Negative_NES_significant, bottom=gsea.Positive_NES_significant, color="#0072B2", label="Negative NES")
axes[1].set_xticks(x, gsea.Database); axes[1].set_title("All-gene GSEA")
axes[1].set_ylabel("Significant pathways (FDR < 0.05)"); axes[1].legend(frameon=False, fontsize=8)

# Integrated all-peak ORA
sets = ["loss_down", "loss_up", "gain_up", "gain_down"]
labels = ["Loss / down", "Loss / up", "Gain / up", "Gain / down"]
x = np.arange(len(sets))
bottom = np.zeros(len(sets))
for db in colors:
    vals = np.array([int(ora[(ora.Gene_set == gene_set) & (ora.Database == db)].Significant_pathways.iloc[0]) for gene_set in sets])
    axes[2].bar(x, vals, bottom=bottom, color=colors[db], label=db)
    bottom += vals
axes[2].set_xticks(x, labels, rotation=25, ha="right"); axes[2].set_title("Integrated ORA: all peaks")
axes[2].set_ylabel("Significant pathways (FDR < 0.05)"); axes[2].legend(frameon=False, fontsize=8)

fig.suptitle("Corrected functional-enrichment summary", fontweight="bold")
for suffix in ("png", "pdf"):
    fig.savefig(OUT / f"enrichment_summary.{suffix}", dpi=300, bbox_inches="tight")
plt.close(fig)
print(OUT / "enrichment_summary.png")
