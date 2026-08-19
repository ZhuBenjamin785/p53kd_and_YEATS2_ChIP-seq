#!/usr/bin/env python3
"""Plot RNA versus H4K16ac effects across the full eligible universes."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "plots"
OUT.mkdir(exist_ok=True)
stats = pd.read_csv(BASE / "integration/tables/full_universe_correlation_statistics.csv")

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
for ax, scope, title in zip(axes, ["all_peaks", "promoter_peaks"], ["All assigned peaks", "Promoter peaks"]):
    data = pd.read_csv(BASE / f"integration/tables/{scope}__full_eligible_gene_integration.csv")
    x = data["H4K16ac_Fold_median"].to_numpy(float)
    y = data["RNA_log2FC"].to_numpy(float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    hb = ax.hexbin(x, y, gridsize=45, mincnt=1, bins="log", cmap="viridis")
    if len(x) > 1:
        slope, intercept = np.polyfit(x, y, 1)
        grid = np.linspace(x.min(), x.max(), 100)
        ax.plot(grid, slope * grid + intercept, color="#D55E00", lw=1.8)
    ax.axhline(0, color="0.5", lw=0.8)
    ax.axvline(0, color="0.5", lw=0.8)
    pearson = stats[(stats.scope == scope) & (stats.chip_gene_summary == "median") & (stats.method == "Pearson")].iloc[0]
    spearman = stats[(stats.scope == scope) & (stats.chip_gene_summary == "median") & (stats.method == "Spearman")].iloc[0]
    ax.set_title(f"{title} (N={len(x):,})\nPearson r={pearson.correlation:.3f}; Spearman ρ={spearman.correlation:.3f}")
    ax.set_xlabel("Median H4K16ac DiffBind fold change")
    ax.set_ylabel("RNA log2 fold change")
fig.colorbar(hb, ax=axes, label="log10 genes per hexagon", shrink=0.85)
fig.suptitle("Full-universe RNA–ChIP effect-size comparison", fontweight="bold")
for suffix in ("png", "pdf"):
    fig.savefig(OUT / f"rna_chip_full_universe_scatter.{suffix}", dpi=300, bbox_inches="tight")
plt.close(fig)
print(OUT / "rna_chip_full_universe_scatter.png")
