#!/usr/bin/env python3
"""Plot the two corrected Fisher contingency tables."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "plots"
OUT.mkdir(exist_ok=True)
df = pd.read_csv(BASE / "integration/tables/fisher_exact_results.csv")

fig, axes = plt.subplots(1, len(df), figsize=(10, 4.3), constrained_layout=True)
for ax, (_, row) in zip(np.atleast_1d(axes), df.iterrows()):
    table = np.array([[row.chip_and_rna, row.chip_and_not_rna],
                      [row.no_chip_and_rna, row.no_chip_and_not_rna]])
    image = ax.imshow(np.log10(table + 1), cmap="Blues", vmin=0)
    for (i, j), value in np.ndenumerate(table):
        ax.text(j, i, f"{value:,}", ha="center", va="center",
                color="white" if image.norm(np.log10(value + 1)) > 0.55 else "black",
                fontsize=12, fontweight="bold")
    ax.set_xticks([0, 1], ["RNA set", "Not RNA set"])
    ax.set_yticks([0, 1], ["ChIP set", "Not ChIP set"])
    ax.set_title(f"{row['test']}\nOR={row.odds_ratio:.3g}; one-sided p={row.fisher_p_value:.3g}\nBH FDR={row.BH_FDR_two_directional_tests:.3g}", fontsize=10)
fig.suptitle("Directional Fisher tests in the 2,100-gene matched universe", fontweight="bold")
for suffix in ("png", "pdf"):
    fig.savefig(OUT / f"fisher_contingency_tables.{suffix}", dpi=300, bbox_inches="tight")
plt.close(fig)
print(OUT / "fisher_contingency_tables.png")
