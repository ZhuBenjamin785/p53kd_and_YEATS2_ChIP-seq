#!/usr/bin/env python3
"""Plot available H4K16ac complexity, insert-size, and FRiP QC tables."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "plots"
OUT.mkdir(exist_ok=True)
QC = BASE / "chipseq/qc"
names = {"Scr_H4K16ac_1_S0_L001":"Scr 1", "Scr_H4K16ac_2_S0_L001":"Scr 2",
         "P53_H4K16ac_1_S0_L001":"p53KD 1", "P53_H4K16ac_2_S0_L001":"p53KD 2"}
colors = ["#4C78A8", "#4C78A8", "#E45756", "#E45756"]

depth = pd.read_csv(QC / "library_complexity_10pct.tsv", sep="\t")
insert = pd.read_csv(QC / "insert_size/insert_size_summary.tsv", sep="\t")
depth["label"] = depth["sample"].map(names); insert["label"] = insert["sample"].map(names)
x = np.arange(len(depth))
fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)

w = 0.36
axes[0,0].bar(x-w/2, depth.NRF, w, color="#4C78A8", label="NRF")
axes[0,0].bar(x+w/2, depth.PBC1, w, color="#F58518", label="PBC1")
axes[0,0].set_xticks(x, depth.label); axes[0,0].set_ylim(0,1)
axes[0,0].set_title("Depth-matched library complexity (10%)"); axes[0,0].legend(frameon=False)

axes[0,1].bar(x, insert.average_insert_size, yerr=insert.insert_size_standard_deviation,
              color=colors, capsize=4)
axes[0,1].set_xticks(x, insert.label); axes[0,1].set_ylabel("Base pairs")
axes[0,1].set_title("Paired-end insert size (mean ± SD)")

axes[1,0].bar(x, 100 * insert.fraction_gt_1000bp, color=colors)
axes[1,0].set_xticks(x, insert.label); axes[1,0].set_ylabel("Fragments >1 kb (%)")
axes[1,0].set_title("Long-fragment fraction")

full = pd.read_csv(QC / "library_complexity_frip.tsv", sep="\t")
if len(full) == 4 and set(full["sample"]) == set(names):
    full["label"] = full["sample"].map(names)
    full = full.set_index("sample").loc[list(names)].reset_index()
    axes[1,1].bar(x, 100 * full.FRIP_read1_fraction, color=colors)
    axes[1,1].set_xticks(x, full.label); axes[1,1].set_ylabel("Read-1 FRiP (%)")
    axes[1,1].set_title("Full-depth FRiP")
else:
    axes[1,1].axis("off")
    axes[1,1].text(0.5, 0.5, "Full-depth QC table is incomplete.\nRerun the overnight workflow, then rerun this script.",
                   ha="center", va="center", transform=axes[1,1].transAxes)

for ax in axes.flat:
    ax.tick_params(axis="x", labelrotation=20)
fig.suptitle("H4K16ac ChIP-seq QC summary", fontweight="bold")
for suffix in ("png", "pdf"):
    fig.savefig(OUT / f"chip_qc_summary.{suffix}", dpi=300, bbox_inches="tight")
plt.close(fig)
print(OUT / "chip_qc_summary.png")
