#!/usr/bin/env python3
"""Join compact H4K16ac QC tables without modifying any primary output."""
from pathlib import Path
import pandas as pd

root = Path(__file__).resolve().parents[3]
base = root / "shared/biological_consensus_repaired"
qc = base / "chipseq/qc"

full = pd.read_csv(qc / "library_complexity_frip.tsv", sep="\t")
depth = pd.read_csv(qc / "library_complexity_10pct.tsv", sep="\t").rename(
    columns={"NRF":"NRF_10pct", "PBC1":"PBC1_10pct", "PBC2":"PBC2_10pct"})
insert = pd.read_csv(qc / "insert_size/insert_size_summary.tsv", sep="\t")
cross = pd.read_csv(qc / "species_cross_mapping.tsv", sep="\t")

out = full.merge(
    depth[["sample","subsampled_fragments","NRF_10pct","PBC1_10pct","PBC2_10pct"]],
    on="sample", validate="one_to_one"
).merge(insert, on="sample", validate="one_to_one").merge(
    cross, on="sample", validate="one_to_one"
)
out["qc_interpretation"] = "provisional"
out.loc[out["sample"].eq("Scr_H4K16ac_2_S0_L001"), "qc_interpretation"] = \
    "provisional; worst coordinate complexity"
out.loc[out["sample"].str.startswith("P53_"), "qc_interpretation"] = \
    "provisional; very low FRiP"

dest = base / "summary/chip_qc_summary.csv"
dest.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(dest, index=False)
print(out.to_string(index=False))
