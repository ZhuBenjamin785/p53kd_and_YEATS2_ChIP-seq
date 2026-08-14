#!/usr/bin/env python3
"""Analyze overlap between MOF-loss and H4K16ac-loss DiffBind peaks."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import pandas as pd
import seaborn as sns


BED_COLUMNS = ["chromosome", "start", "end", "peak_id"]
PAIR_COLUMNS = [
    "MOF_chromosome", "MOF_start", "MOF_end", "MOF_peak_id",
    "H4K16ac_chromosome", "H4K16ac_start", "H4K16ac_end", "H4K16ac_peak_id",
]


def normalize_chromosome(value):
    value = str(value).strip()
    return value if value.startswith("chr") else "chr" + value


def read_diffbind(path):
    table = pd.read_csv(path)
    required = {"seqnames", "start", "end", "Fold", "FDR"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
    for column in ["start", "end", "Fold", "FDR"]:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    table = table.dropna(subset=["seqnames", "start", "end", "Fold", "FDR"]).copy()
    table["seqnames"] = table["seqnames"].map(normalize_chromosome)
    return table


def loss_table(path, prefix):
    table = read_diffbind(path)
    loss = table[(table["FDR"] < 0.05) & (table["Fold"] < 0)].copy()
    loss = loss.sort_values(["seqnames", "start", "end"]).reset_index(drop=True)
    loss["peak_id"] = [f"{prefix}_loss_{index:06d}" for index in range(1, len(loss) + 1)]
    return loss


def write_unsorted_bed(table, path):
    bed = pd.DataFrame({
        "chromosome": table["seqnames"],
        "start": table["start"].astype(int) - 1,
        "end": table["end"].astype(int),
        "peak_id": table["peak_id"],
    })
    bed.to_csv(path, sep="\t", header=False, index=False)


def prepare(args):
    output = Path(args.output)
    beds = output / "beds"
    tables = output / "tables"
    plots = output / "plots"
    for directory in [beds, tables, plots]:
        directory.mkdir(parents=True, exist_ok=True)

    mof = loss_table(args.mof, "MOF")
    h4 = loss_table(args.h4k16ac, "H4K16ac")
    if mof.empty or h4.empty:
        raise ValueError(f"Loss filtering produced MOF={len(mof)} and H4K16ac={len(h4)} peaks")

    write_unsorted_bed(mof, beds / "MOF_loss_unsorted.bed")
    write_unsorted_bed(h4, beds / "H4K16ac_loss_unsorted.bed")
    mof.to_csv(tables / "MOF_loss_DiffBind.csv", index=False)
    h4.to_csv(tables / "H4K16ac_loss_DiffBind.csv", index=False)
    print(f"Prepared {len(mof)} MOF-loss and {len(h4)} H4K16ac-loss peaks", flush=True)


def read_bed(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=BED_COLUMNS)
    return pd.read_csv(path, sep="\t", names=BED_COLUMNS, dtype={"peak_id": str})


def read_pairs(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=PAIR_COLUMNS)
    return pd.read_csv(path, sep="\t", names=PAIR_COLUMNS, dtype=str).assign(
        MOF_start=lambda x: pd.to_numeric(x.MOF_start),
        MOF_end=lambda x: pd.to_numeric(x.MOF_end),
        H4K16ac_start=lambda x: pd.to_numeric(x.H4K16ac_start),
        H4K16ac_end=lambda x: pd.to_numeric(x.H4K16ac_end),
    )


def build_shared(args):
    output = Path(args.output)
    pairs = read_pairs(output / "tables" / "shared_loss_pairs.tsv")
    pairs = pairs.drop_duplicates(["MOF_peak_id", "H4K16ac_peak_id"]).reset_index(drop=True)
    shared = pd.DataFrame({
        "chromosome": pairs["MOF_chromosome"],
        "start": np.maximum(pairs["MOF_start"], pairs["H4K16ac_start"]).astype(int),
        "end": np.minimum(pairs["MOF_end"], pairs["H4K16ac_end"]).astype(int),
        "peak_id": [f"shared_loss_pair_{index:06d}" for index in range(1, len(pairs) + 1)],
    })
    pairs["shared_peak_id"] = shared["peak_id"]
    pairs.to_csv(output / "tables" / "shared_loss_pairs_with_ids.tsv", sep="\t", index=False)
    shared.to_csv(output / "beds" / "shared_loss_intervals.bed", sep="\t", header=False, index=False)
    print(f"Built {len(shared)} shared intervals for ChIPseeker annotation", flush=True)


def read_chrom_sizes(path):
    sizes = pd.read_csv(path, sep="\t", names=["chromosome", "length"])
    sizes["chromosome"] = sizes["chromosome"].map(normalize_chromosome)
    return dict(zip(sizes.chromosome, sizes.length.astype(int)))


def overlap_flags(left, right):
    flags = np.zeros(len(left), dtype=bool)
    right_groups = {
        chrom: values[["start", "end"]].to_numpy(dtype=int)
        for chrom, values in right.groupby("chromosome")
    }
    for index, row in enumerate(left.itertuples(index=False)):
        intervals = right_groups.get(row.chromosome)
        if intervals is not None:
            flags[index] = np.any((intervals[:, 0] < row.end) & (intervals[:, 1] > row.start))
    return flags


def permutation_test(mof, h4, chrom_sizes, iterations, seed):
    valid_mof = mof[mof.chromosome.isin(chrom_sizes)].reset_index(drop=True)
    valid_h4 = h4[h4.chromosome.isin(chrom_sizes)].reset_index(drop=True)
    observed = int(overlap_flags(valid_mof, valid_h4).sum())
    rng = np.random.default_rng(seed)
    null = np.zeros(iterations, dtype=int)
    widths = (valid_mof.end - valid_mof.start).to_numpy(dtype=int)

    for iteration in range(iterations):
        randomized = valid_mof.copy()
        starts = []
        for chrom, width in zip(randomized.chromosome, widths):
            maximum = max(0, chrom_sizes[chrom] - width)
            starts.append(int(rng.integers(0, maximum + 1)))
        randomized["start"] = starts
        randomized["end"] = randomized["start"].to_numpy() + widths
        null[iteration] = int(overlap_flags(randomized, valid_h4).sum())

    expected = float(null.mean())
    empirical_p = float((1 + np.sum(null >= observed)) / (iterations + 1))
    enrichment = float(observed / expected) if expected > 0 else np.nan
    z_score = float((observed - expected) / null.std(ddof=1)) if null.std(ddof=1) > 0 else np.nan
    return {
        "iterations": iterations,
        "seed": seed,
        "MOF_peaks_tested": len(valid_mof),
        "H4K16ac_peaks_tested": len(valid_h4),
        "observed_overlapping_MOF_peaks": observed,
        "mean_random_overlapping_MOF_peaks": expected,
        "enrichment_fold": enrichment,
        "z_score": z_score,
        "empirical_p_value": empirical_p,
    }, null


def make_plots(summary, null, output):
    sns.set_theme(style="whitegrid", context="talk")
    plots = output / "plots"

    fig, ax = plt.subplots(figsize=(8, 6))
    data = pd.DataFrame({
        "Peak set": ["MOF loss", "H4K16ac loss"],
        "Overlapping (%)": [
            summary["MOF_overlap_percentage"],
            summary["H4K16ac_overlap_percentage"],
        ],
    })
    sns.barplot(data=data, x="Peak set", y="Overlapping (%)", hue="Peak set",
                palette=["#B279A2", "#E45756"], legend=False, ax=ax)
    ax.set_ylim(0, max(100, data["Overlapping (%)"].max() * 1.15))
    ax.set_title("Reciprocal overlap of significant loss peaks")
    for patch, value in zip(ax.patches, data["Overlapping (%)"]):
        ax.text(patch.get_x() + patch.get_width() / 2, patch.get_height() + 1,
                f"{value:.1f}%", ha="center", va="bottom", fontweight="bold")
    fig.tight_layout()
    fig.savefig(plots / "loss_peak_overlap_fractions.pdf")
    fig.savefig(plots / "loss_peak_overlap_fractions.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_aspect("equal")
    ax.add_patch(Circle((0.42, 0.5), 0.30, color="#B279A2", alpha=0.55))
    ax.add_patch(Circle((0.64, 0.5), 0.30, color="#E45756", alpha=0.55))
    ax.text(0.26, 0.50, str(summary["MOF_only_peaks"]), ha="center", va="center", fontsize=18)
    ax.text(0.80, 0.50, str(summary["H4K16ac_only_peaks"]), ha="center", va="center", fontsize=18)
    center = (f"{summary['overlapping_MOF_peaks']} MOF peaks\n"
              f"{summary['overlapping_H4K16ac_peaks']} H4K16ac peaks\n"
              f"{summary['overlap_pairs']} pairs")
    ax.text(0.53, 0.50, center, ha="center", va="center", fontsize=11, fontweight="bold")
    ax.text(0.27, 0.84, f"MOF loss\n(n={summary['total_MOF_loss_peaks']})", ha="center")
    ax.text(0.79, 0.84, f"H4K16ac loss\n(n={summary['total_H4K16ac_loss_peaks']})", ha="center")
    ax.set_xlim(0, 1.06)
    ax.set_ylim(0.12, 0.92)
    ax.axis("off")
    ax.set_title("MOF and H4K16ac loss-peak overlap", pad=15, fontweight="bold")
    fig.tight_layout()
    fig.savefig(plots / "loss_peak_overlap_venn.pdf")
    fig.savefig(plots / "loss_peak_overlap_venn.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.histplot(null, bins=min(40, max(10, len(np.unique(null)))), color="#4C78A8", ax=ax)
    ax.axvline(summary["permutation_observed"], color="#D55E00", linewidth=2.5,
               label=f"Observed = {summary['permutation_observed']}")
    ax.set(title="Chromosome-aware permutation test",
           xlabel="Random MOF-loss peaks overlapping H4K16ac-loss peaks",
           ylabel="Permutations")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(plots / "loss_overlap_permutation_test.pdf")
    fig.savefig(plots / "loss_overlap_permutation_test.png", dpi=200)
    plt.close(fig)


def summarize(args):
    output = Path(args.output)
    tables = output / "tables"
    beds = output / "beds"
    mof = read_bed(beds / "MOF_loss.bed")
    h4 = read_bed(beds / "H4K16ac_loss.bed")
    pairs = pd.read_csv(tables / "shared_loss_pairs_with_ids.tsv", sep="\t")

    mof_details = pd.read_csv(tables / "MOF_loss_DiffBind.csv").set_index("peak_id")
    h4_details = pd.read_csv(tables / "H4K16ac_loss_DiffBind.csv").set_index("peak_id")
    annotation_path = beds / "shared_loss_intervals_annotated.tsv"
    if annotation_path.exists() and annotation_path.stat().st_size > 0:
        annotation = pd.read_csv(annotation_path, sep="\t")
        annotation = annotation.drop_duplicates("V4").set_index("V4")
    else:
        annotation = pd.DataFrame(
            columns=["SYMBOL", "geneId", "annotation", "distanceToTSS"]
        )

    nearest_gene = annotation["SYMBOL"].where(
        annotation["SYMBOL"].notna(), annotation["geneId"]
    )

    shared = pd.DataFrame({
        "chromosome": pairs["MOF_chromosome"],
        "start": np.maximum(pairs["MOF_start"], pairs["H4K16ac_start"]).astype(int),
        "end": np.minimum(pairs["MOF_end"], pairs["H4K16ac_end"]).astype(int),
        "MOF_peak_id": pairs["MOF_peak_id"],
        "H4K16ac_peak_id": pairs["H4K16ac_peak_id"],
        "MOF_log2FC": pairs["MOF_peak_id"].map(mof_details["Fold"]),
        "H4K16ac_log2FC": pairs["H4K16ac_peak_id"].map(h4_details["Fold"]),
        "MOF_FDR": pairs["MOF_peak_id"].map(mof_details["FDR"]),
        "H4K16ac_FDR": pairs["H4K16ac_peak_id"].map(h4_details["FDR"]),
        "nearest_gene": pairs["shared_peak_id"].map(nearest_gene),
        "genomic_annotation": pairs["shared_peak_id"].map(annotation["annotation"]),
        "distance_to_TSS": pairs["shared_peak_id"].map(annotation["distanceToTSS"]),
    })
    shared.to_csv(tables / "MOF_H4K16ac_shared_loss_peaks.csv", index=False)

    overlapping_mof = pairs.MOF_peak_id.nunique()
    overlapping_h4 = pairs.H4K16ac_peak_id.nunique()
    summary = {
        "total_MOF_loss_peaks": len(mof),
        "total_H4K16ac_loss_peaks": len(h4),
        "overlapping_MOF_peaks": overlapping_mof,
        "overlapping_H4K16ac_peaks": overlapping_h4,
        "overlap_pairs": len(pairs),
        "MOF_only_peaks": len(mof) - overlapping_mof,
        "H4K16ac_only_peaks": len(h4) - overlapping_h4,
        "MOF_overlap_percentage": 100 * overlapping_mof / len(mof),
        "H4K16ac_overlap_percentage": 100 * overlapping_h4 / len(h4),
    }

    permutation, null = permutation_test(
        mof, h4, read_chrom_sizes(args.chrom_sizes), args.permutations, args.seed
    )
    pd.DataFrame([permutation]).to_csv(tables / "permutation_test_summary.csv", index=False)
    pd.DataFrame({"random_overlapping_MOF_peaks": null}).to_csv(
        tables / "permutation_null_distribution.csv", index=False
    )
    summary["permutation_observed"] = permutation["observed_overlapping_MOF_peaks"]
    pd.DataFrame([summary]).to_csv(tables / "loss_overlap_summary.csv", index=False)
    make_plots(summary, null, output)
    print(f"Results written to {output}", flush=True)


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    prep = subparsers.add_parser("prepare")
    prep.add_argument("--mof", required=True)
    prep.add_argument("--h4k16ac", required=True)
    prep.add_argument("--output", required=True)

    shared = subparsers.add_parser("build-shared")
    shared.add_argument("--output", required=True)

    finish = subparsers.add_parser("summarize")
    finish.add_argument("--output", required=True)
    finish.add_argument("--chrom-sizes", required=True)
    finish.add_argument("--permutations", type=int, default=1000)
    finish.add_argument("--seed", type=int, default=12345)
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    {"prepare": prepare, "build-shared": build_shared, "summarize": summarize}[
        arguments.command
    ](arguments)
