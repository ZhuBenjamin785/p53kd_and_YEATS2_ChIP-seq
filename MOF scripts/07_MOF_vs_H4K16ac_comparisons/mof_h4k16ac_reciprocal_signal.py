#!/usr/bin/env python3
"""Reciprocal signal and proximity analysis for MOF and H4K16ac loss peaks."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
import seaborn as sns
from scipy.stats import binomtest, wilcoxon


THRESHOLDS = (1000, 2000, 5000, 10000)


def chrom_name(value):
    value = str(value)
    return value if value.startswith("chr") else "chr" + value


def losses(path, label):
    data = pd.read_csv(path)
    required = {"seqnames", "start", "end", "Fold", "FDR"}
    if not required.issubset(data):
        raise ValueError(f"{path} lacks {sorted(required.difference(data.columns))}")
    data = data.dropna(subset=list(required)).copy()
    data = data[(data.FDR < 0.05) & (data.Fold < 0)].copy()
    data["chromosome"] = data.seqnames.map(chrom_name)
    data["start0"] = data.start.astype(int) - 1
    data["end"] = data.end.astype(int)
    data["peak_id"] = [f"{label}_loss_{i:06d}" for i in range(1, len(data) + 1)]
    return data.reset_index(drop=True)


def open_bigwigs(paths):
    return [pyBigWig.open(path) for path in paths]


def resolve_chrom(bw, chrom):
    names = bw.chroms()
    if chrom in names:
        return chrom
    alternate = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
    return alternate if alternate in names else None


def bw_mean(bw, chrom, start, end):
    chrom = resolve_chrom(bw, chrom)
    if chrom is None:
        return np.nan
    end = min(end, bw.chroms(chrom))
    if end <= start:
        return np.nan
    value = bw.stats(chrom, max(0, start), end, type="mean", exact=False)[0]
    return np.nan if value is None else float(value)


def reciprocal_signal(regions, control_paths, kd_paths, assay):
    controls, kd = open_bigwigs(control_paths), open_bigwigs(kd_paths)
    rows = []
    try:
        for peak in regions.itertuples(index=False):
            control_values = [bw_mean(x, peak.chromosome, peak.start0, peak.end) for x in controls]
            kd_values = [bw_mean(x, peak.chromosome, peak.start0, peak.end) for x in kd]
            control_mean, kd_mean = np.nanmean(control_values), np.nanmean(kd_values)
            row = {
                "chromosome": peak.chromosome, "start": peak.start0,
                "end": peak.end, "peak_id": peak.peak_id,
            }
            row.update({f"control_rep{i + 1}": value for i, value in enumerate(control_values)})
            row.update({f"KD_rep{i + 1}": value for i, value in enumerate(kd_values)})
            row.update({"control_mean": control_mean, "KD_mean": kd_mean,
                        f"{assay}_KD_minus_control": kd_mean - control_mean})
            rows.append(row)
    finally:
        for bw in controls + kd:
            bw.close()
    return pd.DataFrame(rows)


def signal_summary(table, delta_column, analysis):
    values = table[delta_column].dropna().to_numpy()
    negative = int((values < 0).sum())
    try:
        p_wilcoxon = float(wilcoxon(values, alternative="less").pvalue)
    except ValueError:
        p_wilcoxon = np.nan
    return {
        "analysis": analysis, "regions_tested": len(values),
        "median_KD_minus_control": float(np.median(values)),
        "mean_KD_minus_control": float(np.mean(values)),
        "regions_with_decrease": negative,
        "percent_with_decrease": 100 * negative / len(values),
        "sign_test_p_value": float(binomtest(negative, len(values), 0.5, alternative="greater").pvalue),
        "wilcoxon_less_p_value": p_wilcoxon,
    }


def nearest(query, target):
    groups = {chrom: x.reset_index(drop=True) for chrom, x in target.groupby("chromosome")}
    rows = []
    for peak in query.itertuples(index=False):
        candidates = groups.get(peak.chromosome)
        if candidates is None:
            rows.append((peak.peak_id, peak.chromosome, peak.start0, peak.end, None, np.nan))
            continue
        distances = np.maximum.reduce([
            candidates.start0.to_numpy() - peak.end,
            np.full(len(candidates), peak.start0) - candidates.end.to_numpy(),
            np.zeros(len(candidates), dtype=int),
        ])
        index = int(np.argmin(distances))
        rows.append((peak.peak_id, peak.chromosome, peak.start0, peak.end,
                     candidates.iloc[index].peak_id, int(distances[index])))
    return pd.DataFrame(rows, columns=["MOF_peak_id", "chromosome", "start", "end",
                                      "nearest_H4K16ac_peak_id", "distance_bp"])


def chrom_sizes(path):
    table = pd.read_csv(path, sep="\t", names=["chromosome", "length"])
    table.chromosome = table.chromosome.map(chrom_name)
    return dict(zip(table.chromosome, table.length.astype(int)))


def proximity_permutation(mof, h4, sizes, observed_distances, iterations, seed):
    mof = mof[mof.chromosome.isin(sizes)].reset_index(drop=True)
    h4 = h4[h4.chromosome.isin(sizes)].reset_index(drop=True)
    observed = np.array([(observed_distances.distance_bp <= x).sum() for x in THRESHOLDS])
    null = np.zeros((iterations, len(THRESHOLDS)), dtype=int)
    widths = (mof.end - mof.start0).to_numpy()
    rng = np.random.default_rng(seed)
    for iteration in range(iterations):
        randomized = mof.copy()
        randomized.start0 = [rng.integers(0, max(1, sizes[c] - w + 1))
                             for c, w in zip(randomized.chromosome, widths)]
        randomized.end = randomized.start0 + widths
        distances = nearest(randomized, h4).distance_bp
        null[iteration] = [(distances <= x).sum() for x in THRESHOLDS]
    rows = []
    for index, threshold in enumerate(THRESHOLDS):
        expected = null[:, index].mean()
        rows.append({
            "threshold_bp": threshold, "observed_MOF_peaks": observed[index],
            "observed_percentage": 100 * observed[index] / len(mof),
            "mean_random_MOF_peaks": expected,
            "enrichment_fold": observed[index] / expected if expected else np.nan,
            "empirical_p_value": (1 + (null[:, index] >= observed[index]).sum()) / (iterations + 1),
            "permutations": iterations, "seed": seed,
        })
    return pd.DataFrame(rows), null


def profile(regions, paths, flank, bins):
    bigwigs = open_bigwigs(paths)
    rows = []
    try:
        for peak in regions.itertuples(index=False):
            center = (peak.start0 + peak.end) // 2
            replicate_rows = []
            for bw in bigwigs:
                chrom = resolve_chrom(bw, peak.chromosome)
                if chrom is None or center - flank < 0 or center + flank > bw.chroms(chrom):
                    replicate_rows.append(np.full(bins, np.nan))
                    continue
                values = bw.stats(chrom, center - flank, center + flank,
                                  nBins=bins, exact=False)
                replicate_rows.append(np.array([np.nan if x is None else x for x in values]))
            replicate_rows = np.asarray(replicate_rows, dtype=float)
            rows.append(
                np.full(bins, np.nan)
                if np.isnan(replicate_rows).all()
                else np.nanmean(replicate_rows, axis=0)
            )
    finally:
        for bw in bigwigs:
            bw.close()
    return np.vstack(rows)


def metaplot(regions, control_paths, kd_paths, assay, centered_on, output, flank, bins):
    control = profile(regions, control_paths, flank, bins)
    kd = profile(regions, kd_paths, flank, bins)
    x = np.linspace(-flank, flank, bins)
    profiles = pd.DataFrame({
        "distance_bp": x, "control_mean": np.nanmean(control, axis=0),
        "KD_mean": np.nanmean(kd, axis=0), "KD_minus_control": np.nanmean(kd - control, axis=0),
    })
    profiles.to_csv(output / "tables" / f"{assay}_at_{centered_on}_metaprofile.csv", index=False)

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True,
                             gridspec_kw={"height_ratios": [2, 1]})
    axes[0].plot(x, profiles.control_mean, label="Scramble", color="#4C78A8", lw=2.4)
    axes[0].plot(x, profiles.KD_mean, label="p53 depleted", color="#E45756", lw=2.4)
    axes[0].set(ylabel="Mean log2(ChIP/Input)",
                title=f"{assay} signal centered on {centered_on} loss sites")
    axes[0].legend(frameon=False)
    axes[1].plot(x, profiles.KD_minus_control, color="#B279A2", lw=2.4)
    axes[1].axhline(0, color="0.4", ls="--", lw=1)
    axes[1].set(xlabel="Distance from peak center (bp)", ylabel="KD − Scramble")
    for ax in axes:
        ax.axvline(0, color="0.5", ls=":", lw=1)
        sns.despine(ax=ax)
    fig.tight_layout()
    stem = output / "plots" / f"{assay}_at_{centered_on}_metaplot"
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=200)
    plt.close(fig)


def main(args):
    output = Path(args.output)
    for child in ["tables", "plots"]:
        (output / child).mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="white", context="talk")

    mof_loss, h4_loss = losses(args.mof_diffbind, "MOF"), losses(args.h4k16ac_diffbind, "H4K16ac")
    mof_at_h4 = reciprocal_signal(h4_loss, args.mof_control, args.mof_kd, "MOF")
    h4_at_mof = reciprocal_signal(mof_loss, args.h4_control, args.h4_kd, "H4K16ac")
    mof_at_h4.to_csv(output / "tables" / "MOF_change_at_H4K16ac_loss_sites.csv", index=False)
    h4_at_mof.to_csv(output / "tables" / "H4K16ac_change_at_MOF_loss_sites.csv", index=False)
    pd.DataFrame([
        signal_summary(mof_at_h4, "MOF_KD_minus_control", "MOF at H4K16ac-loss sites"),
        signal_summary(h4_at_mof, "H4K16ac_KD_minus_control", "H4K16ac at MOF-loss sites"),
    ]).to_csv(output / "tables" / "reciprocal_signal_summary.csv", index=False)

    distances = nearest(mof_loss, h4_loss)
    distances.to_csv(output / "tables" / "MOF_to_nearest_H4K16ac_loss_peak.csv", index=False)
    proximity, null = proximity_permutation(
        mof_loss, h4_loss, chrom_sizes(args.chrom_sizes), distances,
        args.permutations, args.seed
    )
    proximity.to_csv(output / "tables" / "distance_threshold_permutation_summary.csv", index=False)
    pd.DataFrame(null, columns=[f"within_{x}bp" for x in THRESHOLDS]).to_csv(
        output / "tables" / "distance_permutation_null.csv", index=False
    )

    fig, ax = plt.subplots(figsize=(8, 6))
    plot_data = proximity.assign(threshold=lambda x: x.threshold_bp / 1000)
    sns.barplot(data=plot_data, x="threshold", y="observed_percentage", color="#B279A2", ax=ax)
    ax.set(title="MOF-loss peaks near H4K16ac-loss peaks", xlabel="Maximum distance (kb)",
           ylabel="MOF-loss peaks within threshold (%)")
    fig.tight_layout()
    fig.savefig(output / "plots" / "distance_threshold_fractions.pdf")
    fig.savefig(output / "plots" / "distance_threshold_fractions.png", dpi=200)
    plt.close(fig)

    metaplot(h4_loss, args.mof_control, args.mof_kd, "MOF", "H4K16ac",
             output, args.flank, args.bins)
    metaplot(mof_loss, args.h4_control, args.h4_kd, "H4K16ac", "MOF",
             output, args.flank, args.bins)
    print(f"Results written to {output}", flush=True)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mof-diffbind", required=True)
    parser.add_argument("--h4k16ac-diffbind", required=True)
    parser.add_argument("--mof-control", nargs=2, required=True)
    parser.add_argument("--mof-kd", nargs=2, required=True)
    parser.add_argument("--h4-control", nargs=2, required=True)
    parser.add_argument("--h4-kd", nargs=2, required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--flank", type=int, default=5000)
    parser.add_argument("--bins", type=int, default=100)
    return parser.parse_args()


if __name__ == "__main__":
    main(parse_args())
