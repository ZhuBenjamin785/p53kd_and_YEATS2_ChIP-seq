#!/usr/bin/env python3
                      
"""Rank protein-coding genes by H4K16ac change from Scramble to p53 KD."""

import gzip
import os
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd

                                                                             
mplconfigdir = Path("/tmp/pybw-mplconfig")
mplconfigdir.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(mplconfigdir))

import seaborn as sns
import pyBigWig

gtf = "gencode.v50.basic.annotation.gtf"
tracks = {
    "Scramble_rep1": "bamcompare/Scramble_rep1_ChIP_vs_Input_log2.bw",
    "Scramble_rep2": "bamcompare/Scramble_rep2_ChIP_vs_Input_log2.bw",
    "KD_rep1": "bamcompare/P53KD_rep1_ChIP_vs_Input_log2.bw",
    "KD_rep2": "bamcompare/P53KD_rep2_ChIP_vs_Input_log2.bw",
}
outdir = Path("shift_results")
promoterbp = 2000
top_n = 50


def load_protein_coding_genes(path):
    """Read only protein-coding gene records; gtf is 1-based inclusive."""
    opener = gzip.open if str(path).endswith(".gz") else open
    genes = []
    wanted = re.compile(r'(gene_id|gene_name|gene_type) "([^"]+)"')
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = dict(wanted.findall(fields[8]))
            if attrs.get("gene_type") != "protein_coding":
                continue
            genes.append((attrs["gene_id"], attrs.get("gene_name", attrs["gene_id"]),
                          fields[0], int(fields[3]) - 1, int(fields[4]), fields[6]))
    return pd.DataFrame(genes, columns=["gene_id", "Gene", "chr", "start", "end", "strand"])


def validate_bigwig_files(paths):
    bad = []
    for path in paths:
        try:
            bw = pyBigWig.open(path)
        except RuntimeError:
            bad.append(path)
            continue
        if bw is None:
            bad.append(path)
            continue
        bw.close()
    if bad:
        raise RuntimeError(
            "These bigWig files cannot be opened and must be regenerated: "
            + ", ".join(bad)
        )


def get_regions(genes):
    tss = np.where(genes["strand"].eq("+"), genes["start"], genes["end"] - 1)
    promoter = list(zip(genes["chr"], np.maximum(0, tss - promoterbp),
                        tss + promoterbp))
    body = list(zip(genes["chr"], genes["start"], genes["end"]))
    return promoter, body


def summarize_track(job):
    """Return exact mean signal for promoter and body regions in one bigWig."""
    name, path, promoter, body = job
    try:
        bigwig = pyBigWig.open(path)
    except RuntimeError as exc:
        raise RuntimeError(f"Could not open bigWig file: {path}") from exc
    chroms = bigwig.chroms()

    def means(regions):
        result = []
        for chrom, start, end in regions:
            if chrom not in chroms:
                result.append(np.nan)
                continue
            end = min(int(end), chroms[chrom])
            value = bigwig.stats(chrom, int(start), end, type="mean", exact=True)[0]
            result.append(0.0 if value is None else value)
        return result

    promoter_signal, body_signal = means(promoter), means(body)
    bigwig.close()
    return name, promoter_signal, body_signal


def filter_low_signal_genes(genes, promoter_signals, body_signals):
    """Keep genes with finite signal in every ChIP/Input log2 track."""
    signals = [*promoter_signals.values(), *body_signals.values()]
    keep = np.logical_and.reduce([np.isfinite(signal) for signal in signals])
    print(f"Finite-signal filter retained {keep.sum():,}/{len(genes):,} genes", flush=True)
    return keep


def rank_changes(genes, signals, region, keep_mask=None):
    table = genes.copy()
    for sample in tracks:
        table[sample] = signals[sample]
    table["Scramble_signal"] = table[["Scramble_rep1", "Scramble_rep2"]].mean(axis=1)
    table["KD_signal"] = table[["KD_rep1", "KD_rep2"]].mean(axis=1)
    if keep_mask is not None:
        table = table.loc[keep_mask].copy()
    table["Difference"] = table["KD_signal"] - table["Scramble_signal"]
    table["log2FC"] = table["Difference"]
    table = table.replace([np.inf, -np.inf], np.nan).dropna()
    table = table[(table["Scramble_signal"] + table["KD_signal"]) > 0].copy()
    table["Direction"] = np.where(table["Difference"] >= 0, "gain", "loss")
    table["Rank_gain"] = table["Difference"].rank(method="min", ascending=False).astype(int)
    table["Rank_loss"] = table["Difference"].rank(method="min", ascending=True).astype(int)
    table["Rank"] = table["Difference"].abs().rank(method="min", ascending=False).astype(int)
    table = table.sort_values("Rank")
    table.to_csv(outdir / f"ranked_{region}_H4K16ac_changes.csv", index=False)
    return table


def tss_profile(rows, paths, bins=100, flank=3000, exact=True):
    profiles = []
    handles = [pyBigWig.open(path) for path in paths]
    for row in rows.itertuples():
        tss = row.start if row.strand == "+" else row.end - 1
        per_sample = []
        for bigwig in handles:
            chrom_size = bigwig.chroms().get(row.chr)
            if chrom_size is None:
                continue
            start, end = max(0, tss - flank), min(chrom_size, tss + flank)
            values = bigwig.stats(row.chr, start, end, nBins=bins, exact=exact)
            per_sample.append(np.array([0.0 if value is None else value
                                        for value in values]))
        if per_sample:
            profiles.append(np.mean(per_sample, axis=0))
    for bigwig in handles:
        bigwig.close()
    return np.mean(profiles, axis=0) if profiles else np.zeros(bins)


def plot_metaplot(rows):
    """Plot mean ChIP/Input signal around TSSs for all retained genes."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    scr_paths = [tracks["Scramble_rep1"], tracks["Scramble_rep2"]]
    kd_paths = [tracks["KD_rep1"], tracks["KD_rep2"]]
    bins = 50
    x = np.linspace(-3, 3, bins)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(x, tss_profile(rows, scr_paths, bins=bins, exact=False),
            label="Scramble", color="blue")
    ax.plot(x, tss_profile(rows, kd_paths, bins=bins, exact=False),
            label="p53 KD", color="red")
    ax.axvline(0, color="grey", ls="--", lw=1)
    ax.set(xlabel="Distance from TSS (kb)", ylabel="Mean log2(ChIP/Input)",
           title=f"H4K16ac metaplot across {len(rows):,} protein-coding genes")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "H4K16ac_metaplot_all_genes.pdf")
    plt.close(fig)


"""
use bam compare for normalizing signal to input for each sample

bpm
bigwig compare

venn diagram of peaks for kd v scramble
only kd
macs3 peak intensity compare
"""





def make_plots(table, region):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.backends.backend_pdf import PdfPages

    gained = table.nlargest(top_n, "Difference")
    lost = table.nsmallest(top_n, "Difference")
    scr_paths = [tracks["Scramble_rep1"], tracks["Scramble_rep2"]]
    kd_paths = [tracks["KD_rep1"], tracks["KD_rep2"]]
    samples = list(tracks)
    top = table.head(top_n)
    heat = top[samples].to_numpy(float)
    row_sd = heat.std(axis=1, keepdims=True)
    heat = (heat - heat.mean(axis=1, keepdims=True)) / np.where(row_sd == 0, 1, row_sd)
    sns.set_theme(style="whitegrid")

    with PdfPages(outdir / f"{region}_shift_plots.pdf") as pdf:
        fig, ax = plt.subplots(figsize=(7, 10))
        sns.heatmap(heat, cmap="vlag", center=0, xticklabels=["Scr1", "Scr2", "KD1", "KD2"],
                    yticklabels=top["Gene"], ax=ax, cbar_kws={"label": "Row z-score"})
        ax.set_title(f"Top {len(top)} absolute {region} shifts")
        ax.tick_params(axis="y", labelsize=5); fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        fig, ax = plt.subplots(figsize=(6, 6))
        sns.scatterplot(data=table, x="Scramble_signal", y="KD_signal", s=10,
                        alpha=.35, linewidth=0, ax=ax)
        high = max(table["Scramble_signal"].max(), table["KD_signal"].max())
        ax.plot([0, high], [0, high], "r--", lw=1); ax.set_title(f"{region}: KD vs Scramble")
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        fig, ax = plt.subplots(figsize=(7, 5))
        average = (table["Scramble_signal"] + table["KD_signal"]) / 2
        ax.scatter(average, table["log2FC"], s=8, alpha=.35)
        ax.axhline(0, color="red", ls="--", lw=1)
        ax.set(xlabel="mean log2(ChIP/Input)", ylabel="log2FC (KD/Scramble)",
               title=f"{region} change plot")
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        x = np.linspace(-3, 3, 100)
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
        for ax, rows, label in [(axes[0], gained, "gained"), (axes[1], lost, "lost")]:
            ax.plot(x, tss_profile(rows, scr_paths), label="Scramble", color="black")
            ax.plot(x, tss_profile(rows, kd_paths), label="p53 KD", color="red")
            ax.axvline(0, color="grey", ls="--", lw=1)
            ax.set(title=f"Top {top_n} {label}", xlabel="Distance from TSS (kb)")
            ax.legend()
        axes[0].set_ylabel("Mean H4K16ac signal")
        fig.suptitle(f"TSS profiles selected by {region} shift")
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        
def ratiofinder(genes, promoter_signals, body_signals, keep_mask=None):
    table = genes.copy()

    for sample in tracks:
        table[f"{sample}_promoter"] = promoter_signals[sample]
        table[f"{sample}_body"] = body_signals[sample]

                        
    table["Scramble_promoter"] = table[
        ["Scramble_rep1_promoter", "Scramble_rep2_promoter"]
    ].mean(axis=1)

    table["KD_promoter"] = table[
        ["KD_rep1_promoter", "KD_rep2_promoter"]
    ].mean(axis=1)

    table["Scramble_body"] = table[
        ["Scramble_rep1_body", "Scramble_rep2_body"]
    ].mean(axis=1)

    table["KD_body"] = table[
        ["KD_rep1_body", "KD_rep2_body"]
    ].mean(axis=1)

    if keep_mask is not None:
        table = table.loc[keep_mask].copy()

    rep1 = (table["KD_rep1_body"] - table["KD_rep1_promoter"]
            - table["Scramble_rep1_body"] + table["Scramble_rep1_promoter"])
    rep2 = (table["KD_rep2_body"] - table["KD_rep2_promoter"]
            - table["Scramble_rep2_body"] + table["Scramble_rep2_promoter"])
    table = table.loc[rep1 * rep2 > 0].copy()
    table["ratio_log2FC"] = ((rep1 + rep2) / 2).loc[table.index]

    table["Rank"] = (
        table["ratio_log2FC"]
        .abs()
        .rank(method="min", ascending=False)
        .astype(int)
    )

    return table.sort_values("Rank")

def plot_ratio_logfc(table):
    import seaborn as sns
    import matplotlib.pyplot as plt
    if table.empty:
        raise ValueError("No genes remain for ratio plotting after low-signal filtering")
    gained = table.nlargest(25, "ratio_log2FC")
    lost = table.nsmallest(25, "ratio_log2FC")

    plot_df = pd.concat([lost, gained]).sort_values("ratio_log2FC")

    fig, ax = plt.subplots(figsize=(8, 10))

    ax.barh(
        plot_df["Gene"],
        plot_df["ratio_log2FC"]
    )

    ax.axvline(0, linestyle="--")

    ax.set_xlabel("log2 fold change of gene-body/promoter ratio")
    ax.set_ylabel("Gene")
    ax.set_title("H4K16ac redistribution after p53 KD")

    plt.tight_layout()
    plt.savefig(
        outdir / "body_promoter_ratio_log2FC.pdf",
        bbox_inches="tight"
    )
    plt.close()

def main():
    for path in [gtf, *tracks.values()]:
        if not Path(path).is_file():
            raise FileNotFoundError(path)
    validate_bigwig_files([*tracks.values()])

    outdir.mkdir(exist_ok=True)
    genes = load_protein_coding_genes(gtf)
    print(f"Loaded {len(genes):,} protein-coding genes", flush=True)
    promoter, body = get_regions(genes)
    jobs = [(name, path, promoter, body) for name, path in tracks.items()]
    promoter_signals, body_signals = {}, {}
    
    try:
        with ProcessPoolExecutor(max_workers=min(4, int(os.getenv("SLURM_CPUS_PER_TASK", "4")))) as pool:
            for name, prom, gene_body in pool.map(summarize_track, jobs):
                promoter_signals[name], body_signals[name] = prom, gene_body
                print(f"Finished {name}", flush=True)
    except (PermissionError, OSError):
        for job in jobs:
            name, prom, gene_body = summarize_track(job)
            promoter_signals[name], body_signals[name] = prom, gene_body
            print(f"Finished {name}", flush=True)
    
    keep_mask = filter_low_signal_genes(genes, promoter_signals, body_signals)
    plot_metaplot(genes.loc[keep_mask])
    for region, signals in [("promoter", promoter_signals), ("genebody", body_signals)]:
        ranked = rank_changes(genes, signals, region, keep_mask)
        make_plots(ranked, region)
        print(f"Wrote {region} results for {len(ranked):,} genes", flush=True)
        
    ratios = ratiofinder(
        genes,
        promoter_signals,
        body_signals,
        keep_mask
    )

    ratios.to_csv(
        outdir / "ranked_body_promoter_ratios.csv",
        index=False
    )
    
    plot_ratio_logfc(ratios)



if __name__ == "__main__":
    main()
