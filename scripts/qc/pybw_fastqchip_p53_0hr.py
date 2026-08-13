#!/usr/bin/env python3
                      
"""Summarize p53 0-hour ChIP/input signal around hg38 gene regions."""

import gzip
import os
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

mplconfigdir = Path("/tmp/pybw-fastqchip-p53-0hr-mplconfig")
mplconfigdir.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(mplconfigdir))

import pyBigWig

gtf = "gencode.v50.basic.annotation.gtf"
tracks = {
    "p53_0hr_ChIP_vs_input":
    "fastqchip_macs3_results/bamcompare/p53_0hr_ChIP_vs_input_log2.bw"
}
outdir = Path("fastqchip_p53_0hr_signal_results")
promoterbp = 2000
top_n = 50


def load_protein_coding_genes(path):
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
            genes.append((
                attrs["gene_id"], attrs.get("gene_name", attrs["gene_id"]),
                fields[0], int(fields[3]) - 1, int(fields[4]), fields[6]
            ))
    return pd.DataFrame(
        genes, columns=["gene_id", "Gene", "chr", "start", "end", "strand"]
    )


def validate_bigwig(path):
    bw = pyBigWig.open(path)
    if bw is None:
        raise RuntimeError(f"Could not open bigWig: {path}")
    bw.close()


def get_regions(genes):
    tss = np.where(genes["strand"].eq("+"), genes["start"], genes["end"] - 1)
    promoter = list(zip(genes["chr"], np.maximum(0, tss - promoterbp), tss + promoterbp))
    body = list(zip(genes["chr"], genes["start"], genes["end"]))
    return promoter, body


def summarize_track(job):
    name, path, promoter, body = job
    bigwig = pyBigWig.open(path)
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


def tss_profile(rows, path, bins=100, flank=3000):
    bigwig = pyBigWig.open(path)
    profiles = []
    for row in rows.itertuples():
        tss = row.start if row.strand == "+" else row.end - 1
        chrom_size = bigwig.chroms().get(row.chr)
        if chrom_size is None:
            continue
        start, end = max(0, tss - flank), min(chrom_size, tss + flank)
        values = bigwig.stats(row.chr, start, end, nBins=bins, exact=False)
        profiles.append(np.array([0.0 if v is None else v for v in values]))
    bigwig.close()
    return np.mean(profiles, axis=0) if profiles else np.zeros(bins)


def make_plots(table, region):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.backends.backend_pdf import PdfPages

    top = table.reindex(table["signal"].abs().sort_values(ascending=False).index).head(top_n)
    matrix = top[["signal"]].to_numpy(float)
    sns.set_theme(style="whitegrid")
    with PdfPages(outdir / f"{region}_signal_plots.pdf") as pdf:
        fig, ax = plt.subplots(figsize=(6, 10))
        sns.heatmap(matrix, cmap="vlag", center=0, xticklabels=["p53 0hr"],
                    yticklabels=top["Gene"], ax=ax,
                    cbar_kws={"label": "Mean log2(ChIP/Input)"})
        ax.set_title(f"Top {len(top)} p53 0-hour {region} signals")
        ax.tick_params(axis="y", labelsize=5)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        fig, ax = plt.subplots(figsize=(7, 5))
        ax.hist(table["signal"].dropna(), bins=100, color="steelblue")
        ax.axvline(0, color="red", ls="--", lw=1)
        ax.set(xlabel="Mean log2(ChIP/Input)", ylabel="Genes",
               title=f"p53 0-hour {region} signal distribution")
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def main():
    for path in [gtf, *tracks.values()]:
        if not Path(path).is_file():
            raise FileNotFoundError(path)
        if path.endswith(".bw"):
            validate_bigwig(path)
    outdir.mkdir(exist_ok=True)
    genes = load_protein_coding_genes(gtf)
    promoter, body = get_regions(genes)
    jobs = [(name, path, promoter, body) for name, path in tracks.items()]
    result = summarize_track(jobs[0])
    _, promoter_signal, body_signal = result
    keep = np.isfinite(promoter_signal) & np.isfinite(body_signal)
    print(f"Finite-signal filter retained {keep.sum():,}/{len(genes):,} genes", flush=True)

    tables = {}
    for region, signal in [("promoter", promoter_signal), ("genebody", body_signal)]:
        table = genes.loc[keep].copy()
        table["signal"] = np.asarray(signal)[keep]
        table["rank"] = table["signal"].abs().rank(method="min", ascending=False).astype(int)
        table = table.sort_values("rank")
        table.to_csv(outdir / f"ranked_{region}_p53_0hr_signal.csv", index=False)
        make_plots(table, region)
        tables[region] = table

    profile = tss_profile(genes.loc[keep], tracks["p53_0hr_ChIP_vs_input"])
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = np.linspace(-3, 3, len(profile))
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(x, profile, color="firebrick")
    ax.axvline(0, color="grey", ls="--", lw=1)
    ax.set(xlabel="Distance from TSS (kb)", ylabel="Mean log2(ChIP/Input)",
           title=f"p53 0-hour signal around TSSs ({keep.sum():,} genes)")
    fig.tight_layout(); fig.savefig(outdir / "p53_0hr_TSS_metaplot.pdf"); plt.close(fig)


if __name__ == "__main__":
    main()
