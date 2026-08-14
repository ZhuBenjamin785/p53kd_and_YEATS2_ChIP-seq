#!/usr/bin/env python3
                      
"""Consistent matplotlib/seaborn figures for the three ChIP-seq datasets."""

import gzip
import re
import sys
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
import seaborn as sns

root = Path.cwd()
out_root = root / "chipseq_summary_plots"
gtf = root / "gencode.v50.basic.annotation.gtf"
control_color = "#4C78A8"
p53_color = "#E45756"
yeats2_color = "#72B7B2"
mof_color = "#B279A2"
gain_color = "#D55E00"
loss_color = "#0072B2"

datasets = {
    "MOF": {
        "title": "MOF occupancy after p53 knockdown",
        "tracks": {
            "Scramble MOF 1": "mof_macs3_results/bamcompare/Scr_MOF_rep1_ChIP_vs_Input_log2.bw",
            "Scramble MOF 2": "mof_macs3_results/bamcompare/Scr_MOF_rep2_ChIP_vs_Input_log2.bw",
            "p53sh MOF 1": "mof_macs3_results/bamcompare/p53sh_MOF_rep1_ChIP_vs_Input_log2.bw",
            "p53sh MOF 2": "mof_macs3_results/bamcompare/p53sh_MOF_rep2_ChIP_vs_Input_log2.bw",
        },
        "groups": {
            "Scramble": ["Scramble MOF 1", "Scramble MOF 2"],
            "p53sh": ["p53sh MOF 1", "p53sh MOF 2"],
        },
        "colors": {"Scramble": control_color, "p53sh": mof_color},
        "diff": "mof_macs3_results/diffbind_results/DiffBind_all_peaks.csv",
    },
    "p53KD": {
        "title": "p53 KD H4K16ac",
        "tracks": {
            "Scramble 1": "bamcompare/Scramble_rep1_ChIP_vs_Input_log2.bw",
            "Scramble 2": "bamcompare/Scramble_rep2_ChIP_vs_Input_log2.bw",
            "p53 KD 1": "bamcompare/P53KD_rep1_ChIP_vs_Input_log2.bw",
            "p53 KD 2": "bamcompare/P53KD_rep2_ChIP_vs_Input_log2.bw",
        },
        "groups": {"Scramble": ["Scramble 1", "Scramble 2"], "p53 KD": ["p53 KD 1", "p53 KD 2"]},
        "colors": {"Scramble": control_color, "p53 KD": p53_color},
        "diff": "diffbind_results/DiffBind_all_peaks.csv",
    },
    "YEATS2KD": {
        "title": "YEATS2 KD H4K16ac",
        "tracks": {
            "Scramble 1": "bamcompareY2/Scramble_YEATS2_rep1_ChIP_vs_Input_log2.bw",
            "Scramble 2": "bamcompareY2/Scramble_YEATS2_rep2_ChIP_vs_Input_log2.bw",
            "YEATS2 KD 1": "bamcompareY2/YEATS2KD_rep1_ChIP_vs_Input_log2.bw",
            "YEATS2 KD 2": "bamcompareY2/YEATS2KD_rep2_ChIP_vs_Input_log2.bw",
        },
        "groups": {"Scramble": ["Scramble 1", "Scramble 2"], "YEATS2 KD": ["YEATS2 KD 1", "YEATS2 KD 2"]},
        "colors": {"Scramble": control_color, "YEATS2 KD": yeats2_color},
        "diff": "diffbind_results_yeats2/DiffBind_YEATS2_all_peaks.csv",
    },
    "p53_0hr": {
        "title": "p53 ChIP, 0 hr TNF-α",
        "tracks": {"p53, 0 hr": "fastqchip_macs3_results/bamcompare/p53_0hr_ChIP_vs_input_log2.bw"},
        "groups": {"p53, 0 hr": ["p53, 0 hr"]},
        "colors": {"p53, 0 hr": p53_color},
        "diff": None,
    },
}

sns.set_theme(style="white", context="notebook", font_scale=1.15)
plt.rcParams.update({"axes.titleweight": "bold", "axes.titlesize": 14,
                     "axes.labelsize": 11, "legend.fontsize": 10,
                     "figure.dpi": 150, "savefig.bbox": "tight"})


def symmetric_bound(values, minimum=1e-9, pad=1.08):
    """Return a symmetric +/- limit for zero-centered signal plots."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    bound = np.nanmax(np.abs(values)) if values.size else minimum
    return max(minimum, bound * pad)


def save_figure(fig, pdf_path):
    """Save the publication PDF and a matching PNG used for contact sheets."""
    pdf_path = Path(pdf_path)
    fig.savefig(pdf_path)
    fig.savefig(pdf_path.with_suffix(".png"), dpi=150)


def make_master_pdf(dataset, outdir, columns=2):
    """Assemble all dataset PDFs into readable multi-panel master pages.

    Every figure created by this program has a same-named PNG preview. The
    previews allow contact-sheet assembly without requiring a PDF renderer.
    The original individual PDFs are retained alongside the master PDF.
    """
    pdfs = sorted(
        p for p in outdir.glob("*.pdf")
        if not p.name.startswith("MASTER_")
    )
    previews = [(pdf, pdf.with_suffix(".png")) for pdf in pdfs if pdf.with_suffix(".png").is_file()]
    if not previews:
        return
    rows = 2
    page_width, page_height = 11.7, 8.3
    master_path = outdir / f"MASTER_{dataset}_figures.pdf"
    with PdfPages(master_path) as master:
        for start in range(0, len(previews), columns * rows):
            page_items = previews[start:start + columns * rows]
            fig, axes = plt.subplots(rows, columns, figsize=(page_width, page_height))
            axes = np.atleast_1d(axes).ravel()
            fig.suptitle(f"{datasets[dataset]['title']} — figure summary", fontsize=16, fontweight="bold")
            for ax, (pdf, png) in zip(axes, page_items):
                image = plt.imread(png)
                ax.imshow(image)
                ax.set_title(pdf.stem.replace("_", " "), fontsize=10)
                ax.axis("off")
            for ax in axes[len(page_items):]:
                ax.axis("off")
            fig.tight_layout(rect=(0, 0, 1, .95))
            master.savefig(fig)
            plt.close(fig)
    print(f"Wrote master PDF: {master_path}", flush=True)


def load_genes():
    genes, pattern = [], re.compile(r'(gene_id|gene_name|gene_type) "([^"]+)"')
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"): continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "gene": continue
            attrs = dict(pattern.findall(fields[8]))
            if attrs.get("gene_type") != "protein_coding": continue
            start, end = int(fields[3]) - 1, int(fields[4])
            tss = start if fields[6] == "+" else end - 1
            genes.append((attrs.get("gene_id"), attrs.get("gene_name"), fields[0],
                          start, end, tss, fields[6]))
    return pd.DataFrame(genes, columns=["gene_id", "Gene", "chr", "start", "end", "tss", "strand"])


def tss_matrices(track_paths, genes, bins=100, flank=2000):
    selected = genes.iloc[:10000].copy()
    matrices = {}
    for name, relative_path in track_paths.items():
        bw = pyBigWig.open(str(root / relative_path)); chroms = bw.chroms(); rows = []
        for row in selected.itertuples():
            if row.chr not in chroms:
                rows.append(np.full(bins, np.nan)); continue
            start, end = max(0, row.tss - flank), min(chroms[row.chr], row.tss + flank)
            if end - start < bins:
                rows.append(np.full(bins, np.nan)); continue
            values = bw.stats(row.chr, start, end, nBins=bins, exact=False)
            rows.append(np.array([0.0 if value is None else value for value in values], dtype=float))
        bw.close(); matrices[name] = np.vstack(rows)
    return selected, matrices


def group_mean(matrices, members):
    return np.nanmean(np.stack([matrices[name] for name in members]), axis=0)


def gene_body_matrices(track_paths, genes, bins=100):
    """Return signal across each gene from TSS to TES in 5-prime to 3-prime order."""
    selected = genes.iloc[:10000].copy()
    matrices = {}
    for name, relative_path in track_paths.items():
        bw = pyBigWig.open(str(root / relative_path))
        chroms = bw.chroms()
        rows = []
        for row in selected.itertuples():
            if row.chr not in chroms or row.end <= row.start:
                rows.append(np.full(bins, np.nan))
                continue
            start, end = max(0, row.start), min(chroms[row.chr], row.end)
            values = bw.stats(row.chr, start, end, nBins=bins, exact=False)
            values = np.array([0.0 if value is None else value for value in values], dtype=float)
            if row.strand == "-":
                values = values[::-1]
            rows.append(values)
        bw.close()
        matrices[name] = np.vstack(rows)
    return selected, matrices


def plot_tss_to_tes(cfg, outdir, genes):
    selected, matrices = gene_body_matrices(cfg["tracks"], genes)
    group_matrices = {name: group_mean(matrices, members)
                      for name, members in cfg["groups"].items()}
    group_names = list(group_matrices)
    x = np.linspace(0, 100, next(iter(group_matrices.values())).shape[1])

                                                                            
                                                                       
    profile_table = pd.DataFrame({"normalized_gene_position_percent": x})
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for group, matrix in group_matrices.items():
        profile = np.nanmean(matrix, axis=0)
        profile_table[group] = profile
        valid = np.isfinite(matrix)
        n = valid.sum(axis=0)
        sem = np.nanstd(matrix, axis=0, ddof=1) / np.sqrt(np.maximum(n, 1))
        sem[n < 2] = 0
        ax.plot(x, profile, lw=2.5, label=group, color=cfg["colors"][group])
        ax.fill_between(x, profile - sem, profile + sem,
                        color=cfg["colors"][group], alpha=0.16, linewidth=0)
    ax.axvline(0, color="0.45", ls="--", lw=1)
    ax.axvline(100, color="0.45", ls="--", lw=1)
    ax.set(title=f"{cfg['title']}: TSS-to-TES metaplot",
           xlabel="Normalized gene position", ylabel="Mean log2(ChIP/Input) signal")
    ax.set_xticks([0, 25, 50, 75, 100], ["TSS", "25%", "50%", "75%", "TES"])
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    save_figure(fig, outdir / "TSS_to_TES_metaplot.pdf")
    plt.close(fig)
    profile_table.to_csv(outdir / "TSS_to_TES_metaplot.csv", index=False)

    if len(group_names) == 2:
        heat = group_matrices[group_names[1]] - group_matrices[group_names[0]]
        color_label = "Δ log2(ChIP/Input)\nKD − Scramble"
        row_label = "Protein-coding genes\nranked by |KD − Scramble|"
        cmap = "RdBu_r"
        vmax = np.nanpercentile(np.abs(heat), 99)
        vmin = -vmax
    else:
        heat = group_matrices[group_names[0]]
        color_label = "log2(p53 ChIP/Input)\ndm6-normalized BAMs"
        row_label = "Protein-coding genes\nranked by |signal|"
        vmin, vmax = np.nanpercentile(heat, [1, 99])
        if vmin < 0:
            vmax = max(abs(vmin), abs(vmax))
            vmin = -vmax
            cmap = "RdBu_r"
        else:
            cmap = "magma"

    score = np.nanmean(np.abs(heat), axis=1)
    keep = np.isfinite(score)
    order = np.argsort(score[keep])[::-1][:3000]
    heat = heat[keep][order]
    selected = selected.loc[keep].iloc[order]
    if not np.isfinite(vmax) or vmax <= vmin:
        vmax = symmetric_bound(heat)
        vmin = -vmax if len(group_names) == 2 else 0
    heat = np.ma.masked_invalid(np.clip(heat, vmin, vmax))

    fig, ax = plt.subplots(figsize=(8.5, 8))
    image = ax.imshow(heat, aspect="auto", interpolation="nearest", cmap=cmap,
                      vmin=vmin, vmax=vmax, extent=(0, 100, heat.shape[0], 0))
    ax.set(title=f"{cfg['title']}: TSS-to-TES metagene signal",
           xlabel="Normalized gene position", ylabel=row_label)
    ax.set_xticks([0, 25, 50, 75, 100], ["TSS", "25%", "50%", "75%", "TES"])
    cbar = fig.colorbar(image, ax=ax, pad=0.02)
    cbar.set_label(color_label)
    ax.set_facecolor("#EEEEEE")
    fig.tight_layout()
    save_figure(fig, outdir / "TSS_to_TES_metagene_heatmap.pdf")
    plt.close(fig)
    selected[["gene_id", "Gene", "chr", "start", "end", "strand"]].to_csv(
        outdir / "TSS_to_TES_heatmap_row_metadata.csv", index=False
    )


def plot_tss(cfg, outdir, genes):
    selected, matrices = tss_matrices(cfg["tracks"], genes)
    group_matrices = {name: group_mean(matrices, members) for name, members in cfg["groups"].items()}
    group_names = list(group_matrices)
    if len(group_names) == 2:
        heat = group_matrices[group_names[1]] - group_matrices[group_names[0]]
        color_label = "delta log2(ChIP/Input)\nKD − Scramble"
        row_label = "Protein-coding gene TSSs\nranked by |KD − Scramble|"
    else:
        heat = group_matrices[group_names[0]]
        color_label = "log2(p53 ChIP/Input)\ndm6-normalized BAMs"
        row_label = "Protein-coding gene TSSs\nranked by |signal|"
    score = np.nanmean(np.abs(heat), axis=1)
    keep = np.isfinite(score)
    order = np.argsort(score[keep])[::-1][:5000]
    heat = np.nan_to_num(heat[keep][order])
    selected = selected.loc[keep].iloc[order]
    vmax = np.nanpercentile(np.abs(heat), 99)
    vmax = vmax if np.isfinite(vmax) and vmax > 0 else 1

    fig, ax = plt.subplots(figsize=(8.5, 8))
    image = ax.imshow(heat, aspect="auto", interpolation="nearest", cmap="vlag",
                      vmin=-vmax, vmax=vmax, extent=(-2000, 2000, heat.shape[0], 0))
    ax.set(title=f"{cfg['title']}: TSS-centered signal", xlabel="Distance from TSS (bp)", ylabel=row_label)
    ax.set_xticks([-2000, -1000, 0, 1000, 2000], ["−2 kb", "−1 kb", "TSS", "+1 kb", "+2 kb"])
    cbar = fig.colorbar(image, ax=ax, pad=0.02); cbar.set_label(color_label)
    fig.tight_layout(); save_figure(fig, outdir / "TSS_centered_heatmap.pdf"); plt.close(fig)

    x = np.linspace(-2000, 2000, heat.shape[1])
    profile_table = pd.DataFrame({"distance_bp": x})
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for group, matrix in group_matrices.items():
        profile = np.nanmean(matrix, axis=0)
        profile_table[group] = profile
        ax.plot(x, profile, lw=2.4, label=group, color=cfg["colors"][group])
    ax.axvline(0, color="0.45", ls="--", lw=1)
    ax.set_xlim(-2000, 2000)
    ax.set(title=f"{cfg['title']}: average TSS profile", xlabel="Distance from TSS (bp)",
           ylabel="Mean log2(ChIP/Input) signal")
    ax.legend(frameon=False); sns.despine(ax=ax)
    fig.tight_layout(); save_figure(fig, outdir / "TSS_average_profile.pdf"); plt.close(fig)
    profile_table.to_csv(outdir / "TSS_average_profiles.csv", index=False)
    selected[["gene_id", "Gene", "chr", "tss"]].to_csv(outdir / "TSS_heatmap_row_metadata.csv", index=False)


def plot_annotation(outdir, cfg):
    annotation = pd.read_csv(outdir / "peak_annotation_summary.csv")
    category_order = annotation.groupby("Category")["Percentage"].sum().sort_values(ascending=True).index
    height = max(5.5, 0.35 * len(category_order) + 2)
    fig, ax = plt.subplots(figsize=(9, height))
    sns.barplot(data=annotation, y="Category", x="Percentage", hue="Sample",
                order=category_order, palette="colorblind", ax=ax)
    ax.set(title=f"{cfg['title']}: detailed peak annotation", xlabel="Peaks (%)", ylabel=None)
    ax.legend(title=None, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    sns.despine(ax=ax); fig.tight_layout()
    fig.savefig(outdir / "peak_annotation_summary_barplot.pdf"); fig.savefig(outdir / "peak_annotation_summary_barplot.png", dpi=150); plt.close(fig)


def plot_diff(diff_path, outdir, cfg):
    if not diff_path or not (root / diff_path).is_file():
        (outdir / "differential_plot_status.txt").write_text(
            "No between-condition differential result exists; descriptive plots only.\n")
        return
    df = pd.read_csv(root / diff_path)
    fold = next((x for x in ["Fold", "log2Fold", "log2FC"] if x in df), None)
    fdr = next((x for x in ["FDR", "fdr", "padj"] if x in df), None)
    mean = next((x for x in ["Conc", "conc", "score"] if x in df), None)
    if not fold or not fdr: return
    df["log2FC"] = pd.to_numeric(df[fold], errors="coerce")
    df["FDR_value"] = pd.to_numeric(df[fdr], errors="coerce")
    df["mean_signal"] = pd.to_numeric(df[mean], errors="coerce") if mean else np.nan
    df["neglog10FDR"] = -np.log10(df["FDR_value"].clip(lower=np.finfo(float).tiny))
    df["Status"] = np.select(
        [(df.FDR_value < .05) & (df.log2FC > 0), (df.FDR_value < .05) & (df.log2FC < 0)],
        ["Significant gain", "Significant loss"], default="Not significant")
    palette = {"Significant gain": gain_color, "Significant loss": loss_color, "Not significant": "#BDBDBD"}
    df.to_csv(outdir / "differential_plot_data.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for status in ["Not significant", "Significant loss", "Significant gain"]:
        part = df[df.Status == status]
        ax.scatter(part.mean_signal, part.log2FC, s=12, alpha=.55, c=palette[status], label=status, edgecolors="none")
    ax.axhline(0, ls="--", color="0.35", lw=1)
    ax.set_ylim(-symmetric_bound(df.log2FC), symmetric_bound(df.log2FC))
    ax.set(title=f"{cfg['title']}: MA plot", xlabel="Mean normalized binding concentration (DiffBind Conc)", ylabel="log2 fold-change")
    ax.legend(frameon=False, markerscale=1.5); sns.despine(ax=ax)
    fig.tight_layout(); save_figure(fig, outdir / "MA_plot.pdf"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for status in ["Not significant", "Significant loss", "Significant gain"]:
        part = df[df.Status == status]
        ax.scatter(part.log2FC, part.neglog10FDR, s=12, alpha=.55, c=palette[status], label=status, edgecolors="none")
    ax.set(title=f"{cfg['title']}: differential binding", xlabel="log2 fold-change", ylabel="−log10(FDR)")
    ax.set_xlim(-symmetric_bound(df.log2FC), symmetric_bound(df.log2FC))
    ax.legend(frameon=False, markerscale=1.5); sns.despine(ax=ax)
    fig.tight_layout(); save_figure(fig, outdir / "volcano_plot.pdf"); plt.close(fig)

    significant = df[df.FDR_value < .05].copy()
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    sns.histplot(data=significant, x="log2FC", hue="Status", bins=50, element="step",
                 palette=palette, common_norm=False, ax=ax)
    ax.axvline(0, ls="--", color="0.35", lw=1)
    ax.set_xlim(-symmetric_bound(significant.log2FC), symmetric_bound(significant.log2FC))
    ax.set(title=f"{cfg['title']}: significant peak effect sizes",
           xlabel="DiffBind log2 fold-change (FDR < 0.05)", ylabel="Peaks")
    sns.despine(ax=ax); fig.tight_layout()
    fig.savefig(outdir / "peak_intensity_log2FC_distribution.pdf"); fig.savefig(outdir / "peak_intensity_log2FC_distribution.png", dpi=150); plt.close(fig)


def main():
    requested = sys.argv[1] if len(sys.argv) > 1 else "all"
    if requested != "all" and requested not in datasets:
        raise SystemExit(f"Unknown dataset: {requested}; choose {', '.join(datasets)}")
    selected = datasets if requested == "all" else {requested: datasets[requested]}
    genes = load_genes()
    for dataset, cfg in selected.items():
        outdir = out_root / dataset; outdir.mkdir(parents=True, exist_ok=True)
        for path in cfg["tracks"].values():
            if not (root / path).is_file(): raise FileNotFoundError(path)
        plot_tss(cfg, outdir, genes)
        plot_tss_to_tes(cfg, outdir, genes)
        plot_annotation(outdir, cfg)
        plot_diff(cfg["diff"], outdir, cfg)
        (outdir / "figure_notes.txt").write_text(
            "Signal tracks are log2(ChIP/Input), generated from dm6 spike-in-normalized BAMs.\n"
            "TSS-to-TES metaplots show the mean gene-body signal with a SEM ribbon.\n"
            "TSS-to-TES heatmap rows are the 3,000 highest-absolute-signal protein-coding genes, not peaks or differential genes.\n"
            "Differential plots use DiffBind-normalized binding results where available.\n")
        make_master_pdf(dataset, outdir)
        print(f"Completed {cfg['title']}", flush=True)


if __name__ == "__main__":
    main()
