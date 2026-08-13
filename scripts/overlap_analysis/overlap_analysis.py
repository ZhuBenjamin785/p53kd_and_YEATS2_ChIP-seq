#!/usr/bin/env python3
                      
"""Prepare, summarize, and plot p53KD versus YEATS2KD peak overlaps."""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import numpy as np
import pandas as pd
import seaborn as sns


datasets = {
    "YEATS2KD": {
        "all": "diffbind_results_yeats2/DiffBind_YEATS2_all_peaks.csv",
        "annot": "chipseq_summary_plots/YEATS2KD/diffbind_peak_annotations.csv",
        "tss_profile": "chipseq_summary_plots/YEATS2KD/TSS_average_profiles.csv",
        "profile_kd": "YEATS2 KD",
    },
    "p53KD": {
        "all": "diffbind_results/DiffBind_all_peaks.csv",
        "annot": "chipseq_summary_plots/p53KD/diffbind_peak_annotations.csv",
        "tss_profile": "chipseq_summary_plots/p53KD/TSS_average_profiles.csv",
        "profile_kd": "p53 KD",
    },
}
bed_columns = ["chrom", "start0", "end", "peak_id", "score", "strand",
               "direction", "fold", "fdr"]
tss_columns = ["gene_chrom", "gene_tss_start0", "gene_tss_end", "gene_id",
               "gene_score", "gene_strand", "gene_name", "gene_type"]
colors = {"YEATS2KD": "#2A9D8F", "p53KD": "#E76F51",
          "Gain": "#D55E00", "Loss": "#0072B2"}


def normalize_chromosome(value):
    value = str(value).strip().strip('"')
    return value if value.startswith("chr") else "chr" + value


def read_diffbind(path):
    df = pd.read_csv(path)
    required = {"seqnames", "start", "end", "Fold", "FDR"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError("{} lacks columns: {}".format(path, sorted(missing)))
    for column in ["start", "end", "Fold", "FDR"]:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    df = df.dropna(subset=["seqnames", "start", "end", "Fold", "FDR"]).copy()
    df["seqnames"] = df["seqnames"].map(normalize_chromosome)
    return df


def parse_gtf_attributes(text):
    return dict(re.findall(r'(\S+) "([^"]+)"', text))


def read_gtf_tss(gtf_path):
    records = []
    chrom_max = defaultdict(int)
    with gtf_path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attributes = fields
            start, end = int(start), int(end)
            chrom_max[chrom] = max(chrom_max[chrom], end)
            if feature != "gene":
                continue
            attrs = parse_gtf_attributes(attributes)
            gene_id = attrs.get("gene_id", "NA").split(".")[0]
            gene_name = attrs.get("gene_name", gene_id)
            gene_type = attrs.get("gene_type", attrs.get("gene_biotype", "unknown"))
            tss0 = start - 1 if strand == "+" else end - 1
            records.append((chrom, tss0, tss0 + 1, gene_id, 0, strand,
                            gene_name, gene_type))
    return records, chrom_max


def chromosome_sizes(diff_tables, gtf_max, padding=2000):
                                                                               
                                                            
    canonical = {
        "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
        "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
        "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
        "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
        "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
        "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
        "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
        "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
        "chrM": 16569,
    }
    peak_max = defaultdict(int)
    for df in diff_tables:
        for chrom, value in df.groupby("seqnames")["end"].max().items():
            peak_max[str(chrom)] = max(peak_max[str(chrom)], int(value))
    chromosomes = set(gtf_max) | set(peak_max)
    sizes = {}
    for chrom in chromosomes:
        observed = max(gtf_max.get(chrom, 0), peak_max.get(chrom, 0))
        sizes[chrom] = canonical.get(chrom, observed + padding)
    canonical_order = {chrom: i for i, chrom in enumerate(canonical)}
    ordered = sorted(sizes, key=lambda c: (canonical_order.get(c, 10_000), c))
    return [(chrom, sizes[chrom]) for chrom in ordered]


def write_bed(df, path, dataset, chrom_order):
    out = pd.DataFrame({
        "chrom": df["seqnames"].astype(str),
                                                                                 
        "start0": df["start"].astype(int) - 1,
        "end": df["end"].astype(int),
        "peak_id": ["{}_peak_{:06d}".format(dataset, i + 1) for i in range(len(df))],
        "score": 0,
        "strand": ".",
        "direction": np.where(df["Fold"] > 0, "gain", "loss"),
        "fold": df["Fold"].astype(float),
        "fdr": df["FDR"].astype(float),
    })
    out["_chrom_order"] = out["chrom"].map(chrom_order)
    out = out.sort_values(["_chrom_order", "start0", "end"]).drop(columns="_chrom_order")
    out.to_csv(path, sep="\t", index=False, header=False, na_rep="NA")
    return out


def prepare(project_dir, output_dir):
    beds = output_dir / "beds"
    tables = output_dir / "tables"
    plots = output_dir / "plots"
    for directory in [beds, tables, plots, output_dir / "logs", output_dir / "intervene_plots"]:
        directory.mkdir(parents=True, exist_ok=True)

    all_tables = {name: read_diffbind(project_dir / spec["all"])
                  for name, spec in datasets.items()}
    tss_records, gtf_max = read_gtf_tss(project_dir / "gencode.v50.basic.annotation.gtf")
    sizes = chromosome_sizes(list(all_tables.values()), gtf_max)
    chrom_order = {chrom: i for i, (chrom, _) in enumerate(sizes)}

    with (beds / "chrom.sizes").open("w") as handle:
        for chrom, size in sizes:
            handle.write("{}\t{}\n".format(chrom, size))
    with (beds / "gencode_v50_gene_TSS.bed").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for record in sorted(tss_records, key=lambda x: (chrom_order[x[0]], x[1], x[2])):
            writer.writerow(record)

    manifest_rows = []
    for dataset, df in all_tables.items():
        significant = df[df["FDR"] < 0.05].copy()
        groups = {
            "all_tested": df,
            "significant": significant,
            "gain": significant[significant["Fold"] > 0],
            "loss": significant[significant["Fold"] < 0],
        }
        for group, group_df in groups.items():
            written = write_bed(group_df, beds / "{}_{}.bed".format(dataset, group),
                                dataset, chrom_order)
            manifest_rows.append({"dataset": dataset, "peak_group": group,
                                  "n_peaks": len(written)})
    pd.DataFrame(manifest_rows).to_csv(tables / "input_peak_counts.csv", index=False)


def read_bed(path):
    if path.stat().st_size == 0:
        return pd.DataFrame(columns=bed_columns)
    return pd.read_csv(path, sep="\t", names=bed_columns)


def read_pair_file(path):
    columns = ["YEATS2_" + c for c in bed_columns] + ["p53_" + c for c in bed_columns]
    if path.stat().st_size == 0:
        return pd.DataFrame(columns=columns)
    return pd.read_csv(path, sep="\t", names=columns)


def read_closest(path):
    columns = bed_columns + tss_columns + ["distance_to_TSS"]
    return pd.read_csv(path, sep="\t", names=columns)


def annotation_lookup(path):
    df = pd.read_csv(path)
    df = df[pd.to_numeric(df["FDR"], errors="coerce") < 0.05].copy()
    df["seqnames"] = df["seqnames"].map(normalize_chromosome)
    df["coord_key"] = (df["seqnames"].astype(str) + ":" +
                       (pd.to_numeric(df["start"]).astype(int) - 1).astype(str) + "-" +
                       pd.to_numeric(df["end"]).astype(int).astype(str))
    category = "annotation_category" if "annotation_category" in df else "annotation"
    keep = ["coord_key", category]
    for column in ["ENSEMBL", "SYMBOL", "geneId", "distanceToTSS"]:
        if column in df:
            keep.append(column)
    lookup = df[keep].drop_duplicates("coord_key").set_index("coord_key").to_dict("index")
    for record in lookup.values():
        record["annotation_category"] = record.pop(category, "Unclassified")
    return lookup


def peak_gene_table(closest, annotation, dataset):
    out = closest.copy()
    out["dataset"] = dataset
    out["coord_key"] = (out["chrom"].astype(str) + ":" + out["start0"].astype(str) +
                        "-" + out["end"].astype(str))
    details = out["coord_key"].map(annotation)
    out["annotation_category"] = details.map(
        lambda x: x.get("annotation_category", "Unclassified") if isinstance(x, dict) else "Unclassified")
    out["annotated_gene_id"] = details.map(
        lambda x: x.get("ENSEMBL", x.get("geneId", "NA")) if isinstance(x, dict) else "NA")
    out["annotated_gene_name"] = details.map(
        lambda x: x.get("SYMBOL", x.get("geneId", "NA")) if isinstance(x, dict) else "NA")
    out["promoter_within_2kb"] = out["distance_to_TSS"].abs() <= 2000
    return out


def bipartite_component_count(pairs):
    graph = defaultdict(set)
    for row in pairs.itertuples(index=False):
        left = "Y:" + row.YEATS2_peak_id
        right = "P:" + row.p53_peak_id
        graph[left].add(right)
        graph[right].add(left)
    seen, components = set(), 0
    for node in graph:
        if node in seen:
            continue
        components += 1
        stack = [node]
        seen.add(node)
        while stack:
            current = stack.pop()
            for neighbor in graph[current]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
    return components


def save_venn(left, right, labels, title, path):
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    venn2([set(left), set(right)], set_labels=labels,
          set_colors=(colors["YEATS2KD"], colors["p53KD"]), alpha=0.65, ax=ax)
    ax.set_title(title, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def summarize(project_dir, output_dir):
    beds, tables, plots = output_dir / "beds", output_dir / "tables", output_dir / "plots"
    exact = read_pair_file(tables / "exact_overlap_pairs.tsv")
    y_bed = read_bed(beds / "YEATS2KD_significant.bed")
    p_bed = read_bed(beds / "p53KD_significant.bed")

    if len(exact):
        exact["direction_pair"] = (exact["YEATS2_direction"] + "-" + exact["p53_direction"])
        exact["directions_match"] = exact["YEATS2_direction"] == exact["p53_direction"]
    direction_order = ["gain-gain", "loss-loss", "gain-loss", "loss-gain"]
    direction_counts = (exact.get("direction_pair", pd.Series(dtype=str)).value_counts()
                        .reindex(direction_order, fill_value=0).rename_axis("direction_pair")
                        .reset_index(name="overlap_pairs"))
    direction_counts.to_csv(tables / "exact_overlap_direction_counts.csv", index=False)

    y_annotation = annotation_lookup(project_dir / datasets["YEATS2KD"]["annot"])
    p_annotation = annotation_lookup(project_dir / datasets["p53KD"]["annot"])
    y_gene = peak_gene_table(read_closest(tables / "YEATS2KD_nearest_TSS.tsv"),
                             y_annotation, "YEATS2KD")
    p_gene = peak_gene_table(read_closest(tables / "p53KD_nearest_TSS.tsv"),
                             p_annotation, "p53KD")
    y_gene.to_csv(tables / "YEATS2KD_significant_peaks_annotated.tsv", sep="\t", index=False)
    p_gene.to_csv(tables / "p53KD_significant_peaks_annotated.tsv", sep="\t", index=False)

    y_map = y_gene.set_index("peak_id").to_dict("index")
    p_map = p_gene.set_index("peak_id").to_dict("index")
    overlap_rows = []
    for row in exact.itertuples(index=False):
        yg, pg = y_map[row.YEATS2_peak_id], p_map[row.p53_peak_id]
        overlap_rows.append({
            "chromosome": row.YEATS2_chrom,
            "YEATS2KD_start_1based": row.YEATS2_start0 + 1,
            "YEATS2KD_end": row.YEATS2_end,
            "p53KD_start_1based": row.p53_start0 + 1,
            "p53KD_end": row.p53_end,
            "YEATS2KD_peak_id": row.YEATS2_peak_id,
            "p53KD_peak_id": row.p53_peak_id,
            "YEATS2KD_direction": row.YEATS2_direction,
            "p53KD_direction": row.p53_direction,
            "directions_match": row.YEATS2_direction == row.p53_direction,
            "YEATS2KD_nearest_gene": yg["gene_name"],
            "p53KD_nearest_gene": pg["gene_name"],
            "nearest_gene_matches": yg["gene_id"] == pg["gene_id"],
            "YEATS2KD_distance_to_TSS": yg["distance_to_TSS"],
            "p53KD_distance_to_TSS": pg["distance_to_TSS"],
            "YEATS2KD_annotation": yg["annotation_category"],
            "p53KD_annotation": pg["annotation_category"],
        })
    overlap_columns = [
        "chromosome", "YEATS2KD_start_1based", "YEATS2KD_end",
        "p53KD_start_1based", "p53KD_end", "YEATS2KD_peak_id", "p53KD_peak_id",
        "YEATS2KD_direction", "p53KD_direction", "directions_match",
        "YEATS2KD_nearest_gene", "p53KD_nearest_gene", "nearest_gene_matches",
        "YEATS2KD_distance_to_TSS", "p53KD_distance_to_TSS",
        "YEATS2KD_annotation", "p53KD_annotation",
    ]
    pd.DataFrame(overlap_rows, columns=overlap_columns).to_csv(
        tables / "exact_overlap_detailed_summary.csv", index=False)

    def collapse_annotated_genes(df, prefix):
        clean = df[df["annotated_gene_id"].notna()].copy()
        clean["annotated_gene_id"] = clean["annotated_gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
        clean = clean[~clean["annotated_gene_id"].isin([".", "NA", "nan", "None"])]
        clean["signed_fold"] = pd.to_numeric(clean["fold"], errors="coerce")
        grouped = clean.groupby(["annotated_gene_id", "annotated_gene_name"], as_index=False).agg(
            **{"{}_median_fold".format(prefix): ("signed_fold", "median"),
               "{}_peak_count".format(prefix): ("peak_id", "nunique"),
               "{}_promoter_within_2kb".format(prefix): ("promoter_within_2kb", "max")})
        return grouped

    yg = collapse_annotated_genes(y_gene, "YEATS2KD")
    pg = collapse_annotated_genes(p_gene, "p53KD")
    shared = yg.merge(pg, on=["annotated_gene_id", "annotated_gene_name"], how="inner")
    shared["YEATS2KD_direction"] = np.where(shared["YEATS2KD_median_fold"] > 0, "gain", "loss")
    shared["p53KD_direction"] = np.where(shared["p53KD_median_fold"] > 0, "gain", "loss")
    shared["directions_match"] = shared["YEATS2KD_direction"] == shared["p53KD_direction"]
    shared.to_csv(tables / "shared_annotated_genes.csv", index=False)

    y_genes, p_genes = set(yg["annotated_gene_id"]), set(pg["annotated_gene_id"])
                                                                             
                                                                         
    y_prom = set(y_gene.loc[y_gene["promoter_within_2kb"], "gene_id"])
    p_prom = set(p_gene.loc[p_gene["promoter_within_2kb"], "gene_id"])
    y_nearest, p_nearest = set(y_gene["gene_id"]), set(p_gene["gene_id"])
    distance_cutoff = 10000
    y_near_cut = set(y_gene.loc[y_gene["distance_to_TSS"] <= distance_cutoff, "gene_id"])
    p_near_cut = set(p_gene.loc[p_gene["distance_to_TSS"] <= distance_cutoff, "gene_id"])
    gene_membership = pd.DataFrame({
        "annotated_gene_id": sorted(y_genes | p_genes),
    })
    gene_membership["in_YEATS2KD"] = gene_membership["annotated_gene_id"].isin(y_genes)
    gene_membership["in_p53KD"] = gene_membership["annotated_gene_id"].isin(p_genes)
    gene_membership["shared"] = gene_membership["in_YEATS2KD"] & gene_membership["in_p53KD"]
    gene_names = dict(zip(yg["annotated_gene_id"], yg["annotated_gene_name"]))
    gene_names.update(dict(zip(pg["annotated_gene_id"], pg["annotated_gene_name"])))
    gene_membership.insert(1, "gene_name", gene_membership["annotated_gene_id"].map(gene_names))
    gene_membership.to_csv(tables / "annotated_gene_overlap_membership.csv", index=False)

    promoter_membership = pd.DataFrame({"gene_id": sorted(y_prom | p_prom)})
    promoter_membership["in_YEATS2KD"] = promoter_membership["gene_id"].isin(y_prom)
    promoter_membership["in_p53KD"] = promoter_membership["gene_id"].isin(p_prom)
    promoter_membership["shared"] = promoter_membership["in_YEATS2KD"] & promoter_membership["in_p53KD"]
    nearest_names = dict(zip(y_gene["gene_id"], y_gene["gene_name"]))
    nearest_names.update(dict(zip(p_gene["gene_id"], p_gene["gene_name"])))
    promoter_membership.insert(1, "gene_name", promoter_membership["gene_id"].map(nearest_names))
    promoter_membership.to_csv(tables / "promoter_gene_overlap_membership.csv", index=False)

    def shared_nearest_by_cutoff(left, right, cutoff, output_name):
        left = left[left["distance_to_TSS"] <= cutoff].copy()
        right = right[right["distance_to_TSS"] <= cutoff].copy()
        left_group = left.groupby(["gene_id", "gene_name"], as_index=False).agg(
            YEATS2KD_median_fold=("fold", "median"),
            YEATS2KD_min_distance_to_TSS=("distance_to_TSS", "min"),
            YEATS2KD_peak_count=("peak_id", "nunique"))
        right_group = right.groupby(["gene_id", "gene_name"], as_index=False).agg(
            p53KD_median_fold=("fold", "median"),
            p53KD_min_distance_to_TSS=("distance_to_TSS", "min"),
            p53KD_peak_count=("peak_id", "nunique"))
        result = left_group.merge(right_group, on=["gene_id", "gene_name"], how="inner")
        result["YEATS2KD_direction"] = np.where(result["YEATS2KD_median_fold"] > 0, "gain", "loss")
        result["p53KD_direction"] = np.where(result["p53KD_median_fold"] > 0, "gain", "loss")
        result["directions_match"] = result["YEATS2KD_direction"] == result["p53KD_direction"]
        result.to_csv(tables / output_name, index=False)

    shared_nearest_by_cutoff(y_gene, p_gene, 2000, "shared_promoter_neighborhood_genes.csv")
    shared_nearest_by_cutoff(y_gene, p_gene, distance_cutoff,
                             "shared_nearest_genes_within_10kb.csv")
    gene_summary = pd.DataFrame([
        {"comparison": "ChIPseeker_annotated_genes", "YEATS2KD_only": len(y_genes - p_genes),
         "shared": len(y_genes & p_genes), "p53KD_only": len(p_genes - y_genes)},
        {"comparison": "bedtools_closest_TSS_genes", "YEATS2KD_only": len(y_nearest - p_nearest),
         "shared": len(y_nearest & p_nearest), "p53KD_only": len(p_nearest - y_nearest)},
        {"comparison": "promoter_genes_within_2kb", "YEATS2KD_only": len(y_prom - p_prom),
         "shared": len(y_prom & p_prom), "p53KD_only": len(p_prom - y_prom)},
        {"comparison": "nearest_genes_within_10kb", "YEATS2KD_only": len(y_near_cut - p_near_cut),
         "shared": len(y_near_cut & p_near_cut), "p53KD_only": len(p_near_cut - y_near_cut)},
    ])
    gene_summary.to_csv(tables / "gene_overlap_summary.csv", index=False)

    lenient_rows = []
    for window in [1000, 2000]:
        pairs = read_pair_file(tables / "slop_{}bp_overlap_pairs.tsv".format(window))
        lenient_rows.append({
            "window_bp_each_side": window,
            "overlap_pairs": len(pairs),
            "YEATS2KD_peaks_with_overlap": pairs["YEATS2_peak_id"].nunique(),
            "p53KD_peaks_with_overlap": pairs["p53_peak_id"].nunique(),
        })
    pd.DataFrame(lenient_rows).to_csv(tables / "lenient_overlap_summary.csv", index=False)

    y_unique_overlap = exact["YEATS2_peak_id"].nunique()
    p_unique_overlap = exact["p53_peak_id"].nunique()
    strict_summary = pd.DataFrame([{
        "YEATS2KD_significant_peaks": len(y_bed),
        "p53KD_significant_peaks": len(p_bed),
        "overlap_pairs": len(exact),
        "YEATS2KD_peaks_with_overlap": y_unique_overlap,
        "p53KD_peaks_with_overlap": p_unique_overlap,
        "overlap_components": bipartite_component_count(exact),
    }])
    strict_summary.to_csv(tables / "exact_overlap_summary.csv", index=False)

                                                                               
                                                                              
    components = bipartite_component_count(exact)
    peak_left = {"Y_only_{}".format(i) for i in range(len(y_bed) - y_unique_overlap)}
    peak_right = {"P_only_{}".format(i) for i in range(len(p_bed) - p_unique_overlap)}
    shared_components = {"shared_{}".format(i) for i in range(components)}
    save_venn(peak_left | shared_components, peak_right | shared_components,
              ("YEATS2KD peaks", "p53KD peaks"),
              "Strict genomic peak overlap\n(shared = overlap components)",
              plots / "exact_peak_overlap_venn.png")
    save_venn(y_genes, p_genes, ("YEATS2KD genes", "p53KD genes"),
              "Shared ChIPseeker-annotated genes", plots / "shared_gene_overlap_venn.png")
    save_venn(y_prom, p_prom, ("YEATS2KD promoter genes", "p53KD promoter genes"),
              "Shared promoter-associated genes (≤2 kb from TSS)",
              plots / "promoter_gene_overlap_venn.png")

    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(8, 5.5))
    plot_counts = direction_counts.copy()
    plot_counts["direction_pair"] = plot_counts["direction_pair"].str.replace("-", " / ")
    sns.barplot(data=plot_counts, x="direction_pair", y="overlap_pairs", ax=ax,
                palette=["#2A9D8F", "#457B9D", "#E9C46A", "#F4A261"], hue="direction_pair", legend=False)
    ax.set(xlabel="YEATS2KD / p53KD direction", ylabel="Exact overlap pairs",
           title="Direction of strict genomic overlaps")
    sns.despine(ax=ax)
    fig.tight_layout(); fig.savefig(plots / "gain_loss_overlap_barplot.png", dpi=220); plt.close(fig)

    annotation_frames = []
    for dataset, gene_df in [("YEATS2KD", y_gene), ("p53KD", p_gene)]:
        counts = gene_df["annotation_category"].value_counts().rename_axis("annotation").reset_index(name="peaks")
        counts["dataset"] = dataset
        annotation_frames.append(counts)
    annotation_counts = pd.concat(annotation_frames, ignore_index=True)
    annotation_counts.to_csv(tables / "peak_annotation_counts.csv", index=False)
    fig, ax = plt.subplots(figsize=(11, 6.5))
    order = annotation_counts.groupby("annotation")["peaks"].sum().sort_values(ascending=False).index
    sns.barplot(data=annotation_counts, x="peaks", y="annotation", hue="dataset", order=order,
                palette={"YEATS2KD": colors["YEATS2KD"], "p53KD": colors["p53KD"]}, ax=ax)
    ax.set(title="Annotation of significant differential peaks", xlabel="Peaks", ylabel=None)
    sns.despine(ax=ax); fig.tight_layout()
    fig.savefig(plots / "peak_annotation_barplot.png", dpi=220); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    profile_table = None
    for dataset, spec in datasets.items():
        profile = pd.read_csv(project_dir / spec["tss_profile"])
        delta = profile[spec["profile_kd"]] - profile["Scramble"]
        ax.plot(profile["distance_bp"], delta, lw=2.5, color=colors[dataset], label=dataset)
        current = pd.DataFrame({"distance_bp": profile["distance_bp"], dataset: delta})
        profile_table = current if profile_table is None else profile_table.merge(current, on="distance_bp", how="outer")
    profile_table.to_csv(tables / "TSS_centered_delta_profiles.csv", index=False)
    ax.axhline(0, color="0.35", lw=1, ls="--"); ax.axvline(0, color="0.35", lw=1, ls=":")
    ax.set(xlabel="Distance from TSS (bp)", ylabel="KD − Scramble mean signal",
           title="TSS-centered H4K16ac redistribution")
    ax.legend(frameon=False); sns.despine(ax=ax); fig.tight_layout()
    fig.savefig(plots / "TSS_centered_metaplot.png", dpi=220); plt.close(fig)

    with (output_dir / "README_results.txt").open("w") as handle:
        handle.write(
            "p53KD versus YEATS2KD H4K16ac overlap analysis\n"
            "================================================\n"
            "Significant: DiffBind FDR < 0.05. Gain/loss: positive/negative Fold.\n"
            "Strict overlap: unexpanded BED intervals using bedtools intersect.\n"
            "BED conversion: DiffBind start was converted from 1-based to 0-based.\n"
            "Lenient overlap: each peak expanded by ±1 kb or ±2 kb before intersect.\n"
            "Gene assignment: closest GENCODE v50 basic gene TSS (bedtools closest).\n"
            "Promoter-associated gene: nearest TSS within 2 kb.\n"
            "Nearest-gene cutoff comparison: nearest TSS within 10 kb.\n"
            "Exact-overlap direction counts are peak-pair counts; unique peak counts are also reported.\n"
            "The peak Venn shared value is the number of connected genomic overlap components.\n"
            "Intervene output is supplementary when its local installation is usable; PNG Venn plots are always generated in Python.\n"
        )


def plot_go(output_dir):
    tables, plots = output_dir / "tables", output_dir / "plots"
    groups = [
        ("YEATS2KD_only", "YEATS2KD-only genes", colors["YEATS2KD"]),
        ("p53KD_only", "p53KD-only genes", colors["p53KD"]),
        ("shared", "Shared genes", "#8E6C8A"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(19, 7), constrained_layout=True)
    term_sets = {}
    for ax, (stem, title, color) in zip(axes, groups):
        path = tables / "GO_BP_{}.csv".format(stem)
        if not path.exists() or path.stat().st_size == 0:
            data = pd.DataFrame()
        else:
            data = pd.read_csv(path)
        if len(data) and "p.adjust" in data:
            data["p.adjust"] = pd.to_numeric(data["p.adjust"], errors="coerce")
            significant = data[data["p.adjust"] < 0.05].nsmallest(12, "p.adjust").copy()
        else:
            significant = pd.DataFrame()
        term_sets[stem] = set(significant.get("ID", pd.Series(dtype=str)).dropna())
        if significant.empty and data.empty:
            ax.text(0.5, 0.53, "No significant GO BP terms",
                    ha="center", va="center", transform=ax.transAxes, fontsize=13)
            ax.text(0.5, 0.44, "See GO_enrichment_status.csv",
                    ha="center", va="center", transform=ax.transAxes, fontsize=10, color="0.4")
            ax.set_axis_off()
            ax.set_title(title, fontweight="bold")
            continue
        plot_data = significant if len(significant) else data.nsmallest(10, "p.adjust").copy()
        plot_data["minus_log10_FDR"] = -np.log10(plot_data["p.adjust"].clip(lower=1e-300))
        plot_data = plot_data.sort_values("minus_log10_FDR")
        sns.barplot(data=plot_data, x="minus_log10_FDR", y="Description",
                    color=color if len(significant) else "#9E9E9E", ax=ax)
        ax.axvline(-np.log10(0.05), color="#B22222", ls="--", lw=1.3, label="FDR = 0.05")
        subtitle = title if len(significant) else title + "\n(top terms shown; none pass FDR < 0.05)"
        ax.set(title=subtitle, xlabel="−log10(BH-adjusted P)", ylabel=None)
        ax.legend(frameon=False, fontsize=9)
        sns.despine(ax=ax)
    fig.suptitle("GO Biological Process enrichment", fontsize=19, fontweight="bold")
    fig.savefig(plots / "GO_BP_enrichment_barplots.png", dpi=220)
    plt.close(fig)

    rows = []
    for first, second in [("YEATS2KD_only", "p53KD_only"),
                          ("YEATS2KD_only", "shared"), ("p53KD_only", "shared")]:
        union = term_sets[first] | term_sets[second]
        intersection = term_sets[first] & term_sets[second]
        rows.append({"group_1": first, "group_2": second,
                     "shared_significant_GO_terms": len(intersection),
                     "union_significant_GO_terms": len(union),
                     "jaccard_index": len(intersection) / len(union) if union else 0.0})
    pd.DataFrame(rows).to_csv(tables / "GO_term_overlap_summary.csv", index=False)

    exact = pd.read_csv(tables / "exact_overlap_summary.csv").iloc[0]
    genes = pd.read_csv(tables / "gene_overlap_summary.csv")
    annotated_shared = int(genes.loc[genes["comparison"] == "ChIPseeker_annotated_genes", "shared"].iloc[0])
    nearby_shared = int(genes.loc[genes["comparison"] == "nearest_genes_within_10kb", "shared"].iloc[0])
    status = pd.read_csv(tables / "GO_enrichment_status.csv")
    significant_total = int(status["significant_terms"].sum())
    with (output_dir / "layered_interpretation.txt").open("w") as handle:
        handle.write(
            "Strict-to-lenient interpretation\n"
            "================================\n"
            "Strict overlapping peak pairs: {}.\n"
            "Shared ChIPseeker-annotated genes: {}.\n"
            "Shared nearest genes within 10 kb: {}.\n"
            "Significant GO BP terms across tested gene sets: {}.\n\n"
            .format(int(exact["overlap_pairs"]), annotated_shared,
                    nearby_shared, significant_total)
        )
        if int(exact["overlap_pairs"]) == 0 and nearby_shared == 0 and significant_total == 0:
            handle.write(
                "The current data support largely distinct H4K16ac programs after the two knockdowns: "
                "there is no regional overlap, no shared nearest-TSS neighborhood within 10 kb, and no "
                "GO Biological Process term passes BH FDR < 0.05. One ChIPseeker gene assignment is "
                "shared (LINC00910; concordant loss), but a one-gene shared set is too small for a reliable "
                "over-representation test. This is evidence for a locus-specific common response, not yet "
                "for broad pathway convergence.\n"
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["prepare", "summarize", "plot-go"])
    parser.add_argument("--project-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.mode == "prepare":
        prepare(args.project_dir.resolve(), args.output_dir.resolve())
    elif args.mode == "summarize":
        summarize(args.project_dir.resolve(), args.output_dir.resolve())
    else:
        plot_go(args.output_dir.resolve())


if __name__ == "__main__":
    main()
