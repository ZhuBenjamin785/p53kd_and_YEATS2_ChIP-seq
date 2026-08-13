#!/usr/bin/env python3
                      
"""Functional enrichment analysis for RNA-seq DESeq2 results."""

from pathlib import Path
import argparse

import pandas as pd


project_dir = Path(__file__).resolve().parent
default_results = project_dir / "rna_seq_dea" / "results.csv"
default_output = project_dir / "rna_seq_fea"


def choose_gene_column(df):
    for column in ("gene_name", "Geneid", "gene", "symbol"):
        if column in df.columns:
            return column
    raise ValueError("DE results must contain gene_name, Geneid, gene, or symbol.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=default_results)
    parser.add_argument("--output-dir", type=Path, default=default_output)
    parser.add_argument("--padj", type=float, default=0.05)
    parser.add_argument("--min-log2fc", type=float, default=1.0)
    parser.add_argument("--direction", choices=("up", "down", "both"), default="both")
    args = parser.parse_args()

    if not args.results.exists():
        raise SystemExit("DESeq2 results not found: {}".format(args.results))
    df = pd.read_csv(args.results)
    required = {"padj", "log2FoldChange"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit("DE results missing required columns: {}".format(", ".join(sorted(missing))))
    gene_column = choose_gene_column(df)
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    significant = df.loc[
        (df["padj"] <= args.padj)
        & (df["log2FoldChange"].abs() >= args.min_log2fc)
    ].copy()
    if args.direction == "up":
        significant = significant[significant["log2FoldChange"] > 0]
    elif args.direction == "down":
        significant = significant[significant["log2FoldChange"] < 0]
    genes = (
        significant[gene_column]
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda values: ~values.isin({"", "NA", "nan", "None"})]
        .drop_duplicates()
        .tolist()
    )
    if not genes:
        raise SystemExit("No genes passed the requested padj/log2FC/direction filters.")

    try:
        from gseapy import barplot, dotplot, enrichr
    except ImportError as exc:
        raise SystemExit("gseapy is unavailable in the active environment.") from exc

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"gene": genes}).to_csv(args.output_dir / "genes_used.tsv", sep="\t", index=False)
    enrichment = enrichr(
        gene_list=genes,
        gene_sets=["KEGG_2021_Human", "MSigDB_Hallmark_2020"],
        organism="human",
        outdir=str(args.output_dir / "enrichr"),
    )
    results = getattr(enrichment, "results", None)
    if results is None or results.empty:
        raise SystemExit("Enrichr returned no enrichment rows.")
    results.to_csv(args.output_dir / "enrichment_results.csv", index=False)
    dotplot(
        results, column="Adjusted P-value", x="Gene_set", size=5, top_term=10,
        figsize=(6, 5), title="RNA-seq functional enrichment", xticklabels_rot=45,
        show_ring=True, marker="o", ofname=str(args.output_dir / "enrichment_dotplot.png"),
    )
    barplot(
        results, column="Adjusted P-value", group="Gene_set", size=5, top_term=10,
        figsize=(8, 5), ofname=str(args.output_dir / "enrichment_barplot.png"),
    )
    print("Enrichment used {} genes; results written to {}".format(len(genes), args.output_dir))


if __name__ == "__main__":
    main()
