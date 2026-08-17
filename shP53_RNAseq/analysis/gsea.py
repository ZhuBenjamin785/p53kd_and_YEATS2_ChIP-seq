
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import gseapy as gp


projectdir = Path(__file__).resolve().parents[2]
defaultresults = projectdir / "shared" / "rna_seq_dea" / "shp53_vs_shLacZ_0hr" / "significant_results.csv"
defaultoutdir = projectdir / "shared" / "rna_seq_dea" / "gsea_out"
defaultgeneset = "KEGG_2016"


def build_ranked_list(results_csv: Path, padj: float = 0.05, min_log2fc: float = 1.0) -> pd.DataFrame:
    """Load and rank only significant DESeq2 results."""
    if not results_csv.exists():
        raise FileNotFoundError(f"Missing input file: {results_csv}")

    df = pd.read_csv(results_csv)
    gene_column = next(
        (column for column in ("gene_name", "Geneid", "gene", "symbol", "Unnamed: 0")
         if column in df.columns),
        None,
    )
    required = {"log2FoldChange", "padj"}
    missing = required - set(df.columns)
    if gene_column is None:
        missing.add("gene_name/Geneid")
    if missing:
        raise ValueError(f"{results_csv} is missing required columns: {sorted(missing)}")

    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    df = df.loc[
        (df["padj"] <= padj)
        & (df["log2FoldChange"].abs() >= min_log2fc)
    ].copy()
    if df.empty:
        raise ValueError("No significant genes passed the padj/log2FC filters.")

    rnk = (
        df.loc[:, [gene_column, "log2FoldChange"]]
        .rename(columns={gene_column: "gene", "log2FoldChange": "score"})
        .assign(gene=lambda frame: frame["gene"].astype(str).str.upper())
        .dropna()
        .drop_duplicates(subset="gene")
        .sort_values("score", ascending=False)
        .reset_index(drop=True)
    )
    return rnk



def run_prerank(
    rnk: pd.DataFrame,
    gene_sets: str,
    outdir: Path,
    permutation_num: int = 1000,
    min_size: int = 1,
    max_size: int = 1000,
    seed: int = 6,
) -> gp.prerank:
    """Run GSEApy preranked analysis and return the result object."""
    outdir.mkdir(parents=True, exist_ok=True)

    return gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        permutation_num=permutation_num,
        outdir=str(outdir),
        seed=seed,
        min_size=min_size,
        max_size=max_size,
        method="multilevel",
        verbose=True,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run preranked GSEA on DESeq2 results with gseapy."
    )
    parser.add_argument(
        "--results",
        type=Path,
        default=defaultresults,
        help="Path to DESeq2 results.csv",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=defaultoutdir,
        help="Output directory for GSEA results",
    )
    parser.add_argument(
        "--gene-sets",
        default=defaultgeneset,
        help="Enrichr gene set name or path to a .gmt file",
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=1000,
        help="Number of permutations",
    )
    parser.add_argument("--padj", type=float, default=0.05)
    parser.add_argument("--min-log2fc", type=float, default=1.0)
    parser.add_argument(
        "--min-size",
        type=int,
        default=1,
        help="Minimum gene set size",
    )
    parser.add_argument(
        "--max-size",
        type=int,
        default=1000,
        help="Maximum gene set size",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=6,
        help="Random seed",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ranked = build_ranked_list(args.results, args.padj, args.min_log2fc)
    pre_res = run_prerank(
        ranked,
        gene_sets=args.gene_sets,
        outdir=args.outdir,
        permutation_num=args.permutations,
        min_size=args.min_size,
        max_size=args.max_size,
        seed=args.seed,
    )

    print(pre_res.res2d.head())
    print(f"Saved GSEA output to {args.outdir}")


if __name__ == "__main__":
    main()
