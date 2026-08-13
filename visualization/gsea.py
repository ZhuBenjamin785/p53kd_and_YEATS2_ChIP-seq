
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import gseapy as gp


projectdir = Path(__file__).resolve().parent
defaultresults = projectdir / "output_files" / "significant_results.csv"
defaultoutdir = projectdir / "gsea_out"
defaultgeneset = "KEGG_2016"


def build_ranked_list(results_csv: Path) -> pd.DataFrame:
    """Load DESeq2 results and return a 2-column ranked gene table."""
    if not results_csv.exists():
        raise FileNotFoundError(f"Missing input file: {results_csv}")

    df = pd.read_csv(results_csv)
    required = {"Geneid", "log2FoldChange"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"{results_csv} is missing required columns: {sorted(missing)}"
        )

    rnk = (
        df.loc[:, ["Geneid", "log2FoldChange"]]
        .rename(columns={"Geneid": "gene", "log2FoldChange": "score"})
        .assign(gene=lambda frame: frame["gene"].astype(str).str.upper())
        .dropna()
        .drop_duplicates(subset="gene")
        .sort_values("score", ascending=False)
        .reset_index(drop=True)
    )



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
    ranked = build_ranked_list(args.results)
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
