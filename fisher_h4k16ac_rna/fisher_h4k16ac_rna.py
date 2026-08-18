#!/usr/bin/env python3
"""One-sided Fisher tests for directional H4K16ac/RNA gene overlap."""

import argparse
import csv
from pathlib import Path

from scipy.stats import fisher_exact


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--loss", required=True, type=Path)
    parser.add_argument("--gain", type=Path, help="Optional H4K16ac gain CSV")
    parser.add_argument("--rna", required=True, type=Path)
    parser.add_argument("--background", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--chip-gene-column", default="SYMBOL")
    parser.add_argument("--rna-gene-column", default="gene_name")
    parser.add_argument("--rna-fold-change-column", default="log2FoldChange")
    parser.add_argument("--background-gene-column", default="Gene")
    return parser.parse_args()


def clean_gene(value):
    value = (value or "").strip().upper()
    return value if value and value not in {"NA", "."} else None


def read_gene_set(path, column):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if column not in (reader.fieldnames or []):
            raise ValueError(f"Column '{column}' not found in {path}")
        return {
            gene
            for row in reader
            if (gene := clean_gene(row[column])) is not None
        }


def read_rna_sets(path, gene_column, fold_change_column):
    down, up = set(), set()
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for column in (gene_column, fold_change_column):
            if column not in (reader.fieldnames or []):
                raise ValueError(f"Column '{column}' not found in {path}")
        for row in reader:
            gene = clean_gene(row[gene_column])
            if gene is None:
                continue
            try:
                fold_change = float(row[fold_change_column])
            except (TypeError, ValueError):
                continue
            if fold_change < 0:
                down.add(gene)
            elif fold_change > 0:
                up.add(gene)
    return down, up


def run_test(name, chip_genes, rna_genes, universe):
    # Fisher's test must contain only genes eligible for the comparison.
    chip_genes = chip_genes & universe
    rna_genes = rna_genes & universe

    a = len(chip_genes & rna_genes)          # ChIP change + RNA change
    b = len(chip_genes - rna_genes)          # ChIP change + not RNA change
    c = len(rna_genes - chip_genes)          # no ChIP change + RNA change
    d = len(universe - chip_genes - rna_genes)  # neither
    odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")

    print(f"\n{name}")
    print("Contingency table: [[overlap, ChIP only], [RNA only, neither]]")
    print(f"[[{a}, {b}], [{c}, {d}]]")
    print(f"H4K16ac-associated genes in universe: {len(chip_genes)}")
    print(f"Directional RNA genes in universe: {len(rna_genes)}")
    print(f"Overlapping genes: {a}")
    print(f"Odds ratio: {odds_ratio:.6g}")
    print(f"One-sided Fisher exact p-value: {p_value:.6g}")

    return {
        "test": name,
        "universe_genes": len(universe),
        "chip_genes": len(chip_genes),
        "rna_genes": len(rna_genes),
        "overlap": a,
        "chip_and_rna": a,
        "chip_and_not_rna": b,
        "no_chip_and_rna": c,
        "no_chip_and_not_rna": d,
        "odds_ratio": odds_ratio,
        "fisher_p_value_greater": p_value,
    }


def main():
    args = arguments()
    universe = read_gene_set(args.background, args.background_gene_column)
    loss_genes = read_gene_set(args.loss, args.chip_gene_column)
    rna_down, rna_up = read_rna_sets(
        args.rna, args.rna_gene_column, args.rna_fold_change_column
    )

    if not universe:
        raise ValueError("The background universe is empty")

    results = [run_test("H4K16ac loss vs RNA down", loss_genes, rna_down, universe)]

    if args.gain is not None:
        gain_genes = read_gene_set(args.gain, args.chip_gene_column)
        results.append(run_test("H4K16ac gain vs RNA up", gain_genes, rna_up, universe))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    print(f"\nResults saved to: {args.output}")


if __name__ == "__main__":
    main()
