#!/usr/bin/env python3
"""Python entry point for the clusterProfiler H4K16ac/RNA ORA workflow."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


DEFAULT_RSCRIPT = Path("/home/nqp9093/.conda/envs/chipseeker/bin/Rscript")
R_WORKFLOW = Path(__file__).with_name("h4k16ac_rna_ora.R")
GENE_COLUMNS = {"gene", "symbol", "gene_name", "gene_symbol"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run GO BP, KEGG, and Reactome hypergeometric ORA using "
            "clusterProfiler and a user-supplied gene universe."
        )
    )
    parser.add_argument("--background", required=True, type=Path)
    parser.add_argument("--loss-down", required=True, type=Path)
    parser.add_argument("--gain-up", required=True, type=Path)
    parser.add_argument("--loss-up", required=True, type=Path)
    parser.add_argument("--gain-down", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--show", type=int, default=15)
    parser.add_argument("--min-gs-size", type=int, default=10)
    parser.add_argument("--max-gs-size", type=int, default=500)
    parser.add_argument(
        "--rscript", type=Path, default=DEFAULT_RSCRIPT,
        help=f"Rscript executable (default: {DEFAULT_RSCRIPT})",
    )
    return parser.parse_args()


def inspect_gene_table(path: Path) -> tuple[str, int]:
    """Return the recognized symbol column and the number of unique symbols."""
    if not path.is_file():
        raise ValueError(f"Input does not exist or is not a file: {path}")

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        first_line = handle.readline()
        handle.seek(0)
        delimiter = "\t" if "\t" in first_line else ","
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Input has no header: {path}")
        by_lower = {name.strip().lower(): name for name in reader.fieldnames}
        matches = [by_lower[name] for name in GENE_COLUMNS if name in by_lower]
        if not matches:
            raise ValueError(
                f"No recognized gene-symbol column in {path}; found: "
                + ", ".join(reader.fieldnames)
            )
        column = matches[0]
        genes = {
            row[column].strip()
            for row in reader
            if row.get(column) is not None and row[column].strip()
        }
    return column, len(genes)


def main() -> int:
    args = parse_args()
    if not 0 < args.fdr <= 1:
        raise SystemExit("--fdr must be in (0, 1].")
    if args.show < 1 or args.min_gs_size < 1 or args.max_gs_size < args.min_gs_size:
        raise SystemExit("Plot count and gene-set size limits must be positive and ordered.")
    if not args.rscript.is_file():
        raise SystemExit(f"Rscript executable not found: {args.rscript}")
    if not R_WORKFLOW.is_file():
        raise SystemExit(f"R workflow not found: {R_WORKFLOW}")

    inputs = {
        "background": args.background,
        "loss-down": args.loss_down,
        "gain-up": args.gain_up,
        "loss-up": args.loss_up,
        "gain-down": args.gain_down,
    }
    print("Validated inputs:", flush=True)
    for label, path in inputs.items():
        try:
            column, count = inspect_gene_table(path)
        except ValueError as error:
            raise SystemExit(str(error)) from error
        print(f"  {label}: {count:,} unique symbols ({column}; {path})", flush=True)

    command = [
        str(args.rscript), str(R_WORKFLOW),
        "--background", str(args.background),
        "--loss-down", str(args.loss_down),
        "--gain-up", str(args.gain_up),
        "--loss-up", str(args.loss_up),
        "--gain-down", str(args.gain_down),
        "--outdir", str(args.outdir),
        "--fdr", str(args.fdr),
        "--show", str(args.show),
        "--min-gs-size", str(args.min_gs_size),
        "--max-gs-size", str(args.max_gs_size),
    ]
    print("\nStarting clusterProfiler ORA...", flush=True)
    try:
        completed = subprocess.run(command, check=False)
    except OSError as error:
        raise SystemExit(f"Could not start R workflow: {error}") from error
    if completed.returncode != 0:
        print(f"ORA failed with exit status {completed.returncode}.", file=sys.stderr)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
