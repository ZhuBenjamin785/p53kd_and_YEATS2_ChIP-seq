#!/usr/bin/env python3
"""Build scope-matched universes and run all-peak and promoter-peak ORA."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
SINGLE_RUNNER = SCRIPT_DIR / "run_h4k16ac_rna_ora.py"
DEFAULT_BASE = Path("/gpfs/projects/b1042/LauberthLab/BenFolder")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Construct separate all-peak and promoter-peak gene universes from "
            "tested RNA/ChIP results, then run clusterProfiler ORA for both."
        )
    )
    parser.add_argument(
        "--rna-results", type=Path,
        default=DEFAULT_BASE / "shared/rna_seq_dea/shp53_vs_shLacZ_0hr/results.csv",
        help="Complete RNA differential-expression results, including nonsignificant genes.",
    )
    parser.add_argument(
        "--chip-annotations", type=Path,
        default=DEFAULT_BASE / "shared/chipseq_summary_plots/p53KD/diffbind_peak_annotations.csv",
        help="Complete annotated DiffBind results, including nonsignificant peaks.",
    )
    parser.add_argument(
        "--integration-root", type=Path,
        default=DEFAULT_BASE / "shared/rna_chip_integration/p53KD_H4K16ac_vs_RNAseq",
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=DEFAULT_BASE / "shared/rna_chip_integration/p53KD_H4K16ac_vs_RNAseq/ora_matched_universes",
    )
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--show", type=int, default=15)
    parser.add_argument("--min-gs-size", type=int, default=10)
    parser.add_argument("--max-gs-size", type=int, default=500)
    parser.add_argument(
        "--rscript", type=Path,
        default=Path("/home/nqp9093/.conda/envs/chipseeker/bin/Rscript"),
    )
    return parser.parse_args()


def finite(value: str | None) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def clean_symbol(value: str | None) -> str | None:
    if value is None:
        return None
    value = value.strip().upper()
    return value if value and value not in {"NA", "."} else None


def require_columns(reader: csv.DictReader, needed: set[str], path: Path) -> None:
    found = set(reader.fieldnames or [])
    missing = needed - found
    if missing:
        raise SystemExit(f"Missing columns in {path}: {', '.join(sorted(missing))}")


def read_rna_eligible(path: Path) -> set[str]:
    """Genes eligible for RNA significance testing in the integration."""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        require_columns(reader, {"gene_name", "log2FoldChange", "padj"}, path)
        return {
            gene
            for row in reader
            if (gene := clean_symbol(row["gene_name"])) is not None
            and finite(row["log2FoldChange"])
            and finite(row["padj"])
        }


def read_chip_eligible(path: Path) -> tuple[set[str], set[str]]:
    """Return genes linked to any tested peak and to any tested promoter peak."""
    all_peak_genes: set[str] = set()
    promoter_peak_genes: set[str] = set()
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        require_columns(reader, {"SYMBOL", "Fold", "FDR", "annotation"}, path)
        for row in reader:
            gene = clean_symbol(row["SYMBOL"])
            if gene is None or not finite(row["Fold"]) or not finite(row["FDR"]):
                continue
            all_peak_genes.add(gene)
            # This deliberately matches grepl("promoter", Peak_annotation,
            # ignore.case=TRUE) in the upstream integration script.
            if "promoter" in (row["annotation"] or "").lower():
                promoter_peak_genes.add(gene)
    return all_peak_genes, promoter_peak_genes


def write_universe(path: Path, genes: set[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Gene"])
        writer.writerows((gene,) for gene in sorted(genes))


def write_audit(path: Path, rows: list[dict[str, str | int]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def foreground_paths(integration_root: Path, scope: str) -> dict[str, Path]:
    table_dir = integration_root / scope / "tables"
    return {
        "loss-down": table_dir / "loss_down.csv",
        "gain-up": table_dir / "gain_up.csv",
        "loss-up": table_dir / "loss_up.csv",
        "gain-down": table_dir / "gain_down.csv",
    }


def run_scope(args: argparse.Namespace, scope: str, universe: Path) -> int:
    inputs = foreground_paths(args.integration_root, scope)
    command = [
        "python", str(SINGLE_RUNNER),
        "--background", str(universe),
        "--loss-down", str(inputs["loss-down"]),
        "--gain-up", str(inputs["gain-up"]),
        "--loss-up", str(inputs["loss-up"]),
        "--gain-down", str(inputs["gain-down"]),
        "--outdir", str(args.outdir / scope),
        "--fdr", str(args.fdr),
        "--show", str(args.show),
        "--min-gs-size", str(args.min_gs_size),
        "--max-gs-size", str(args.max_gs_size),
        "--rscript", str(args.rscript),
    ]
    print(f"\nRunning {scope} ORA with {universe}...", flush=True)
    return subprocess.run(command, check=False).returncode


def main() -> int:
    args = parse_args()
    for path in (args.rna_results, args.chip_annotations, args.integration_root):
        if not path.exists():
            raise SystemExit(f"Required input does not exist: {path}")

    rna_genes = read_rna_eligible(args.rna_results)
    chip_all, chip_promoter = read_chip_eligible(args.chip_annotations)
    universes = {
        "all_peaks": rna_genes & chip_all,
        "promoter_peaks": rna_genes & chip_promoter,
    }
    universe_dir = args.outdir / "universes"
    universe_paths = {
        scope: universe_dir / f"{scope}_eligible_gene_universe.csv"
        for scope in universes
    }
    for scope, genes in universes.items():
        if not genes:
            raise SystemExit(f"Derived {scope} universe is empty; refusing to run ORA.")
        write_universe(universe_paths[scope], genes)

    args.outdir.mkdir(parents=True, exist_ok=True)
    audit_rows = [
        {
            "Scope": "all_peaks",
            "RNA_eligible_genes": len(rna_genes),
            "ChIP_eligible_genes": len(chip_all),
            "Matched_universe_genes": len(universes["all_peaks"]),
            "Promoter_rule": "not applicable",
        },
        {
            "Scope": "promoter_peaks",
            "RNA_eligible_genes": len(rna_genes),
            "ChIP_eligible_genes": len(chip_promoter),
            "Matched_universe_genes": len(universes["promoter_peaks"]),
            "Promoter_rule": "annotation contains promoter (case-insensitive)",
        },
    ]
    write_audit(universe_dir / "matched_universe_audit.csv", audit_rows)
    print("Matched universes:", flush=True)
    for row in audit_rows:
        print(
            f"  {row['Scope']}: {row['Matched_universe_genes']:,} genes "
            f"(RNA eligible={row['RNA_eligible_genes']:,}; "
            f"ChIP eligible={row['ChIP_eligible_genes']:,})",
            flush=True,
        )

    return_codes = [
        run_scope(args, scope, universe_paths[scope])
        for scope in ("all_peaks", "promoter_peaks")
    ]
    return max(return_codes)


if __name__ == "__main__":
    raise SystemExit(main())
