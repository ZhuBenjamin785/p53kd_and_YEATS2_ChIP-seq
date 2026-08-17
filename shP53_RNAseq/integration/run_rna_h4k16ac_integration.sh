#!/usr/bin/env bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=rna_h4k16ac_integrate
#SBATCH -t 02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH -N 1
#SBATCH --output=shared/log/rna_h4k16ac_integration-%j.out
#SBATCH --error=shared/log/rna_h4k16ac_integration-%j.err

set -euo pipefail

PROJECT=/gpfs/projects/b1042/LauberthLab/BenFolder
RNA_FILE="${1:-$PROJECT/shared/rna_seq_dea/shp53_vs_shLacZ_0hr/significant_results.csv}"
CHIP_FILE="${2:-$PROJECT/p53kdH4K16ac/diffbind_results/split_genes/significant_annotated.csv}"
OUTPUT_DIR="${3:-$PROJECT/shared/rna_chip_integration/p53KD_H4K16ac_vs_RNAseq}"
RSCRIPT=/home/nqp9093/.conda/envs/chipseeker/bin/Rscript

mkdir -p "$PROJECT/shared/log" "$OUTPUT_DIR"

if [[ ! -x "$RSCRIPT" ]]; then
    echo "ERROR: Rscript not found: $RSCRIPT" >&2
    exit 1
fi

"$RSCRIPT" "$PROJECT/shP53_RNAseq/integration/integrate_rna_h4k16ac.R" \
    "$RNA_FILE" "$CHIP_FILE" "$OUTPUT_DIR"
