#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 02:00:00
#SBATCH --mem=16G
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --job-name=split_p53kd_genes
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

project_dir="${PROJECT_DIR:-/gpfs/projects/b1042/LauberthLab/BenFolder}"
input_csv="${DIFFBIND_INPUT:-${project_dir}/p53kdH4K16ac/diffbind_results/DiffBind_all_peaks.csv}"
output_dir="${DIFFBIND_OUTPUT:-${project_dir}/p53kdH4K16ac/diffbind_results/split_genes}"
fdr_cutoff="${FDR_CUTOFF:-0.05}"
logfc_cutoff="${LOGFC_CUTOFF:-0.05}"
rscript="${RSCRIPT:-/gpfs/home/nqp9093/.conda/envs/chipseeker/bin/Rscript}"

[[ -s "$input_csv" ]] || { echo "Missing or empty DiffBind table: $input_csv" >&2; exit 1; }
[[ -x "$rscript" ]] || { echo "Rscript not found or not executable: $rscript" >&2; exit 1; }

mkdir -p "$output_dir"
"$rscript" --vanilla \
  "$project_dir/shared/scripts/overlap_analysis/split_diffbind_genes.R" \
  "$input_csv" "$output_dir" "$fdr_cutoff" "$logfc_cutoff"
