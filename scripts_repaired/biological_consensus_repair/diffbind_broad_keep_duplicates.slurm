#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=diffbind_broad_keepdup
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=shared/biological_consensus_repaired/logs/diffbind_broad_keepdup_%j.out
#SBATCH --error=shared/biological_consensus_repaired/logs/diffbind_broad_keepdup_%j.err
set -euo pipefail
ROOT="/gpfs/projects/b1042/LauberthLab/BenFolder"
RS="/home/nqp9093/.conda/envs/chipseeker/bin/Rscript"
cd "$ROOT"
"$RS" --vanilla "$ROOT/shared/scripts/biological_consensus_repair/diffbind_broad_sensitivity.R" \
 "$ROOT/shared/biological_consensus_repaired/chipseq/diffbind_broad_keep_duplicates" FALSE
