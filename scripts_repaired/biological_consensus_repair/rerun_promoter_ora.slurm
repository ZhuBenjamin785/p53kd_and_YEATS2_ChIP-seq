#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=promoter_ora_verified
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=shared/biological_consensus_repaired/logs/promoter_ora_%j.out
#SBATCH --error=shared/biological_consensus_repaired/logs/promoter_ora_%j.err
set -euo pipefail
ROOT="/gpfs/projects/b1042/LauberthLab/BenFolder"
BASE="$ROOT/shared/biological_consensus_repaired"
PY="/home/nqp9093/.conda/envs/pydeseq2/bin/python"
RS="/home/nqp9093/.conda/envs/chipseeker/bin/Rscript"
cd "$ROOT"
"$PY" "$ROOT/shared/scripts/over_representation_analysis/run_h4k16ac_rna_ora.py" \
 --background "$BASE/integration/gene_sets/promoter_peaks__eligible_universe.csv" \
 --loss-down "$BASE/integration/gene_sets/promoter_peaks__loss_down.csv" \
 --gain-up "$BASE/integration/gene_sets/promoter_peaks__gain_up.csv" \
 --loss-up "$BASE/integration/gene_sets/promoter_peaks__loss_up.csv" \
 --gain-down "$BASE/integration/gene_sets/promoter_peaks__gain_down.csv" \
 --outdir "$BASE/ora/promoter_peaks" --rscript "$RS" --fdr 0.05 \
 --min-gs-size 10 --max-gs-size 500
