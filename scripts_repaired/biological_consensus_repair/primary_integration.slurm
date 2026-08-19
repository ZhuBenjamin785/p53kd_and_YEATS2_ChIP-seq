#!/bin/bash
# Rebuild the required 2,100-gene primary integration before enrichment.
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=consensus_primary_integration
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=shared/biological_consensus_repaired/logs/primary_integration_%j.out
#SBATCH --error=shared/biological_consensus_repaired/logs/primary_integration_%j.err
set -euo pipefail
ROOT="/gpfs/projects/b1042/LauberthLab/BenFolder"
PY="/home/nqp9093/.conda/envs/pydeseq2/bin/python"
mkdir -p "$ROOT/shared/biological_consensus_repaired/logs"
cd "$ROOT"
"$PY" shared/scripts/biological_consensus_repair/repair_integration.py \
  --rna shared/rna_seq_dea/shp53_vs_shLacZ_0hr/results.csv \
  --chip shared/chipseq_summary_plots/p53KD/diffbind_peak_annotations.csv \
  --outdir shared/biological_consensus_repaired/integration \
  --expected-all-universe-size 2100
