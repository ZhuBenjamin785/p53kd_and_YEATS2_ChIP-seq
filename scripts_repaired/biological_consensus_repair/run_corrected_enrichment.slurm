#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=consensus_enrichment_repair
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=shared/biological_consensus_repaired/logs/enrichment_%j.out
#SBATCH --error=shared/biological_consensus_repaired/logs/enrichment_%j.err

set -euo pipefail
ROOT="/gpfs/projects/b1042/LauberthLab/BenFolder"
PY="/home/nqp9093/.conda/envs/pydeseq2/bin/python"
RS="/home/nqp9093/.conda/envs/chipseeker/bin/Rscript"
ORA="$ROOT/shared/scripts/over_representation_analysis/run_h4k16ac_rna_ora.py"
BASE="$ROOT/shared/biological_consensus_repaired"
mkdir -p "$BASE/logs" "$BASE/ora" "$BASE/fea" "$BASE/gsea"
cd "$ROOT"

run_ora_scope() {
  local scope="$1"
  "$PY" "$ORA" \
    --background "$BASE/integration/gene_sets/${scope}__eligible_universe.csv" \
    --loss-down "$BASE/integration/gene_sets/${scope}__loss_down.csv" \
    --gain-up "$BASE/integration/gene_sets/${scope}__gain_up.csv" \
    --loss-up "$BASE/integration/gene_sets/${scope}__loss_up.csv" \
    --gain-down "$BASE/integration/gene_sets/${scope}__gain_down.csv" \
    --outdir "$BASE/ora/$scope" --rscript "$RS" --fdr 0.05 \
    --min-gs-size 10 --max-gs-size 500
}

run_ora_scope all_peaks
run_ora_scope promoter_peaks

"$RS" "$ROOT/shared/scripts/biological_consensus_repair/directional_rna_fea.R" \
  --background "$BASE/fea/rna_eligible_universe.csv" \
  --up "$BASE/fea/rna_up.csv" --down "$BASE/fea/rna_down.csv" \
  --outdir "$BASE/fea/directional" --fdr 0.05 --min-size 10 --max-size 500

"$RS" "$ROOT/shared/scripts/biological_consensus_repair/gsea_all_eligible.R" \
  --results "$ROOT/shared/rna_seq_dea/shp53_vs_shLacZ_0hr/results.csv" \
  --outdir "$BASE/gsea/all_eligible" --fdr 0.05 --min-size 15 --max-size 500 --seed 20260818

MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/matplotlib_consensus_repair" \
  "$PY" "$ROOT/shared/scripts/biological_consensus_repair/plot_corrected_enrichment.py"
