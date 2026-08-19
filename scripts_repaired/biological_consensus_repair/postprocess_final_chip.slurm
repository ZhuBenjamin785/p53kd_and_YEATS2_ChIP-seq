#!/bin/bash
# Compact downstream processing after the corrected ChIP sensitivity finishes.
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=h4k16ac_postprocess
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=shared/biological_consensus_repaired/logs/chip_postprocess_%j.out
#SBATCH --error=shared/biological_consensus_repaired/logs/chip_postprocess_%j.err
set -euo pipefail
ROOT="/gpfs/projects/b1042/LauberthLab/BenFolder"
BASE="$ROOT/shared/biological_consensus_repaired"
PY="/home/nqp9093/.conda/envs/pydeseq2/bin/python"
RS="/home/nqp9093/.conda/envs/chipseeker/bin/Rscript"
cd "$ROOT"

FINAL="$BASE/chipseq/diffbind_mapq30_broad_remove_duplicates"
[[ -s "$FINAL/DiffBind_all_peaks.csv" && -s "$FINAL/diffbind_peak_annotations.csv" && -s "$FINAL/DiffBind_object.rds" ]] || {
  echo "Missing final ChIP outputs" >&2; exit 1;
}

# Recreate the broad-region (legacy BAM) sensitivity integration used to
# distinguish peak-window effects from the final coordinate-deduplicated run.
BROAD="$BASE/chipseq/diffbind_broad_remove_duplicates"
"$RS" --vanilla shared/scripts/biological_consensus_repair/extract_diffbind_replicate_qc.R \
  "$BROAD/DiffBind_object.rds" "$BROAD/replicate_correlations.csv"
"$PY" shared/scripts/biological_consensus_repair/repair_integration.py \
  --chip "$BROAD/diffbind_peak_annotations.csv" \
  --outdir "$BASE/integration_sensitivity_broad_remove_duplicates" \
  --expected-all-universe-size 0

"$RS" --vanilla shared/scripts/biological_consensus_repair/extract_diffbind_replicate_qc.R \
  "$FINAL/DiffBind_object.rds" "$FINAL/replicate_correlations.csv"

"$PY" shared/scripts/biological_consensus_repair/repair_integration.py \
  --chip "$FINAL/diffbind_peak_annotations.csv" \
  --outdir "$BASE/integration_sensitivity_mapq30_deduplicated" \
  --expected-all-universe-size 0

CONS="$BASE/chipseq/consensus_mapq30_deduplicated"
PEAKS="$BASE/chipseq/mapq30_peaks"
mkdir -p "$CONS"
for fraction in 0.25 0.50 0.75; do
  label=${fraction/./p}
  "$PY" shared/scripts/biological_consensus_repair/build_reciprocal_consensus.py \
    --rep1 "$PEAKS/Scramble_H4K16ac_rep1/Scramble_H4K16ac_rep1_peaks.broadPeak" \
    --rep2 "$PEAKS/Scramble_H4K16ac_rep2/Scramble_H4K16ac_rep2_peaks.broadPeak" \
    --reciprocal-fraction "$fraction" --output "$CONS/Scramble_reciprocal_${label}.bed"
  "$PY" shared/scripts/biological_consensus_repair/build_reciprocal_consensus.py \
    --rep1 "$PEAKS/p53KD_H4K16ac_rep1/p53KD_H4K16ac_rep1_peaks.broadPeak" \
    --rep2 "$PEAKS/p53KD_H4K16ac_rep2/p53KD_H4K16ac_rep2_peaks.broadPeak" \
    --reciprocal-fraction "$fraction" --output "$CONS/p53KD_reciprocal_${label}.bed"
done

"$PY" shared/scripts/biological_consensus_repair/summarize_chip_qc.py
"$PY" shared/scripts/biological_consensus_repair/compare_old_new.py \
  > "$BASE/summary/old_vs_new_numeric_comparison.txt"
MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/matplotlib_simple_plots" \
  "$PY" "$BASE/plotscripts/make_all_plots.py"
"$PY" shared/scripts/biological_consensus_repair/validate_repaired_outputs.py \
  > "$BASE/provenance/VALIDATION_RESULTS.txt"

# Compact, deterministic checksums. Large source BAM hashes are produced by
# checksum_chip_inputs.slurm and are not repeated here.
{
  sha256sum shared/rna_seq_dea/shp53_vs_shLacZ_0hr/results.csv
  sha256sum shared/chipseq_summary_plots/p53KD/diffbind_peak_annotations.csv
  find shared/scripts/biological_consensus_repair -maxdepth 1 -type f -print0 | sort -z | xargs -0 sha256sum
  sha256sum shared/scripts/over_representation_analysis/h4k16ac_rna_ora.R
} > "$BASE/provenance/SHA256SUMS.txt.tmp"
mv "$BASE/provenance/SHA256SUMS.txt.tmp" "$BASE/provenance/SHA256SUMS.txt"
find "$BASE" -type f ! -path '*/logs/*' ! -name '*SHA256SUMS*' -size -500M -print0 |
  sort -z | xargs -0 sha256sum > "$BASE/provenance/FINAL_OUTPUT_SHA256SUMS.txt.tmp"
mv "$BASE/provenance/FINAL_OUTPUT_SHA256SUMS.txt.tmp" "$BASE/provenance/FINAL_OUTPUT_SHA256SUMS.txt"
echo "Final ChIP postprocessing completed"
