#!/bin/bash
# Record the exact immutable BAM/peak inputs used by the repaired ChIP checks.
# This reads large files but does not modify them.
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=h4k16ac_input_hashes
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=shared/biological_consensus_repaired/logs/chip_hashes_%j.out
#SBATCH --error=shared/biological_consensus_repaired/logs/chip_hashes_%j.err
set -euo pipefail
ROOT="/gpfs/projects/b1042/LauberthLab/BenFolder"
OUT="$ROOT/shared/biological_consensus_repaired/provenance/CHIP_INPUT_SHA256SUMS.txt"
TMP="${OUT}.tmp"
cd "$ROOT"
samples=(Scr_H4K16ac_1_S0_L001 Scr_H4K16ac_2_S0_L001 P53_H4K16ac_1_S0_L001 P53_H4K16ac_2_S0_L001 Scr_Input1_S0_L001 Scr_Input2_S0_L001 MutP53_Input1_S0_L001 MutP53_Input2_S0_L001)
: > "$TMP"
for species in human dm6; do
  for sample in "${samples[@]}"; do
    sha256sum "p53kdH4K16ac/p53kdbamfiles/$species/${sample}.sorted.bam" >> "$TMP"
  done
done
for peak in \
  Scramble_H4K16ac_rep1/Scramble_H4K16ac_rep1_peaks.broadPeak \
  Scramble_H4K16ac_rep2/Scramble_H4K16ac_rep2_peaks.broadPeak \
  p53KD_H4K16ac_rep1/p53KD_H4K16ac_rep1_peaks.broadPeak \
  p53KD_H4K16ac_rep2/p53KD_H4K16ac_rep2_peaks.broadPeak; do
  sha256sum "p53kdH4K16ac/macs3_results_p53kd/peaks/$peak" >> "$TMP"
done
mv "$TMP" "$OUT"
echo "Wrote $OUT"
