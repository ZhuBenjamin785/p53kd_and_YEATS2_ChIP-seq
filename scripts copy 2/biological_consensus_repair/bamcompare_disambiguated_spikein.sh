#!/bin/bash
# Generate ChIP/Input log2 tracks without per-library BPM renormalization.
# This script deliberately requires explicitly validated per-library scale
# factors; it will not silently reuse the unfiltered legacy counts. Validation
# may be from competitive alignment or a clearly labelled MAPQ sensitivity.
set -euo pipefail
ROOT="/gpfs/projects/b1042/LauberthLab/BenFolder"
BAMS="${HUMAN_VALIDATED_BAM_DIR:-${HUMAN_DISAMBIGUATED_BAM_DIR:-$ROOT/shared/biological_consensus_repaired/chipseq/disambiguated_bams/human}}"
FACTORS="${VALIDATED_FACTOR_FILE:-${DISAMBIGUATED_FACTOR_FILE:-$ROOT/shared/biological_consensus_repaired/chipseq/disambiguated_spikein_scale_factors.tsv}}"
OUT="${CORRECTED_TRACK_DIR:-$ROOT/shared/biological_consensus_repaired/chipseq/tracks}"
DEEPTOOLS="${DEEPTOOLS_BIN:-/software/deeptools/3.5.6/bin}"
[[ -s "$FACTORS" ]] || { echo "UNRESOLVED: missing validated factor file: $FACTORS" >&2; exit 2; }
mkdir -p "$OUT"

factor() { awk -F '\t' -v s="$1" 'NR==1{for(i=1;i<=NF;i++){if($i=="sample")a=i;if($i=="scale_factor")b=i}} NR>1&&$a==s{print $b}' "$FACTORS"; }
compare() {
  local chip="$1" input="$2" output="$3" sf_chip sf_input
  sf_chip=$(factor "$chip"); sf_input=$(factor "$input")
  [[ -n "$sf_chip" && -n "$sf_input" ]] || { echo "Missing scale factor for $chip or $input" >&2; exit 1; }
  [[ -s "$BAMS/$chip.sorted.bam" && -s "$BAMS/$input.sorted.bam" ]] || { echo "Missing validated BAM" >&2; exit 1; }
  "$DEEPTOOLS/bamCompare" -b1 "$BAMS/$chip.sorted.bam" -b2 "$BAMS/$input.sorted.bam" \
    --operation log2 --scaleFactors "$sf_chip:$sf_input" --normalizeUsing None \
    --binSize 25 --numberOfProcessors "${SLURM_CPUS_PER_TASK:-4}" --outFileFormat bigwig -o "$OUT/$output"
}
compare P53_H4K16ac_1_S0_L001 MutP53_Input1_S0_L001 P53KD_rep1_ChIP_vs_Input_log2.spikein.bw
compare P53_H4K16ac_2_S0_L001 MutP53_Input2_S0_L001 P53KD_rep2_ChIP_vs_Input_log2.spikein.bw
compare Scr_H4K16ac_1_S0_L001 Scr_Input1_S0_L001 Scramble_rep1_ChIP_vs_Input_log2.spikein.bw
compare Scr_H4K16ac_2_S0_L001 Scr_Input2_S0_L001 Scramble_rep2_ChIP_vs_Input_log2.spikein.bw
