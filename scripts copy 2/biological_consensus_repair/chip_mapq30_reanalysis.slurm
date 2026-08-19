#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=h4k16ac_mapq30_repair
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=shared/biological_consensus_repaired/logs/chip_mapq30_%j.out
#SBATCH --error=shared/biological_consensus_repaired/logs/chip_mapq30_%j.err
set -euo pipefail
ROOT="/gpfs/projects/b1042/LauberthLab/BenFolder"
BASE="$ROOT/shared/biological_consensus_repaired/chipseq"
QC="$BASE/qc/species_cross_mapping.tsv"
HUMAN_IN="$ROOT/p53kdH4K16ac/p53kdbamfiles/human"
DM6_IN="$ROOT/p53kdH4K16ac/p53kdbamfiles/dm6"
TMPBASE="${SLURM_TMPDIR:-/tmp}/h4k16ac_mapq30_${SLURM_JOB_ID:-local}"
# Corrected BAMs are large, reproducible intermediates and stay in job-local
# scratch. Permanent outputs retain factors, duplicate metrics, peaks, tracks,
# DiffBind results, exact commands, and source-BAM checksums.
HUMAN_OUT="$TMPBASE/mapq30_deduplicated_bams/human"
DM6_OUT="$TMPBASE/mapq30_deduplicated_bams/dm6"
PEAKS="$BASE/mapq30_peaks"
DIFF="$BASE/diffbind_mapq30_broad_remove_duplicates"
FACTORS="$BASE/mapq30_spikein_scale_factors.tsv"
TRACK_TMP="$TMPBASE/tracks_mapq30"
TRACK_MANIFEST="$BASE/mapq30_track_validation_manifest.tsv"
mkdir -p "$HUMAN_OUT" "$DM6_OUT" "$PEAKS" "$DIFF" "$TRACK_TMP"
mkdir -p "$TMPBASE" "$BASE/duplicate_metrics_mapq30"
cd "$ROOT"
SAMTOOLS="/home/nqp9093/.conda/envs/rseqc_env/bin/samtools"
module load MACS3/3.0.2
module load deeptools/3.5.6
"$SAMTOOLS" --version

[[ -s "$QC" ]] || { echo "Missing cross-mapping QC: $QC" >&2; exit 1; }
rows=$(awk 'NR>1{n++} END{print n+0}' "$QC")
shared=$(awk 'NR>1{s+=$4} END{print s+0}' "$QC")
[[ "$rows" -eq 8 ]] || {
  echo "STOP: require all 8 libraries in the cross-mapping screen; rows=$rows" >&2
  exit 2
}
echo "MAPQ30 shared query names detected by the preflight screen: $shared; retained shared names will be removed from both species."

samples=(Scr_H4K16ac_1_S0_L001 Scr_H4K16ac_2_S0_L001 P53_H4K16ac_1_S0_L001 P53_H4K16ac_2_S0_L001 Scr_Input1_S0_L001 Scr_Input2_S0_L001 MutP53_Input1_S0_L001 MutP53_Input2_S0_L001)
for sample in "${samples[@]}"; do
  for species in human dm6; do
    in_dir="$HUMAN_IN"; out_dir="$HUMAN_OUT"; [[ "$species" == dm6 ]] && { in_dir="$DM6_IN"; out_dir="$DM6_OUT"; }
    # The source BAMs have zero duplicate flags. Reconstruct mate tags and
    # mark/remove coordinate duplicates before applying the MAPQ filter.
    # This is a sensitivity analysis: without UMIs, true histone signal can
    # be indistinguishable from PCR duplicates at identical coordinates.
    # Three piped processes run concurrently, so four helper threads each
    # stays within the 16-CPU allocation.
    "$SAMTOOLS" view -@ 4 -u -f 2 -F 2828 "$in_dir/${sample}.sorted.bam" |
      "$SAMTOOLS" sort -n -@ 4 -O bam - |
      "$SAMTOOLS" fixmate -@ 4 -m -O bam - - |
      "$SAMTOOLS" sort -@ 4 -o "$TMPBASE/${sample}.${species}.fixmate.sorted.bam" -
    "$SAMTOOLS" markdup -@ 8 -r -s "$TMPBASE/${sample}.${species}.fixmate.sorted.bam" \
      "$TMPBASE/${sample}.${species}.deduplicated.bam" \
      2> "$BASE/duplicate_metrics_mapq30/${sample}.${species}.txt"
    "$SAMTOOLS" view -@ 8 -bh -f 2 -F 2828 -q 30 \
      -o "$out_dir/${sample}.sorted.bam" "$TMPBASE/${sample}.${species}.deduplicated.bam"
    rm -f "$TMPBASE/${sample}.${species}.fixmate.sorted.bam" \
      "$TMPBASE/${sample}.${species}.deduplicated.bam"
  done
done

# Explicitly exclude any retained query name found in both independently
# aligned species. This is conservative post-alignment disambiguation, not a
# claim that separate alignment is equivalent to a competitive reference.
printf 'sample\tshared_mapq30_qnames_removed\n' > "$BASE/mapq30_removed_shared_qnames.tsv"
for sample in "${samples[@]}"; do
  "$SAMTOOLS" view -@ 4 -f 66 -F 2828 -q 30 "$HUMAN_OUT/${sample}.sorted.bam" |
    cut -f1 | LC_ALL=C sort -u > "$TMPBASE/${sample}.human.qnames"
  "$SAMTOOLS" view -@ 4 -f 66 -F 2828 -q 30 "$DM6_OUT/${sample}.sorted.bam" |
    cut -f1 | LC_ALL=C sort -u > "$TMPBASE/${sample}.dm6.qnames"
  comm -12 "$TMPBASE/${sample}.human.qnames" "$TMPBASE/${sample}.dm6.qnames" \
    > "$TMPBASE/${sample}.shared.qnames"
  nshared=$(wc -l < "$TMPBASE/${sample}.shared.qnames")
  printf '%s\t%s\n' "$sample" "$nshared" >> "$BASE/mapq30_removed_shared_qnames.tsv"
  for species in human dm6; do
    out_dir="$HUMAN_OUT"; [[ "$species" == dm6 ]] && out_dir="$DM6_OUT"
    bam="$out_dir/${sample}.sorted.bam"
    if [[ "$nshared" -gt 0 ]]; then
      "$SAMTOOLS" view -@ 8 -bh -N "^$TMPBASE/${sample}.shared.qnames" \
        -o "$TMPBASE/${sample}.${species}.unambiguous.bam" "$bam"
      mv "$TMPBASE/${sample}.${species}.unambiguous.bam" "$bam"
    fi
    "$SAMTOOLS" index -@ 8 "$bam"
    "$SAMTOOLS" quickcheck -v "$bam"
  done
done

printf 'sample\tdm6_fragments_mapq30\tscale_factor\n' > "$FACTORS.tmp"
counts=()
minimum=""
for sample in "${samples[@]}"; do
  n=$("$SAMTOOLS" view -@ 8 -c -f 66 -F 2828 -q 30 "$DM6_OUT/${sample}.sorted.bam")
  counts+=("$n"); [[ -z "$minimum" || "$n" -lt "$minimum" ]] && minimum="$n"
done
for i in "${!samples[@]}"; do
  sf=$(awk -v m="$minimum" -v n="${counts[$i]}" 'BEGIN{printf "%.12f",m/n}')
  printf '%s\t%s\t%s\n' "${samples[$i]}" "${counts[$i]}" "$sf" >> "$FACTORS.tmp"
done
mv "$FACTORS.tmp" "$FACTORS"

chip=(Scr_H4K16ac_1_S0_L001 Scr_H4K16ac_2_S0_L001 P53_H4K16ac_1_S0_L001 P53_H4K16ac_2_S0_L001)
input=(Scr_Input1_S0_L001 Scr_Input2_S0_L001 MutP53_Input1_S0_L001 MutP53_Input2_S0_L001)
peak=(Scramble_H4K16ac_rep1 Scramble_H4K16ac_rep2 p53KD_H4K16ac_rep1 p53KD_H4K16ac_rep2)
for i in "${!chip[@]}"; do
  mkdir -p "$PEAKS/${peak[$i]}"
  macs3 callpeak -t "$HUMAN_OUT/${chip[$i]}.sorted.bam" -c "$HUMAN_OUT/${input[$i]}.sorted.bam" \
    -f BAMPE -g hs -n "${peak[$i]}" --outdir "$PEAKS/${peak[$i]}" --broad --broad-cutoff 0.1 --keep-dup 1
done

RS="/home/nqp9093/.conda/envs/chipseeker/bin/Rscript"
export CHIP_RESULT_STATUS="PROVISIONAL — retained MAPQ30 cross-species query names excluded; separate-alignment and replicate-QC limitations remain"
"$RS" --vanilla "$ROOT/shared/scripts/biological_consensus_repair/diffbind_broad_sensitivity.R" \
 "$DIFF" TRUE "$HUMAN_OUT" "$DM6_OUT" "$PEAKS"

export HUMAN_VALIDATED_BAM_DIR="$HUMAN_OUT"
export VALIDATED_FACTOR_FILE="$FACTORS"
export CORRECTED_TRACK_DIR="$TRACK_TMP"
bash "$ROOT/shared/scripts/biological_consensus_repair/bamcompare_disambiguated_spikein.sh"
printf 'track\tbytes\tsha256\n' > "$TRACK_MANIFEST.tmp"
for bw in "$TRACK_TMP"/*.bw; do
  bytes=$(stat -c '%s' "$bw")
  digest=$(sha256sum "$bw" | awk '{print $1}')
  printf '%s\t%s\t%s\n' "$(basename "$bw")" "$bytes" "$digest" >> "$TRACK_MANIFEST.tmp"
done
mv "$TRACK_MANIFEST.tmp" "$TRACK_MANIFEST"
rm -r "$TMPBASE"
echo "MAPQ30 ChIP repair completed: $BASE"
