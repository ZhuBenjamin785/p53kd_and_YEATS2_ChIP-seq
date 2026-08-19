#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=h4k16ac_qc_repair
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=shared/biological_consensus_repaired/logs/chip_qc_%j.out
#SBATCH --error=shared/biological_consensus_repaired/logs/chip_qc_%j.err

set -euo pipefail
ROOT="/gpfs/projects/b1042/LauberthLab/BenFolder"
OUT="$ROOT/shared/biological_consensus_repaired/chipseq/qc"
TMPBASE="${SLURM_TMPDIR:-/tmp}/h4k16ac_qc_${SLURM_JOB_ID:-local}"
mkdir -p "$OUT/flagstat" "$TMPBASE"
cd "$ROOT"
module load samtools
command -v samtools >/dev/null || { echo "samtools unavailable" >&2; exit 1; }

samples=(Scr_H4K16ac_1_S0_L001 Scr_H4K16ac_2_S0_L001 P53_H4K16ac_1_S0_L001 P53_H4K16ac_2_S0_L001 Scr_Input1_S0_L001 Scr_Input2_S0_L001 MutP53_Input1_S0_L001 MutP53_Input2_S0_L001)
peak_names=(Scramble_H4K16ac_rep1 Scramble_H4K16ac_rep2 p53KD_H4K16ac_rep1 p53KD_H4K16ac_rep2)

printf 'sample\thuman_primary_proper_fragments_mapq30\tdm6_primary_proper_fragments_mapq30\tshared_qnames_mapq30\tshared_fraction_of_smaller_species_set\n' > "$OUT/species_cross_mapping.tsv"
printf 'sample\ttotal_primary_proper_fragments_mapq30\tdistinct_fragment_coordinates\tsingleton_coordinates\tdoubleton_coordinates\tNRF\tPBC1\tPBC2\tFRiP_read1_fraction\n' > "$OUT/library_complexity_frip.tsv"

for i in "${!samples[@]}"; do
  sample="${samples[$i]}"
  human="$ROOT/p53kdH4K16ac/p53kdbamfiles/human/${sample}.sorted.bam"
  fly="$ROOT/p53kdH4K16ac/p53kdbamfiles/dm6/${sample}.sorted.bam"
  [[ -s "$human" && -s "$human.bai" && -s "$fly" && -s "$fly.bai" ]] || { echo "Missing BAM/index for $sample" >&2; exit 1; }
  samtools quickcheck -v "$human" "$fly"
  samtools flagstat -@ 4 "$human" > "$OUT/flagstat/${sample}.human.flagstat.txt"
  samtools flagstat -@ 4 "$fly" > "$OUT/flagstat/${sample}.dm6.flagstat.txt"

  # One primary, properly paired, MAPQ>=30 first mate represents each fragment.
  samtools view -@ 4 -f 66 -F 2816 -q 30 "$human" | cut -f1 | sort -u > "$TMPBASE/${sample}.human.qnames"
  samtools view -@ 4 -f 66 -F 2816 -q 30 "$fly" | cut -f1 | sort -u > "$TMPBASE/${sample}.dm6.qnames"
  nh=$(wc -l < "$TMPBASE/${sample}.human.qnames")
  nf=$(wc -l < "$TMPBASE/${sample}.dm6.qnames")
  no=$(comm -12 "$TMPBASE/${sample}.human.qnames" "$TMPBASE/${sample}.dm6.qnames" | wc -l)
  shared_frac=$(awk -v o="$no" -v h="$nh" -v f="$nf" 'BEGIN{m=(h<f?h:f); if(m) printf "%.6f",o/m; else print "NA"}')
  printf '%s\t%s\t%s\t%s\t%s\n' "$sample" "$nh" "$nf" "$no" "$shared_frac" >> "$OUT/species_cross_mapping.tsv"

  if (( i < 4 )); then
    coords="$TMPBASE/${sample}.fragment_coordinates"
    samtools view -@ 4 -f 66 -F 2816 -q 30 "$human" |
      awk 'BEGIN{OFS="\t"} $7=="=" {a=$4; b=$8; if(a>b){t=a;a=b;b=t}; print $3,a,b}' |
      sort | uniq -c > "$coords"
    total=$(awk '{s+=$1} END{print s+0}' "$coords")
    distinct=$(wc -l < "$coords")
    n1=$(awk '$1==1{n++} END{print n+0}' "$coords")
    n2=$(awk '$1==2{n++} END{print n+0}' "$coords")
    nrf=$(awk -v d="$distinct" -v t="$total" 'BEGIN{if(t)printf "%.6f",d/t;else print "NA"}')
    pbc1=$(awk -v n="$n1" -v d="$distinct" 'BEGIN{if(d)printf "%.6f",n/d;else print "NA"}')
    pbc2=$(awk -v n="$n1" -v d="$n2" 'BEGIN{if(d)printf "%.6f",n/d;else print "Inf"}')
    peak="$ROOT/p53kdH4K16ac/macs3_results_p53kd/peaks/${peak_names[$i]}/${peak_names[$i]}_peaks.broadPeak"
    in_peak=$(samtools view -@ 4 -c -f 66 -F 2816 -q 30 -L "$peak" "$human")
    frip=$(awk -v n="$in_peak" -v t="$total" 'BEGIN{if(t)printf "%.6f",n/t;else print "NA"}')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$total" "$distinct" "$n1" "$n2" "$nrf" "$pbc1" "$pbc2" "$frip" >> "$OUT/library_complexity_frip.tsv"
  fi
done

rm -r "$TMPBASE"
echo "ChIP QC completed: $OUT"

