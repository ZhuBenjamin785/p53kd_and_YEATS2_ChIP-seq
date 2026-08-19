#!/bin/bash
# Paired-end insert-size QC for the four H4K16ac libraries. Existing BAMs are
# read only; output is a compact samtools-stats table.
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=h4k16ac_insert_qc
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=shared/biological_consensus_repaired/logs/chip_insert_qc_%j.out
#SBATCH --error=shared/biological_consensus_repaired/logs/chip_insert_qc_%j.err
set -euo pipefail
ROOT="/gpfs/projects/b1042/LauberthLab/BenFolder"
OUT="$ROOT/shared/biological_consensus_repaired/chipseq/qc/insert_size"
mkdir -p "$OUT"
cd "$ROOT"
module load samtools
samples=(Scr_H4K16ac_1_S0_L001 Scr_H4K16ac_2_S0_L001 P53_H4K16ac_1_S0_L001 P53_H4K16ac_2_S0_L001)
printf 'sample\taverage_insert_size\tinsert_size_standard_deviation\tinward_oriented_pairs\toutward_oriented_pairs\tother_oriented_pairs\tmapq30_fragments\tfragments_gt_1000bp\tfraction_gt_1000bp\n' > "$OUT/insert_size_summary.tsv"
for sample in "${samples[@]}"; do
  stats="$OUT/${sample}.samtools_stats.txt"
  samtools stats -@ 4 -f 2 -F 2816 \
    "$ROOT/p53kdH4K16ac/p53kdbamfiles/human/${sample}.sorted.bam" > "$stats"
  avg=$(awk -F '\t' '$1=="SN" && $2=="insert size average:"{print $3}' "$stats")
  sd=$(awk -F '\t' '$1=="SN" && $2=="insert size standard deviation:"{print $3}' "$stats")
  inward=$(awk -F '\t' '$1=="SN" && $2=="inward oriented pairs:"{print $3}' "$stats")
  outward=$(awk -F '\t' '$1=="SN" && $2=="outward oriented pairs:"{print $3}' "$stats")
  other=$(awk -F '\t' '$1=="SN" && $2=="pairs with other orientation:"{print $3}' "$stats")
  read -r total long fraction < <(
    samtools view -@ 4 -f 66 -F 2816 -q 30 \
      "$ROOT/p53kdH4K16ac/p53kdbamfiles/human/${sample}.sorted.bam" |
      awk '{n++; t=$9; if(t<0)t=-t; if(t>1000)l++} END{printf "%d %d %.8f\n",n+0,l+0,n?l/n:0}'
  )
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$avg" "$sd" "$inward" "$outward" "$other" "$total" "$long" "$fraction" >> "$OUT/insert_size_summary.tsv"
done
echo "Wrote $OUT/insert_size_summary.tsv"
