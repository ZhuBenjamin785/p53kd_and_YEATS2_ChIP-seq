#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=h4k16ac_complexity_10pct
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=shared/biological_consensus_repaired/logs/chip_complexity_10pct_%j.out
#SBATCH --error=shared/biological_consensus_repaired/logs/chip_complexity_10pct_%j.err
set -euo pipefail
ROOT="/gpfs/projects/b1042/LauberthLab/BenFolder"
OUT="$ROOT/shared/biological_consensus_repaired/chipseq/qc"
TMPBASE="${SLURM_TMPDIR:-/tmp}/h4k16ac_complexity_${SLURM_JOB_ID:-local}"
mkdir -p "$OUT" "$TMPBASE"
module load samtools
samples=(Scr_H4K16ac_1_S0_L001 Scr_H4K16ac_2_S0_L001 P53_H4K16ac_1_S0_L001 P53_H4K16ac_2_S0_L001)

one_sample() {
 local sample="$1" bam="$ROOT/p53kdH4K16ac/p53kdbamfiles/human/${sample}.sorted.bam"
 local coords="$TMPBASE/${sample}.coords" result="$TMPBASE/${sample}.result"
 # Deterministic 10% subsample after common primary/proper-pair/MAPQ filters.
 samtools view -@ 2 -s 20260818.10 -f 66 -F 2816 -q 30 "$bam" |
   awk 'BEGIN{OFS="\t"} $7=="=" {a=$4;b=$8;if(a>b){t=a;a=b;b=t};print $3,a,b}' |
   sort --parallel=2 -S 4G | uniq -c > "$coords"
 total=$(awk '{s+=$1}END{print s+0}' "$coords"); distinct=$(wc -l < "$coords")
 n1=$(awk '$1==1{n++}END{print n+0}' "$coords"); n2=$(awk '$1==2{n++}END{print n+0}' "$coords")
 nrf=$(awk -v d="$distinct" -v t="$total" 'BEGIN{printf "%.6f",d/t}')
 pbc1=$(awk -v n="$n1" -v d="$distinct" 'BEGIN{printf "%.6f",n/d}')
 pbc2=$(awk -v n="$n1" -v d="$n2" 'BEGIN{if(d)printf "%.6f",n/d;else print "Inf"}')
 printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$total" "$distinct" "$n1" "$n2" "$nrf" "$pbc1" "$pbc2" > "$result"
}
export -f one_sample
for sample in "${samples[@]}"; do one_sample "$sample" & done
wait
printf 'sample\tsubsampled_fragments\tdistinct_coordinates\tsingletons\tdoubletons\tNRF\tPBC1\tPBC2\n' > "$OUT/library_complexity_10pct.tsv"
for sample in "${samples[@]}"; do cat "$TMPBASE/${sample}.result" >> "$OUT/library_complexity_10pct.tsv"; done
rm -r "$TMPBASE"
cat "$OUT/library_complexity_10pct.tsv"
