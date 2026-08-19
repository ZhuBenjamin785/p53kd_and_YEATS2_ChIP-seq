#!/usr/bin/env bash
set -euo pipefail
usage(){ echo "Usage: $0 --human-bam FILE --dm6-bam FILE --sample NAME --outdir DIR" >&2; exit 2; }
human=; dm6=; sample=; outdir=
while (($#)); do case "$1" in --human-bam) human=$2; shift 2;; --dm6-bam) dm6=$2; shift 2;; --sample) sample=$2; shift 2;; --outdir) outdir=$2; shift 2;; *) usage;; esac; done
[[ -n "$human" && -n "$dm6" && -n "$sample" && -n "$outdir" ]] || usage
command -v samtools >/dev/null || { echo 'samtools is required' >&2; exit 1; }
mkdir -p "$outdir"; tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
samtools view -@ 2 -f 3 -F 3852 -q 30 "$human" | awk '{print $1}' | sort -u > "$tmp/human.qnames"
samtools view -@ 2 -f 3 -F 3852 -q 30 "$dm6" | awk '{print $1}' | sort -u > "$tmp/dm6.qnames"
shared=$(comm -12 "$tmp/human.qnames" "$tmp/dm6.qnames" | tee "$outdir/${sample}.shared_qnames.tsv" | wc -l)
human_n=$(wc -l < "$tmp/human.qnames"); dm6_n=$(wc -l < "$tmp/dm6.qnames"); smaller=$(( human_n < dm6_n ? human_n : dm6_n ))
fraction=$(awk -v s="$shared" -v n="$smaller" 'BEGIN{if(n==0) print "0.000000"; else printf "%.6f",s/n}')
printf 'sample\thuman_primary_proper_fragments_mapq30\tdm6_primary_proper_fragments_mapq30\tshared_qnames_mapq30\tshared_fraction_of_smaller_species_set\n%s\t%s\t%s\t%s\t%s\n' "$sample" "$human_n" "$dm6_n" "$shared" "$fraction" > "$outdir/${sample}.species_cross_mapping.tsv"
echo "Wrote $outdir/${sample}.species_cross_mapping.tsv"
