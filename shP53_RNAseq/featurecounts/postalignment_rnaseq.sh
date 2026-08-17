#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=rnaseq_counts
#SBATCH -t 12:00:00
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH -N 1
#SBATCH --cpus-per-task=16

set -euo pipefail

cd /gpfs/projects/b1042/LauberthLab/BenFolder

module load anaconda3
source /software/anaconda3/2018.12/etc/profile.d/conda.sh
conda activate rseqc_env
module load samtools
module load subread/2.0.3

bam_dir="shared/BAMfiles/rnaseq_bamfiles/human"
annotation="shared/beds/all.bed"
output_dir="shared/rna_seq_featurecounts"
strandedness_dir="${output_dir}/strandedness"
reference_dir="${output_dir}/reference"
rseqc_bed="${reference_dir}/all_rseqc.bed12"
featurecounts_saf="${reference_dir}/all_exons.saf"
count_file="${output_dir}/rna_seq_featureCounts.txt"
clean_count_file="${output_dir}/rna_seq_featureCounts_cleaned.txt"
threads="${SLURM_CPUS_PER_TASK:-16}"
FEATURECOUNTS_STRAND="${FEATURECOUNTS_STRAND:-0}"

if [[ ! -d "$bam_dir" ]]; then
    echo "ERROR: bam directory not found: $bam_dir" >&2
    exit 1
fi
if [[ ! -s "$annotation" ]]; then
    echo "ERROR: annotation file not found or empty: $annotation" >&2
    exit 1
fi
if ! command -v infer_experiment.py >/dev/null 2>&1; then
    echo "ERROR: infer_experiment.py is not available in the active environment." >&2
    exit 1
fi
if ! command -v featureCounts >/dev/null 2>&1; then
    echo "ERROR: featureCounts is not available in PATH." >&2
    exit 1
fi

mkdir -p "$strandedness_dir" "$reference_dir"

# Keep each SRR as an independent biological replicate. No BAM merging is
# performed; the metadata maps each SRR to its experimental condition.
mapfile -t bams < <(find "$bam_dir" -maxdepth 1 -type f -name 'SRR*.sorted.bam' | sort)
if ((${#bams[@]} == 0)); then
    echo "ERROR: no SRR BAM files found in $bam_dir" >&2
    exit 1
fi
for input_bam in "${bams[@]}"; do
    if [[ ! -s "$input_bam" || ! -s "$input_bam.bai" ]]; then
        echo "ERROR: missing BAM or BAM index: $input_bam" >&2
        exit 1
    fi
done

awk -F '\t' 'BEGIN { OFS="\t" }
    $8 == "exon" && NF >= 6 {
        name = ($4 == "" || $4 == ".") ? "exon_" NR : $4
        print $1, $2, $3, name, 0, $6, $2, $3, 0, 1, $3 - $2, 0
    }' "$annotation" > "$rseqc_bed"

awk -F '\t' 'BEGIN { OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand" }
    $8 == "exon" && NF >= 6 {
        name = ($4 == "" || $4 == ".") ? "exon_" NR : $4
        print name, $1, $2 + 1, $3, $6
    }' "$annotation" > "$featurecounts_saf"

if [[ ! -s "$rseqc_bed" || ! -s "$featurecounts_saf" ]]; then
    echo "ERROR: failed to derive RNA-seq annotation files from $annotation" >&2
    exit 1
fi

for bam in "${bams[@]}"; do
    sample=$(basename "$bam" .sorted.bam)
    echo "Running RSeQC strandedness for $sample"
    infer_experiment.py -r "$rseqc_bed" -i "$bam" \
        > "$strandedness_dir/${sample}_infer_experiment.txt"
done

echo "Running featureCounts on ${#bams[@]} single-end bam files"
featureCounts \
    -T "$threads" \
    -a "$featurecounts_saf" \
    -F SAF \
    -s "$FEATURECOUNTS_STRAND" \
    -o "$count_file" \
    "${bams[@]}"

python3 shP53_RNAseq/featurecounts/featureCountscleanup.py "$count_file" "$clean_count_file" "$annotation"

echo "RNA-seq strandedness reports: $strandedness_dir"
echo "Raw featureCounts output: $count_file"
echo "Cleaned featureCounts output: $clean_count_file"
