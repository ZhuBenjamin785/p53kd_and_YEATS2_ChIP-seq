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

bam_dir="rnaseq_bamfiles/human"
annotation="beds/all.bed"
output_dir="rna_seq_featurecounts"
strandedness_dir="${output_dir}/strandedness"
reference_dir="${output_dir}/reference"
merged_bam_dir="${output_dir}/merged_bams"
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

mkdir -p "$strandedness_dir" "$reference_dir" "$merged_bam_dir"

sample_names=(
    "GSM2746544_SW480_shLacZ_0hr"
    "GSM2746545_SW480_shLacZ_16hr_TNF"
    "GSM2746546_SW480_shp53_0hr"
    "GSM2746547_SW480_shp53_16hr_TNF"
)
sample_runs=(
    "SRR5944073 SRR5944074"
    "SRR5944075 SRR5944076"
    "SRR5944077 SRR5944078"
    "SRR5944079 SRR5944080"
)

bams=()
for index in "${!sample_names[@]}"; do
    sample="${sample_names[$index]}"
    merged_bam="$merged_bam_dir/${sample}.sorted.bam"
    read -r run1 run2 <<< "${sample_runs[$index]}"
    input1="$bam_dir/${run1}.sorted.bam"
    input2="$bam_dir/${run2}.sorted.bam"
    for input_bam in "$input1" "$input2"; do
        if [[ ! -s "$input_bam" || ! -s "$input_bam.bai" ]]; then
            echo "ERROR: missing bam or bam index: $input_bam" >&2
            exit 1
        fi
    done
    echo "Merging $run1 and $run2 into $sample"
    samtools merge -f -@ "$threads" "$merged_bam" "$input1" "$input2"
    samtools index -@ "$threads" "$merged_bam"
    bams+=("$merged_bam")
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

python3 scripts/preprocessing/featureCountscleanup.py "$count_file" "$clean_count_file" "$annotation"

echo "RNA-seq strandedness reports: $strandedness_dir"
echo "Raw featureCounts output: $count_file"
echo "Cleaned featureCounts output: $clean_count_file"
