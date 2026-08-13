#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 10:00:00
#SBATCH --mem=64G
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=bamcompare_p53_0hr
#SBATCH --output=slurm-%j.out

set -euo pipefail

cd /projects/b1042/LauberthLab/BenFolder || exit 1

unset PYTHONPATH PYTHONHOME
module load deeptools/3.5.6

deeptools_bin=/software/deeptools/3.5.6/bin
[[ -x "$deeptools_bin/bamCompare" ]] || { echo "Missing deepTools bamCompare" >&2; exit 1; }
[[ -x "$deeptools_bin/bigwigAverage" ]] || { echo "Missing deepTools bigwigAverage" >&2; exit 1; }
[[ -x "$deeptools_bin/bigwigCompare" ]] || { echo "Missing deepTools bigwigCompare" >&2; exit 1; }

input_dir="fastqchip_macs3_results/normalized_bams"
output_dir="fastqchip_macs3_results/bamcompare"
mkdir -p "$output_dir"

compare() {
    local chip="$1"
    local input="$2"
    local output="$3"

    [[ -s "$input_dir/$chip" ]] || { echo "Missing ChIP bam: $input_dir/$chip" >&2; exit 1; }
    [[ -s "$input_dir/$input" ]] || { echo "Missing Input bam: $input_dir/$input" >&2; exit 1; }

    "$deeptools_bin/bamCompare" \
        -b1 "$input_dir/$chip" \
        -b2 "$input_dir/$input" \
        --operation log2 \
        --scaleFactorsMethod None \
        --normalizeUsing BPM \
        --binSize 25 \
        --numberOfProcessors "${SLURM_CPUS_PER_TASK:-32}" \
        --outFileFormat bigwig \
        -o "$output_dir/$output"
}

compare \
    p53_0hr.spikein_normalized.bam \
    input_0hr.spikein_normalized.bam \
    p53_0hr_ChIP_vs_input_log2.bw
