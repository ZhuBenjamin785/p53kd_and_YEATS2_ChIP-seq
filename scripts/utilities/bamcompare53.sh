#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 10:00:00
#SBATCH --mem=64G
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=bamcompare_p53kd
#SBATCH --output=slurm-%j.out

set -euo pipefail

cd /projects/b1042/LauberthLab/BenFolder || exit 1

unset PYTHONPATH PYTHONHOME
module load deeptools/3.5.6

deeptools_bin=/software/deeptools/3.5.6/bin
[[ -x "$deeptools_bin/bamCompare" ]] || { echo "Missing deepTools bamCompare" >&2; exit 1; }
[[ -x "$deeptools_bin/bigwigAverage" ]] || { echo "Missing deepTools bigwigAverage" >&2; exit 1; }
[[ -x "$deeptools_bin/bigwigCompare" ]] || { echo "Missing deepTools bigwigCompare" >&2; exit 1; }

input_dir="macs3_results_p53kd/normalized_bam"
output_dir="bamcompare"
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
    P53_H4K16ac_1_S0_L001.spikein_normalized.bam \
    MutP53_Input1_S0_L001.spikein_normalized.bam \
    P53KD_rep1_ChIP_vs_Input_log2.bw

compare \
    P53_H4K16ac_2_S0_L001.spikein_normalized.bam \
    MutP53_Input2_S0_L001.spikein_normalized.bam \
    P53KD_rep2_ChIP_vs_Input_log2.bw

compare \
    Scr_H4K16ac_1_S0_L001.spikein_normalized.bam \
    Scr_Input1_S0_L001.spikein_normalized.bam \
    Scramble_rep1_ChIP_vs_Input_log2.bw

compare \
    Scr_H4K16ac_2_S0_L001.spikein_normalized.bam \
    Scr_Input2_S0_L001.spikein_normalized.bam \
    Scramble_rep2_ChIP_vs_Input_log2.bw
