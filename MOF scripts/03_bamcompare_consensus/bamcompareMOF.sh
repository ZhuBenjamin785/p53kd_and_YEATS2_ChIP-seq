#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 10:00:00
#SBATCH --mem=64G
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=bamcompare_mof
#SBATCH --output=slurm-%j.out

set -euo pipefail

cd /projects/b1042/LauberthLab/BenFolder || exit 1

unset PYTHONPATH PYTHONHOME
module load deeptools/3.5.6
module load bedtools/2.31.1

deeptools_bin=/software/deeptools/3.5.6/bin
[[ -x "$deeptools_bin/bamCompare" ]] || { echo "Missing deepTools bamCompare" >&2; exit 1; }
[[ -x "$deeptools_bin/bigwigAverage" ]] || { echo "Missing deepTools bigwigAverage" >&2; exit 1; }
[[ -x "$deeptools_bin/bigwigCompare" ]] || { echo "Missing deepTools bigwigCompare" >&2; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "Missing bedtools" >&2; exit 1; }

input_dir="mof_macs3_results/normalized_bam"
output_dir="mof_macs3_results/bamcompare"
peak_dir="mof_macs3_results/peaks"
consensus_dir="mof_macs3_results/consensus_peaks"
mkdir -p "$output_dir" "$consensus_dir"

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

average_replicates() {
    local replicate_1="$1"
    local replicate_2="$2"
    local output="$3"

    [[ -s "$output_dir/$replicate_1" ]] || { echo "Missing replicate bigWig: $output_dir/$replicate_1" >&2; exit 1; }
    [[ -s "$output_dir/$replicate_2" ]] || { echo "Missing replicate bigWig: $output_dir/$replicate_2" >&2; exit 1; }

    "$deeptools_bin/bigwigAverage" \
        -b "$output_dir/$replicate_1" "$output_dir/$replicate_2" \
        --binSize 25 \
        --numberOfProcessors "${SLURM_CPUS_PER_TASK:-32}" \
        -o "$output_dir/$output"
}

make_consensus_peaks() {
    local replicate_1="$1"
    local replicate_2="$2"
    local output="$3"

    [[ -s "$peak_dir/$replicate_1/$replicate_1"_peaks.broadPeak ]] || { echo "Missing replicate peaks: $peak_dir/$replicate_1/$replicate_1"_peaks.broadPeak >&2; exit 1; }
    [[ -s "$peak_dir/$replicate_2/$replicate_2"_peaks.broadPeak ]] || { echo "Missing replicate peaks: $peak_dir/$replicate_2/$replicate_2"_peaks.broadPeak >&2; exit 1; }

    bedtools intersect \
        -nonamecheck \
        -a "$peak_dir/$replicate_1/$replicate_1"_peaks.broadPeak \
        -b "$peak_dir/$replicate_2/$replicate_2"_peaks.broadPeak \
        > "$consensus_dir/$output"

    [[ -s "$consensus_dir/$output" ]] || { echo "Consensus peak file is empty: $consensus_dir/$output" >&2; exit 1; }
}

compare \
    p53sh_MOF_1.spikein_normalized.bam \
    MutP53_Input1_S0_L001.spikein_normalized.bam \
    p53sh_MOF_rep1_ChIP_vs_Input_log2.bw

compare \
    p53sh_MOF_2.spikein_normalized.bam \
    MutP53_Input2_S0_L001.spikein_normalized.bam \
    p53sh_MOF_rep2_ChIP_vs_Input_log2.bw

compare \
    Scr_MOF_1.spikein_normalized.bam \
    Scr_Input1_S0_L001.spikein_normalized.bam \
    Scr_MOF_rep1_ChIP_vs_Input_log2.bw

compare \
    Scr_MOF_2.spikein_normalized.bam \
    Scr_Input2_S0_L001.spikein_normalized.bam \
    Scr_MOF_rep2_ChIP_vs_Input_log2.bw

average_replicates \
    Scr_MOF_rep1_ChIP_vs_Input_log2.bw \
    Scr_MOF_rep2_ChIP_vs_Input_log2.bw \
    Scr_MOF_consensus_ChIP_vs_Input_log2.bw

average_replicates \
    p53sh_MOF_rep1_ChIP_vs_Input_log2.bw \
    p53sh_MOF_rep2_ChIP_vs_Input_log2.bw \
    p53sh_MOF_consensus_ChIP_vs_Input_log2.bw

make_consensus_peaks \
    Scr_MOF_rep1 \
    Scr_MOF_rep2 \
    Scr_MOF_consensus.bed

make_consensus_peaks \
    p53sh_MOF_rep1 \
    p53sh_MOF_rep2 \
    p53sh_MOF_consensus.bed
