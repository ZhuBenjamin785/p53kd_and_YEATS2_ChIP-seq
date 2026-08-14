#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 03:00:00
#SBATCH --mem=32G
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=mof_h4_reciprocal
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

cd /gpfs/projects/b1042/LauberthLab/BenFolder || exit 1

python_bin=/gpfs/home/nqp9093/.conda/envs/pybw/bin/python
output_dir="mof_h4k16ac_reciprocal_signal"
mkdir -p "$output_dir/matplotlib_cache"
export MPLCONFIGDIR="$PWD/$output_dir/matplotlib_cache"

"$python_bin" mof_h4k16ac_reciprocal_signal.py \
    --mof-diffbind mof_macs3_results/diffbind_results/DiffBind_all_peaks.csv \
    --h4k16ac-diffbind diffbind_results/DiffBind_all_peaks.csv \
    --mof-control \
        mof_macs3_results/bamcompare/Scr_MOF_rep1_ChIP_vs_Input_log2.bw \
        mof_macs3_results/bamcompare/Scr_MOF_rep2_ChIP_vs_Input_log2.bw \
    --mof-kd \
        mof_macs3_results/bamcompare/p53sh_MOF_rep1_ChIP_vs_Input_log2.bw \
        mof_macs3_results/bamcompare/p53sh_MOF_rep2_ChIP_vs_Input_log2.bw \
    --h4-control \
        bamcompare/Scramble_rep1_ChIP_vs_Input_log2.bw \
        bamcompare/Scramble_rep2_ChIP_vs_Input_log2.bw \
    --h4-kd \
        bamcompare/P53KD_rep1_ChIP_vs_Input_log2.bw \
        bamcompare/P53KD_rep2_ChIP_vs_Input_log2.bw \
    --chrom-sizes p53KD_YEATS2KD_overlap_analysis/beds/chrom.sizes \
    --output "$output_dir" \
    --permutations 1000 \
    --seed 12345 \
    --flank 5000 \
    --bins 100
