#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 03:00:00
#SBATCH --mem=32G
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=mof_h4_loss_overlap
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

cd /gpfs/projects/b1042/LauberthLab/BenFolder || exit 1

module load anaconda3/2022.05
module load bedtools/2.31.1
source /software/anaconda3/2022.05/etc/profile.d/conda.sh

python_script="mof_h4k16ac_loss_overlap.py"
annotation_script="chipseekerannotation.r"
output_dir="mof_h4k16ac_loss_overlap"
mof_diffbind="mof_macs3_results/diffbind_results/DiffBind_all_peaks.csv"
h4k16ac_diffbind="diffbind_results/DiffBind_all_peaks.csv"
chrom_sizes="p53KD_YEATS2KD_overlap_analysis/beds/chrom.sizes"

for file in "$python_script" "$annotation_script" "$mof_diffbind" \
            "$h4k16ac_diffbind" "$chrom_sizes"; do
    [[ -s "$file" ]] || { echo "Missing or empty input: $file" >&2; exit 1; }
done
command -v bedtools >/dev/null 2>&1 || { echo "bedtools is unavailable" >&2; exit 1; }

mkdir -p "$output_dir/beds" "$output_dir/tables" "$output_dir/plots" \
             "$output_dir/matplotlib_cache"
export MPLCONFIGDIR="$PWD/$output_dir/matplotlib_cache"

conda activate pybw
python -E "$python_script" prepare \
    --mof "$mof_diffbind" \
    --h4k16ac "$h4k16ac_diffbind" \
    --output "$output_dir"

bedtools sort \
    -i "$output_dir/beds/MOF_loss_unsorted.bed" \
    > "$output_dir/beds/MOF_loss.bed"
bedtools sort \
    -i "$output_dir/beds/H4K16ac_loss_unsorted.bed" \
    > "$output_dir/beds/H4K16ac_loss.bed"

bedtools intersect -nonamecheck -u \
    -a "$output_dir/beds/MOF_loss.bed" \
    -b "$output_dir/beds/H4K16ac_loss.bed" \
    > "$output_dir/beds/MOF_loss_overlapping_H4K16ac_loss.bed"
bedtools intersect -nonamecheck -u \
    -a "$output_dir/beds/H4K16ac_loss.bed" \
    -b "$output_dir/beds/MOF_loss.bed" \
    > "$output_dir/beds/H4K16ac_loss_overlapping_MOF_loss.bed"
bedtools intersect -nonamecheck -wa -wb \
    -a "$output_dir/beds/MOF_loss.bed" \
    -b "$output_dir/beds/H4K16ac_loss.bed" \
    > "$output_dir/tables/shared_loss_pairs.tsv"

python -E "$python_script" build-shared --output "$output_dir"

chipseeker_rscript=/gpfs/home/nqp9093/.conda/envs/chipseeker/bin/Rscript
[[ -x "$chipseeker_rscript" ]] || { echo "Missing ChIPseeker Rscript" >&2; exit 1; }
if [[ -s "$output_dir/beds/shared_loss_intervals.bed" ]]; then
    "$chipseeker_rscript" --vanilla "$annotation_script" \
        "$output_dir/beds/shared_loss_intervals.bed"
else
    echo "No direct shared-loss intervals; skipping ChIPseeker annotation."
fi

python -E "$python_script" summarize \
    --output "$output_dir" \
    --chrom-sizes "$chrom_sizes" \
    --permutations 1000 \
    --seed 12345
