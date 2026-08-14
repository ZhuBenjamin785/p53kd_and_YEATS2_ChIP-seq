#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 3:00:00
#SBATCH --mem=32G
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=plots_MOF
#SBATCH --output=log/slurm-%j.out
#SBATCH --error=log/slurm-%j.err

set -euo pipefail

cd /gpfs/projects/b1042/LauberthLab/BenFolder || exit 1
mkdir -p log

module load anaconda3/2022.05
source /software/anaconda3/2022.05/etc/profile.d/conda.sh

conda activate chipseeker
Rscript --vanilla scripts/chipseq_summary/chipseq_summary_plots.r MOF

conda activate pybw
python -E scripts/chipseq_summary/chipseq_summary_plots.py MOF
