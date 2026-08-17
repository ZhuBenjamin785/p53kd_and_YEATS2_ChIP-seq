#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=rnaseq_qc_trim
#SBATCH -t 12:00:00
#SBATCH --mem=25G
#SBATCH --ntasks=1
#SBATCH -N 1
#SBATCH --cpus-per-task=16
set -euo pipefail
cd /gpfs/projects/b1042/LauberthLab/BenFolder
module load anaconda3
source /software/anaconda3/2018.12/etc/profile.d/conda.sh
conda activate rseqc_env
module load fastqc multiqc TrimGalore/0.6.10
python3 -u shP53_RNAseq/qctrimqc/qctrimqc_rnaseq.py
