#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=rnaseq_align
#SBATCH -t 24:00:00
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH -N 1
#SBATCH --cpus-per-task=18
set -euo pipefail
cd /gpfs/projects/b1042/LauberthLab/BenFolder
module load hisat2 samtools/1.14
python3 -u shP53_RNAseq/alignment/rnaseqSE_alignment.py
