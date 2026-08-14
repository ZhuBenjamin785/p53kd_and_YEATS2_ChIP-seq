#!/usr/bin/env bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=mof_align
#SBATCH -t 24:00:00
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH -N 1
#SBATCH --cpus-per-task=18

set -euo pipefail

cd /projects/b1042/LauberthLab/BenFolder || exit 1

module load bowtie2/2.5.4
module load samtools/1.14

python3 -u alignment.py
