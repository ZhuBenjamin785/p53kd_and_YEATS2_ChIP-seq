#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 10:00:00
#SBATCH --mem=32G
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --cpus-per-task=32

set -euo pipefail

cd /projects/b1042/LauberthLab/BenFolder || exit 1

module load anaconda3
source /software/anaconda3/2018.12/etc/profile.d/conda.sh
conda activate /gpfs/home/nqp9093/.conda/envs/sratoolkit
module load MACS3/3.0.2
module load samtools
module load deeptools/3.5.6

mkdir -p fastqchip tmp
  awk 'NF {print $1}' tmp/SRR_Acc_List.txt | while read -r srr; do
      fasterq-dump "$srr" --split-files --threads 4 --temp tmp --outdir fastqchip || exit 1
done
