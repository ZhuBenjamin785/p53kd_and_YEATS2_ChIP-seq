#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 00:45:00
#SBATCH --mem=20G
#SBATCH -n 2
#SBATCH -N 1

#SBATCH --cpus-per-task=30
#SBATCH --output=shared/log/job-%j.out

set -euo pipefail

cd /projects/b1042/LauberthLab/BenFolder || exit 1

module load anaconda3
source /software/anaconda3/2018.12/etc/profile.d/conda.sh
conda activate pydeseq2

# Do not load another Python/module stack after activating pydeseq2; that
# replaces the conda interpreter and hides pandas/PyDESeq2.
export MPLCONFIGDIR="/projects/b1042/LauberthLab/BenFolder/shared/tmp/matplotlib_${SLURM_JOB_ID:-local}"
mkdir -p "$MPLCONFIGDIR"
python -c 'import pandas, numpy, pydeseq2'
python -u shP53_RNAseq/analysis/dea.py --deseq2 \
    --counts shared/rna_seq_featurecounts/rna_seq_featureCounts_cleaned.txt \
    --metadata rna_seq_metadata.tsv \
    --case shp53_0hr --control shLacZ_0hr \
    --output-dir shared/rna_seq_dea/shp53_vs_shLacZ_0hr
