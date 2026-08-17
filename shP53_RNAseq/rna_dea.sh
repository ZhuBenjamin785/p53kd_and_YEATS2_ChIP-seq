#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=rnaseq_dea
#SBATCH -t 02:00:00
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --output=shared/log/job-%j.out

set -euo pipefail
cd /gpfs/projects/b1042/LauberthLab/BenFolder

module load anaconda3
source /software/anaconda3/2018.12/etc/profile.d/conda.sh
conda activate pydeseq2
if (( $# == 0 )); then
    python3 -u shP53_RNAseq/analysis/dea.py --case shp53_0hr --control shLacZ_0hr \
        --counts shared/rna_seq_featurecounts/rna_seq_featureCounts_cleaned.txt \
        --metadata rna_seq_metadata.tsv --output-dir shared/rna_seq_dea/exploratory_shp53_vs_shLacZ_0hr
    python3 -u shP53_RNAseq/analysis/dea.py --case shp53_16hr_TNF --control shLacZ_16hr_TNF \
        --counts shared/rna_seq_featurecounts/rna_seq_featureCounts_cleaned.txt \
        --metadata rna_seq_metadata.tsv --output-dir shared/rna_seq_dea/exploratory_shp53_vs_shLacZ_16hr_TNF
    python3 -u shP53_RNAseq/analysis/dea.py --case shLacZ_16hr_TNF --control shLacZ_0hr \
        --counts shared/rna_seq_featurecounts/rna_seq_featureCounts_cleaned.txt \
        --metadata rna_seq_metadata.tsv --output-dir shared/rna_seq_dea/exploratory_TNF_shLacZ
    python3 -u shP53_RNAseq/analysis/dea.py --case shp53_16hr_TNF --control shp53_0hr \
        --counts shared/rna_seq_featurecounts/rna_seq_featureCounts_cleaned.txt \
        --metadata rna_seq_metadata.tsv --output-dir shared/rna_seq_dea/exploratory_TNF_shp53
else
    python3 -u shP53_RNAseq/analysis/dea.py \
        --counts shared/rna_seq_featurecounts/rna_seq_featureCounts_cleaned.txt \
        --metadata rna_seq_metadata.tsv "$@"
fi
