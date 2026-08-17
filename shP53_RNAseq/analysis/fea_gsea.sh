#!/usr/bin/env bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=shp53_fea_gsea
#SBATCH -t 04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH -N 1
#SBATCH --output=shared/log/fea_gsea-%j.out
#SBATCH --error=shared/log/fea_gsea-%j.err

set -euo pipefail

cd /projects/b1042/LauberthLab/BenFolder

module load anaconda3
source /software/anaconda3/2018.12/etc/profile.d/conda.sh
conda activate pydeseq2

export MPLCONFIGDIR="/projects/b1042/LauberthLab/BenFolder/shared/tmp/matplotlib_${SLURM_JOB_ID:-local}"
mkdir -p "$MPLCONFIGDIR" shared/log

# Select a PyDESeq2 comparison with the first argument, or use the default
# comparison produced by shP53_RNAseq/analysis/pydeseq2.sh.
comparison="${1:-shp53_vs_shLacZ_0hr}"
results="${DEA_RESULTS:-shared/rna_seq_dea/${comparison}/significant_results.csv}"
output_root="${DEA_OUTPUT_DIR:-shared/rna_seq_dea/${comparison}}"

if [[ ! -s "$results" ]]; then
    echo "ERROR: DESeq2 results file not found or empty: $results" >&2
    exit 1
fi

if ! python -c 'import pandas, gseapy' >/dev/null 2>&1; then
    echo "ERROR: this environment needs both pandas and gseapy." >&2
    echo "Install gseapy in the pydeseq2 environment, then resubmit:" >&2
    echo "  conda install -n pydeseq2 -c conda-forge -c bioconda gseapy" >&2
    exit 1
fi

python -u shP53_RNAseq/analysis/fea.py \
    --results "$results" \
    --output-dir "$output_root/fea" \
    --padj 0.05 \
    --min-log2fc 1.0 \
    --direction both

python -u shP53_RNAseq/analysis/gsea.py \
    --results "$results" \
    --outdir "$output_root/gsea_out" \
    --padj 0.05 \
    --min-log2fc 1.0 \
    --gene-sets KEGG_2021_Human

echo "FEA output: $output_root/fea"
echo "GSEA output: $output_root/gsea_out"
