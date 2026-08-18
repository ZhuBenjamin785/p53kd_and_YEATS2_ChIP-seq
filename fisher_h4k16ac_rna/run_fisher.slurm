#!/bin/bash
#SBATCH --job-name=H4K16ac_RNA_Fisher
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=H4K16ac_RNA_Fisher_%j.out
#SBATCH --error=H4K16ac_RNA_Fisher_%j.err

set -euo pipefail

BASE="/gpfs/projects/b1042/LauberthLab/BenFolder"
SCRIPT="$BASE/shared/scripts/fisher_h4k16ac_rna/fisher_h4k16ac_rna.py"

# Change these input filenames if needed.
LOSS="$BASE/p53kdH4K16ac/diffbind_results/H4K16ac_lost_annotated.csv"
GAIN="$BASE/p53kdH4K16ac/diffbind_results/H4K16ac_gained_annotated.csv"
RNA="$BASE/shared/rna_seq_dea/shp53_vs_shLacZ_0hr/significant_results.csv"
BACKGROUND="$BASE/shared/rna_chip_integration/p53KD_H4K16ac_vs_RNAseq/ora_matched_universes/universes/all_peaks_eligible_gene_universe.csv"
OUTPUT="$BASE/shared/rna_chip_integration/p53KD_H4K16ac_vs_RNAseq/fisher_exact_results.csv"

# Activate the conda environment containing scipy.
source /software/anaconda3/2022.05/etc/profile.d/conda.sh
conda activate /home/nqp9093/.conda/envs/pydeseq2

# Change these column names if your input headers differ:
# ChIP gene = SYMBOL; RNA gene = gene_name; RNA fold change = log2FoldChange;
# background gene = Gene. Remove the --gain line to skip the reverse test.
python "$SCRIPT" \
  --loss "$LOSS" \
  --gain "$GAIN" \
  --rna "$RNA" \
  --background "$BACKGROUND" \
  --output "$OUTPUT" \
  --chip-gene-column SYMBOL \
  --rna-gene-column gene_name \
  --rna-fold-change-column log2FoldChange \
  --background-gene-column Gene
