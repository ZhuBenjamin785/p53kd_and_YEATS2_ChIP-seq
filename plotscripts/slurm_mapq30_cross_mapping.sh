#!/usr/bin/env bash
#SBATCH --job-name=mapq30_crossmap
#SBATCH --partition=genomics
#SBATCH --account=b1042
#SBATCH --output=shared/biological_consensus_repaired/plotscripts/log_%x_%j.out
#SBATCH --error=shared/biological_consensus_repaired/plotscripts/log_%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail
cd /gpfs/projects/b1042/LauberthLab/BenFolder

# Edit these three values before submission.
HUMAN_BAM="/path/to/human.bam"
DM6_BAM="/path/to/dm6.bam"
SAMPLE="sample_name"
OUTDIR="shared/biological_consensus_repaired/chipseq/qc"

bash shared/biological_consensus_repaired/plotscripts/run_mapq30_cross_mapping.sh \
  --human-bam "$HUMAN_BAM" \
  --dm6-bam "$DM6_BAM" \
  --sample "$SAMPLE" \
  --outdir "$OUTDIR"
