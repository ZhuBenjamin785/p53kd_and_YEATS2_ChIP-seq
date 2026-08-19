#!/usr/bin/env bash
#SBATCH --job-name=consensus_plots
#SBATCH --partition=genomics
#SBATCH --account=b1042
#SBATCH --output=shared/biological_consensus_repaired/plotscripts/log_%x_%j.out
#SBATCH --error=shared/biological_consensus_repaired/plotscripts/log_%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

set -euo pipefail
cd /gpfs/projects/b1042/LauberthLab/BenFolder
export MPLCONFIGDIR=/tmp/mpl_consensus_repair_${SLURM_JOB_ID}
mkdir -p "$MPLCONFIGDIR"
python shared/biological_consensus_repaired/plotscripts/make_all_plots.py
