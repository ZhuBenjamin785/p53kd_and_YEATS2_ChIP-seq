#!/usr/bin/env bash
#SBATCH --job-name=validate_plot_scripts
#SBATCH --partition=genomics
#SBATCH --account=b1042
#SBATCH --output=shared/biological_consensus_repaired/plotscripts/log_%x_%j.out
#SBATCH --error=shared/biological_consensus_repaired/plotscripts/log_%x_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

set -euo pipefail
cd /gpfs/projects/b1042/LauberthLab/BenFolder
python -m py_compile shared/biological_consensus_repaired/plotscripts/*.py
echo "All plotting scripts compiled successfully."
