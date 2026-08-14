#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=mof_diffbind
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

trap 'status=$?; echo "DiffBind exited with status ${status}" >&2; exit ${status}' EXIT



cd /projects/b1042/LauberthLab/BenFolder || exit 1

source /software/anaconda3/2022.05/etc/profile.d/conda.sh
module load anaconda3/2022.05

conda activate chipseeker




echo "Using R: $(command -v Rscript)"
Rscript --vanilla -e 'if (!requireNamespace("DiffBind", quietly=TRUE)) stop("DiffBind is not installed in the active R environment")'
exec Rscript --vanilla diffbindMOF.r
