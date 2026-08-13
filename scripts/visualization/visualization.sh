#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 05:00:00
#SBATCH --mem=25G
#SBATCH -N 1
#SBATCH --cpus-per-task=16

set -euo pipefail

cd /projects/b1042/LauberthLab/BenFolder || exit 1

rscript=/gpfs/home/nqp9093/.conda/envs/chipseeker/bin/Rscript
if [[ ! -x "$rscript" ]]; then
  echo "Rscript not found or not executable: $rscript" >&2
  exit 1
fi

"$rscript" scripts/visualization/chipseekervisualization_fasterqchip.r
