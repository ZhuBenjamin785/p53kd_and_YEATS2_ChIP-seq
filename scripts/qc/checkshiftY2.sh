#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 18:00:00
#SBATCH --mem=32G
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=h4k16ac_shift
#SBATCH --output=slurm-%j.out

set -euo pipefail



cd /projects/b1042/LauberthLab/BenFolder || exit 1

source /software/anaconda3/2022.05/etc/profile.d/conda.sh
unset PYTHONPATH PYTHONHOME
conda activate pybw
env -u PYTHONPATH -u PYTHONHOME python -E scripts/qc/pybwY2.py

echo "Finished. Results are in shift_resultsY2/"
