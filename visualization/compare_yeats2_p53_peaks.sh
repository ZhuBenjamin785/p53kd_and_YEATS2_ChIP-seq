#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=yeats2_p53_overlap
#SBATCH -t 01:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

set -euo pipefail
cd /gpfs/projects/b1042/LauberthLab/BenFolder

module load python
python3 scripts/visualization/compare_yeats2_p53_peaks.py --permutations 10000
