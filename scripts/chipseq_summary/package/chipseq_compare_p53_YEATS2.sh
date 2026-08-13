#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 02:00:00
#SBATCH --mem=16G
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=compare_p53_YEATS2
#SBATCH --output=log/slurm-%j.out
#SBATCH --error=log/slurm-%j.err

set -euo pipefail
cd /gpfs/projects/b1042/LauberthLab/BenFolder
module load anaconda3/2022.05
source /software/anaconda3/2022.05/etc/profile.d/conda.sh
conda activate pybw
package_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python -E "${package_dir}/chipseq_compare_p53_YEATS2.py"
