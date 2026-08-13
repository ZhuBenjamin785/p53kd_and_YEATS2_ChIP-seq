#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 3:00:00
#SBATCH --mem=32G
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=plots_p53_0hr
#SBATCH --output=log/slurm-%j.out
#SBATCH --error=log/slurm-%j.err
set -euo pipefail
cd /gpfs/projects/b1042/LauberthLab/BenFolder
module load anaconda3/2022.05
source /software/anaconda3/2022.05/etc/profile.d/conda.sh
conda activate chipseeker
package_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
Rscript --vanilla "${package_dir}/chipseq_summary_plots.r" p53_0hr
conda activate pybw
python -E "${package_dir}/chipseq_summary_plots.py" p53_0hr
