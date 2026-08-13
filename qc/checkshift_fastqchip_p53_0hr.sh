#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 18:00:00
#SBATCH --mem=32G
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=p53_0hr_signal
#SBATCH --output=log/slurm-%j.out
#SBATCH --error=log/slurm-%j.err

set -euo pipefail
cd /gpfs/projects/b1042/LauberthLab/BenFolder

source /software/anaconda3/2022.05/etc/profile.d/conda.sh
unset PYTHONPATH PYTHONHOME
conda activate pybw
env -u PYTHONPATH -u PYTHONHOME python -E scripts/qc/pybw_fastqchip_p53_0hr.py

echo "Finished. Results are in fastqchip_p53_0hr_signal_results/"
