#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=fasterqchip_dm6_align
#SBATCH -t 24:00:00
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=log/slurm-%j.out

set -euo pipefail
cd "${FASTQCHIP_PROJECT_DIR:-/gpfs/projects/b1042/LauberthLab/BenFolder}"
module load bowtie2/2.5.4
module load samtools/1.14
python3 -u scripts/alignment/alignment_fasterqchip_dm6.py
