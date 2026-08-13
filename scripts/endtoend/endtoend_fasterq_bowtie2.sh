#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=fasterq_bowtie2_align
#SBATCH -t 24:00:00
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH -N 1
#SBATCH --cpus-per-task=18
#SBATCH --output=log/slurm-%j.out

set -euo pipefail
cd /gpfs/projects/b1042/LauberthLab/BenFolder

module load bowtie2/2.5.4
module load samtools/1.14
python3 -u scripts/alignment/alignment_fasterq_bowtie2.py
