#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 00:45:00
#SBATCH --mem=20G
#SBATCH -n 2
#SBATCH -N 1

#SBATCH --cpus-per-task=30
#SBATCH --output=log/job-%j.out
cd /projects/b1042/LauberthLab/BenFolder || exit 1

module load anaconda3
source /software/anaconda3/2018.12/etc/profile.d/conda.sh
conda activate pydeseq2
module load samtools
module load python
module load subread/2.0.3
module load numpy



python3 scripts/differential_expression/deseq2test.py
