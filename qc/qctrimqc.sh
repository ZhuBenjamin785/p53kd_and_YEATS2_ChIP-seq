#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 05:00:00
#SBATCH --mem=25G
#SBATCH -n 2
#SBATCH -N 1

#SBATCH --cpus-per-task=16


cd /projects/b1042/LauberthLab/BenFolder || exit 1

module load anaconda3
source /software/anaconda3/2018.12/etc/profile.d/conda.sh
conda activate rseqc_env
module load samtools
module load hisat2
module load fastqc
module load multiqc
module load python
module load TrimGalore/0.6.10
module load subread/2.0.3



python scripts/qc/qctrimqc.py
