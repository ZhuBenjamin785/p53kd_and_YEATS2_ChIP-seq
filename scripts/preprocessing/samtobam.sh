#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 00:45:00
#SBATCH --mem=20G
#SBATCH -n 2
#SBATCH -N 1

#SBATCH --cpus-per-task=30


cd /projects/b1042/LauberthLab/BenFolder || exit 1


module load samtools


samtools index finalresult.sorted.bam
    

