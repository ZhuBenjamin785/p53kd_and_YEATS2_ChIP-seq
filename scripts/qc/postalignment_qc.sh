#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH --mem=16G
#SBATCH -n 2
#SBATCH -N 1

#SBATCH --cpus-per-task=40


cd /projects/b1042/LauberthLab/BenFolder || exit 1


module load bedops/2.4.40 
module load kentUtils/302
module load anaconda3
source /software/anaconda3/2018.12/etc/profile.d/conda.sh
conda activate rseqc_env

python scripts/alignment/alignment.py

for species in human dm6; do
    for bam in BAMfiles/$species/*.bam; do
        [[ "$bam" == *_sorted.bam ]] && continue

        outfile="BAMfiles/$species/$(basename "${bam%.bam}")_sorted.bam"

        samtools sort -@ 40 -m 4G -o "$outfile" "$bam"
    done
done
