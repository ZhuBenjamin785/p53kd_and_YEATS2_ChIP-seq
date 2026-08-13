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
conda activate chipseeker
module load samtools
module load hisat2
module load fastqc
module load multiqc
module load python
module load TrimGalore/0.6.10
module load subread/2.0.3
module load anaconda3/2022.05
module load bedtools/2.31.1


bedtools intersect \
    -a macs3_results_p53kd/peaks/Scramble_H4K16ac_rep1/Scramble_H4K16ac_rep1_peaks.broadPeak \
    -b macs3_results_p53kd/peaks/Scramble_H4K16ac_rep2/Scramble_H4K16ac_rep2_peaks.broadPeak \
    > scramble_consensus.bed

bedtools intersect \
    -a  macs3_results_p53kd/peaks/p53KD_H4K16ac_rep1/p53KD_H4K16ac_rep1_peaks.broadPeak \
    -b  macs3_results_p53kd/peaks/p53KD_H4K16ac_rep2/p53KD_H4K16ac_rep2_peaks.broadPeak \
    > p53kd_consensus.bed

bedtools intersect \
    -a scramble_consensus.bed \
    -b p53kd_consensus.bed \
    > final_consensus.bed