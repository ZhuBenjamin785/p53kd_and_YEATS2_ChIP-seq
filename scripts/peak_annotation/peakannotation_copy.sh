#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 18:00:00
#SBATCH --mem=64G
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --cpus-per-task=32

set -euo pipefail

cd /projects/b1042/LauberthLab/BenFolder || exit 1

module load anaconda3
source /software/anaconda3/2018.12/etc/profile.d/conda.sh
conda activate chipseeker
module load MACS3/3.0.2
module load samtools
module load deeptools/3.5.6

inputdir="${inputdir:-BAMfiles/human}"
outdir="${outdir:-macs3_results}"
sample="YEATS2KD_H4K16ac_rep1"

chip_bam="${inputdir}/H4K16ac_YEATS2_shRNA#3_ChIP1_S157_L003_sorted.bam"
input_bam="${inputdir}/YEATS2_shRNA#3_Input1_S153_L003_sorted.bam"


mkdir -p "${outdir}/peaks/${sample}" "${outdir}/tracks"





samtools index "${chip_bam}"
samtools index "${input_bam}"

bamCoverage -b "${chip_bam}" -o "${outdir}/tracks/$(basename "${chip_bam%.bam}").bw" -of bigwig -p 8
bamCoverage -b "${input_bam}" -o "${outdir}/tracks/$(basename "${input_bam%.bam}").bw" -of bigwig -p 8

macs3 callpeak \
  -t "${chip_bam}" \
  -c "${input_bam}" \
  -f BAMPE \
  -g hs \
  -n "${sample}" \
  --outdir "${outdir}/peaks/${sample}" \
  --broad \
  --broad-cutoff 0.1

Rscript scripts/peak_annotation/chipseekerannotation.r \
  "${outdir}/peaks/${sample}/${sample}_peaks.broadPeak"
