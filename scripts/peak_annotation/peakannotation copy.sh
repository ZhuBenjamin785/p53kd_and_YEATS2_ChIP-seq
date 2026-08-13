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

mkdir -p "${outdir}/peaks" "${outdir}/dedup_bam" "${outdir}/tracks"

H4K16ac_Scrameble_ChIP1_scale="${H4K16ac_Scrameble_ChIP1_scale:-0.785}"
H4K16ac_Scrameble_ChIP2_scale="${H4K16ac_Scrameble_ChIP2_scale:-0.912}"
H4K16ac_YEATS2_ChIP1_scale="${H4K16ac_YEATS2_ChIP1_scale:-1}"
H4K16ac_YEATS2_ChIP2_scale="${H4K16ac_YEATS2_ChIP2_scale:-0.685}"


samtools view -@ 8 -s "34.${H4K16ac_Scrameble_ChIP1_scale}" -b \
  "${inputdir}/H4K16ac_Scrameble_ChIP1_S155_L003_sorted.bam" \
  | samtools sort -@ 8 -o "${outdir}/dedup_bam/H4K16ac_Scrameble_ChIP1_S155_L003_sorted_downsampled.bam"
samtools index "${outdir}/dedup_bam/H4K16ac_Scrameble_ChIP1_S155_L003_sorted_downsampled.bam"

samtools view -@ 8 -s "34.${H4K16ac_Scrameble_ChIP2_scale}" -b \
  "${inputdir}/H4K16ac_Scrameble_ChIP2_S156_L003_sorted.bam" \
  | samtools sort -@ 8 -o "${outdir}/dedup_bam/H4K16ac_Scrameble_ChIP2_S156_L003_sorted_downsampled.bam"
samtools index "${outdir}/dedup_bam/H4K16ac_Scrameble_ChIP2_S156_L003_sorted_downsampled.bam"

samtools view -@ 8 -s "34.${H4K16ac_YEATS2_ChIP1_scale}" -b \
  "${inputdir}/H4K16ac_YEATS2_shRNA#3_ChIP1_S157_L003_sorted.bam" \
  | samtools sort -@ 8 -o "${outdir}/dedup_bam/H4K16ac_YEATS2_shRNA#3_ChIP1_S157_L003_sorted_downsampled.bam"
samtools index "${outdir}/dedup_bam/H4K16ac_YEATS2_shRNA#3_ChIP1_S157_L003_sorted_downsampled.bam"

samtools view -@ 8 -s "34.${H4K16ac_YEATS2_ChIP2_scale}" -b \
  "${inputdir}/H4K16ac_YEATS2_shRNA#3_ChIP2_S158_L003_sorted.bam" \
  | samtools sort -@ 8 -o "${outdir}/dedup_bam/H4K16ac_YEATS2_shRNA#3_ChIP2_S158_L003_sorted_downsampled.bam"
samtools index "${outdir}/dedup_bam/H4K16ac_YEATS2_shRNA#3_ChIP2_S158_L003_sorted_downsampled.bam"

samtools view -@ 8 -s "34.${H4K16ac_Scrameble_ChIP1_scale}" -b \
  "${inputdir}/scramble_Input_rep1_S151_L003_sorted.bam" \
  | samtools sort -@ 8 -o "${outdir}/dedup_bam/scramble_Input_rep1_S151_L003_sorted_downsampled.bam"
samtools index "${outdir}/dedup_bam/scramble_Input_rep1_S151_L003_sorted_downsampled.bam"

samtools view -@ 8 -s "34.${H4K16ac_Scrameble_ChIP2_scale}" -b \
  "${inputdir}/scramble_Input_rep2_S152_L003_sorted.bam" \
  | samtools sort -@ 8 -o "${outdir}/dedup_bam/scramble_Input_rep2_S152_L003_sorted_downsampled.bam"
samtools index "${outdir}/dedup_bam/scramble_Input_rep2_S152_L003_sorted_downsampled.bam"

samtools view -@ 8 -s "34.${H4K16ac_YEATS2_ChIP1_scale}" -b \
  "${inputdir}/YEATS2_shRNA#3_Input1_S153_L003_sorted.bam" \
  | samtools sort -@ 8 -o "${outdir}/dedup_bam/YEATS2_shRNA#3_Input1_S153_L003_sorted_downsampled.bam"
samtools index "${outdir}/dedup_bam/YEATS2_shRNA#3_Input1_S153_L003_sorted_downsampled.bam"

samtools view -@ 8 -s "34.${H4K16ac_YEATS2_ChIP2_scale}" -b \
  "${inputdir}/YEATS2_shRNA#3_Input2_S154_L003_sorted.bam" \
  | samtools sort -@ 8 -o "${outdir}/dedup_bam/YEATS2_shRNA#3_Input2_S154_L003_sorted_downsampled.bam"
samtools index "${outdir}/dedup_bam/YEATS2_shRNA#3_Input2_S154_L003_sorted_downsampled.bam"


downsampled_bams=(
  "${outdir}/dedup_bam/H4K16ac_Scrameble_ChIP1_S155_L003_sorted_downsampled.bam"
  "${outdir}/dedup_bam/scramble_Input_rep1_S151_L003_sorted_downsampled.bam"
  "${outdir}/dedup_bam/H4K16ac_Scrameble_ChIP2_S156_L003_sorted_downsampled.bam"
  "${outdir}/dedup_bam/scramble_Input_rep2_S152_L003_sorted_downsampled.bam"
  "${outdir}/dedup_bam/H4K16ac_YEATS2_shRNA#3_ChIP1_S157_L003_sorted_downsampled.bam"
  "${outdir}/dedup_bam/YEATS2_shRNA#3_Input1_S153_L003_sorted_downsampled.bam"
  "${outdir}/dedup_bam/H4K16ac_YEATS2_shRNA#3_ChIP2_S158_L003_sorted_downsampled.bam"
  "${outdir}/dedup_bam/YEATS2_shRNA#3_Input2_S154_L003_sorted_downsampled.bam"
)

for bam in "${downsampled_bams[@]}"; do
  [[ -f "${bam}" ]] || {
    echo "Missing downsampled sorted bam: ${bam}" >&2
    exit 1
  }
done

for bam in "${downsampled_bams[@]}"; do
  sample="$(basename "${bam%.bam}")"
  bamCoverage \
    -b "${bam}" \
    -o "${outdir}/tracks/${sample}.bw" \
    -of bigwig \
    -p 8
done

call_peaks() {
  local chip_bam="$1"
  local input_bam="$2"
  local sample="$3"

  mkdir -p "${outdir}/peaks/${sample}"

  macs3 callpeak \
    -t "${chip_bam}" \
    -c "${input_bam}" \
    -f BAMPE \
    -g hs \
    -n "${sample}" \
    --outdir "${outdir}/peaks/${sample}" \
    --broad \
    --broad-cutoff 0.1
}

call_peaks \
  "${downsampled_bams[0]}" \
  "${downsampled_bams[1]}" \
  "Scramble_H4K16ac_rep1"

call_peaks \
  "${downsampled_bams[2]}" \
  "${downsampled_bams[3]}" \
  "Scramble_H4K16ac_rep2"

call_peaks \
  "${downsampled_bams[4]}" \
  "${downsampled_bams[5]}" \
  "YEATS2KD_H4K16ac_rep1"

call_peaks \
  "${downsampled_bams[6]}" \
  "${downsampled_bams[7]}" \
  "YEATS2KD_H4K16ac_rep2"

Rscript scripts/peak_annotation/chipseekerannotation.r \
  "${outdir}/peaks/Scramble_H4K16ac_rep1/Scramble_H4K16ac_rep1_peaks.broadPeak" \
  "${outdir}/peaks/Scramble_H4K16ac_rep2/Scramble_H4K16ac_rep2_peaks.broadPeak" \
  "${outdir}/peaks/YEATS2KD_H4K16ac_rep1/YEATS2KD_H4K16ac_rep1_peaks.broadPeak" \
  "${outdir}/peaks/YEATS2KD_H4K16ac_rep2/YEATS2KD_H4K16ac_rep2_peaks.broadPeak"
