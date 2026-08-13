#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=yeats2_peaks
#SBATCH -t 18:00:00
#SBATCH --mem=64G
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task=32

set -euo pipefail

project_dir="${YEATS2_PROJECT_DIR:-/projects/b1042/LauberthLab/BenFolder}"
[[ -d "${project_dir}" ]] || {
  echo "Project directory not found: ${project_dir}" >&2
  exit 1
}
cd "${project_dir}" || exit 1

module load anaconda3
source /software/anaconda3/2018.12/etc/profile.d/conda.sh
conda activate chipseeker
module load MACS3/3.0.2
module load samtools
module load deeptools/3.5.6

inputdir="${inputdir:-BAMfiles/human}"
outdir="${outdir:-macs3_results_yeats2}"

for program in samtools bamCoverage macs3 Rscript; do
  command -v "${program}" >/dev/null 2>&1 || {
    echo "Required program is unavailable: ${program}" >&2
    exit 1
  }
done

annotation_script="${project_dir}/scripts/peak_annotation/chipseekerannotation.r"
[[ -s "${annotation_script}" ]] || {
  echo "Missing or empty annotation script: ${annotation_script}" >&2
  exit 1
}

mkdir -p "${outdir}/peaks" "${outdir}/dedup_bam" "${outdir}/tracks"

H4K16ac_Scrameble_ChIP1_scale="${H4K16ac_Scrameble_ChIP1_scale:-785}"
H4K16ac_Scrameble_ChIP2_scale="${H4K16ac_Scrameble_ChIP2_scale:-912}"
H4K16ac_YEATS2_ChIP2_scale="${H4K16ac_YEATS2_ChIP2_scale:-685}"

source_bams=(
  "${inputdir}/H4K16ac_Scrameble_ChIP1_S155_L003_sorted.bam"
  "${inputdir}/H4K16ac_Scrameble_ChIP2_S156_L003_sorted.bam"
  "${inputdir}/H4K16ac_YEATS2_shRNA#3_ChIP1_S157_L003_sorted.bam"
  "${inputdir}/H4K16ac_YEATS2_shRNA#3_ChIP2_S158_L003_sorted.bam"
  "${inputdir}/scramble_Input_rep1_S151_L003_sorted.bam"
  "${inputdir}/scramble_Input_rep2_S152_L003_sorted.bam"
  "${inputdir}/YEATS2_shRNA#3_Input1_S153_L003_sorted.bam"
  "${inputdir}/YEATS2_shRNA#3_Input2_S154_L003_sorted.bam"
)

for bam in "${source_bams[@]}"; do
  [[ -s "${bam}" ]] || {
    echo "Missing or empty source bam: ${bam}" >&2
    exit 1
  }
done

for scale in \
  "${H4K16ac_Scrameble_ChIP1_scale}" \
  "${H4K16ac_Scrameble_ChIP2_scale}" \
  "${H4K16ac_YEATS2_ChIP2_scale}"; do
  [[ "${scale}" =~ ^[0-9]{1,3}$ ]] || {
    echo "Invalid subsampling fraction '${scale}'; use 1-3 digits after the decimal point" >&2
    exit 1
  }
done


samtools view -@ 8 -s "34.${H4K16ac_Scrameble_ChIP1_scale}" -b \
  "${inputdir}/H4K16ac_Scrameble_ChIP1_S155_L003_sorted.bam" \
  | samtools sort -@ 8 -o "${outdir}/dedup_bam/H4K16ac_Scrameble_ChIP1_S155_L003_sorted_downsampled.bam"
samtools index "${outdir}/dedup_bam/H4K16ac_Scrameble_ChIP1_S155_L003_sorted_downsampled.bam"

samtools view -@ 8 -s "34.${H4K16ac_Scrameble_ChIP2_scale}" -b \
  "${inputdir}/H4K16ac_Scrameble_ChIP2_S156_L003_sorted.bam" \
  | samtools sort -@ 8 -o "${outdir}/dedup_bam/H4K16ac_Scrameble_ChIP2_S156_L003_sorted_downsampled.bam"
samtools index "${outdir}/dedup_bam/H4K16ac_Scrameble_ChIP2_S156_L003_sorted_downsampled.bam"

samtools sort -@ 8 \
  -o "${outdir}/dedup_bam/H4K16ac_YEATS2_shRNA#3_ChIP1_S157_L003_sorted_downsampled.bam" \
  "${inputdir}/H4K16ac_YEATS2_shRNA#3_ChIP1_S157_L003_sorted.bam"
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

samtools sort -@ 8 \
  -o "${outdir}/dedup_bam/YEATS2_shRNA#3_Input1_S153_L003_sorted_downsampled.bam" \
  "${inputdir}/YEATS2_shRNA#3_Input1_S153_L003_sorted.bam"
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
  [[ -s "${bam}" ]] || {
    echo "Missing or empty normalized bam: ${bam}" >&2
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

Rscript "${annotation_script}" \
  "${outdir}/peaks/Scramble_H4K16ac_rep1/Scramble_H4K16ac_rep1_peaks.broadPeak" \
  "${outdir}/peaks/Scramble_H4K16ac_rep2/Scramble_H4K16ac_rep2_peaks.broadPeak" \
  "${outdir}/peaks/YEATS2KD_H4K16ac_rep1/YEATS2KD_H4K16ac_rep1_peaks.broadPeak" \
  "${outdir}/peaks/YEATS2KD_H4K16ac_rep2/YEATS2KD_H4K16ac_rep2_peaks.broadPeak"
