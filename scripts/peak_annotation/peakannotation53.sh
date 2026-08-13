#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=p53kd_peaks
#SBATCH -t 18:00:00
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH -N 1
#SBATCH --cpus-per-task=32

set -euo pipefail

project_dir="${P53KD_PROJECT_DIR:-/projects/b1042/LauberthLab/BenFolder}"
[[ -d "${project_dir}" ]] || {
  echo "Project directory not found: ${project_dir}" >&2
  exit 1
}
cd "${project_dir}"

module load anaconda3
source /software/anaconda3/2018.12/etc/profile.d/conda.sh
conda activate chipseeker
module load MACS3/3.0.2
module load samtools
module load deeptools/3.5.6

human_bam_dir="${human_bam_dir:-p53kdbamfiles/human}"
dm6_bam_dir="${dm6_bam_dir:-p53kdbamfiles/dm6}"
outdir="${outdir:-macs3_results_p53kd}"
samtools_threads=8

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

chip_samples=(
  "Scr_H4K16ac_1_S0_L001"
  "Scr_H4K16ac_2_S0_L001"
  "P53_H4K16ac_1_S0_L001"
  "P53_H4K16ac_2_S0_L001"
)
input_samples=(
  "Scr_Input1_S0_L001"
  "Scr_Input2_S0_L001"
  "MutP53_Input1_S0_L001"
  "MutP53_Input2_S0_L001"
)
peak_samples=(
  "Scramble_H4K16ac_rep1"
  "Scramble_H4K16ac_rep2"
  "p53KD_H4K16ac_rep1"
  "p53KD_H4K16ac_rep2"
)

required_bams=()
for sample in "${chip_samples[@]}" "${input_samples[@]}"; do
  required_bams+=("${human_bam_dir}/${sample}.sorted.bam")
done
for sample in "${chip_samples[@]}"; do
  required_bams+=("${dm6_bam_dir}/${sample}.sorted.bam")
done

for bam in "${required_bams[@]}"; do
  [[ -s "${bam}" ]] || {
    echo "Missing or empty bam: ${bam}" >&2
    echo "Submit this job with an afterok dependency on the alignment job." >&2
    exit 1
  }
  [[ -s "${bam}.bai" ]] || {
    echo "Missing or empty bam index: ${bam}.bai" >&2
    echo "The alignment job may not have completed successfully." >&2
    exit 1
  }
done

if ! samtools quickcheck -v "${required_bams[@]}"; then
  echo "One or more alignment BAMs failed samtools quickcheck." >&2
  exit 1
fi

mkdir -p "${outdir}/peaks" "${outdir}/normalized_bam" "${outdir}/tracks"

dm6_counts=()
minimum_dm6_count=""
for sample in "${chip_samples[@]}"; do
  bam="${dm6_bam_dir}/${sample}.sorted.bam"
  count="$(samtools view -@ "${samtools_threads}" -c -f 2 -F 3844 "${bam}")"
  (( count > 0 )) || {
    echo "No primary, properly paired dm6 alignments found in ${bam}" >&2
    exit 1
  }
  dm6_counts+=("${count}")
  if [[ -z "${minimum_dm6_count}" ]] || (( count < minimum_dm6_count )); then
    minimum_dm6_count="${count}"
  fi
done

normalize_bam() {
  local source_bam="$1"
  local destination_bam="$2"
  local source_dm6_count="$3"
  local fraction

  if (( source_dm6_count == minimum_dm6_count )); then
    samtools view -@ "${samtools_threads}" -bh \
      -o "${destination_bam}" "${source_bam}"
  else
    fraction="$(awk -v target="${minimum_dm6_count}" -v observed="${source_dm6_count}" \
      'BEGIN { printf "%.12f", target / observed }')"
    samtools view -@ "${samtools_threads}" -bh \
      -s "53${fraction#0}" -o "${destination_bam}" "${source_bam}"
  fi
  samtools index -@ "${samtools_threads}" "${destination_bam}"
}

factor_file="${outdir}/spikein_normalization_factors.tsv"
printf 'peak_sample\tdm6_alignment_count\ttarget_count\tsubsampling_fraction\n' > "${factor_file}"

normalized_chip_bams=()
normalized_input_bams=()
all_normalized_bams=()
for index in "${!chip_samples[@]}"; do
  fraction="$(awk -v target="${minimum_dm6_count}" -v observed="${dm6_counts[index]}" \
    'BEGIN { printf "%.12f", target / observed }')"
  printf '%s\t%s\t%s\t%s\n' \
    "${peak_samples[index]}" "${dm6_counts[index]}" \
    "${minimum_dm6_count}" "${fraction}" >> "${factor_file}"

  chip_output="${outdir}/normalized_bam/${chip_samples[index]}.spikein_normalized.bam"
  input_output="${outdir}/normalized_bam/${input_samples[index]}.spikein_normalized.bam"

  normalize_bam \
    "${human_bam_dir}/${chip_samples[index]}.sorted.bam" \
    "${chip_output}" \
    "${dm6_counts[index]}"
  normalize_bam \
    "${human_bam_dir}/${input_samples[index]}.sorted.bam" \
    "${input_output}" \
    "${dm6_counts[index]}"

  normalized_chip_bams+=("${chip_output}")
  normalized_input_bams+=("${input_output}")
  all_normalized_bams+=("${chip_output}" "${input_output}")
done

for bam in "${all_normalized_bams[@]}"; do
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

peak_files=()
for index in "${!peak_samples[@]}"; do
  call_peaks \
    "${normalized_chip_bams[index]}" \
    "${normalized_input_bams[index]}" \
    "${peak_samples[index]}"
  peak_files+=(
    "${outdir}/peaks/${peak_samples[index]}/${peak_samples[index]}_peaks.broadPeak"
  )
done

Rscript "${annotation_script}" "${peak_files[@]}"
