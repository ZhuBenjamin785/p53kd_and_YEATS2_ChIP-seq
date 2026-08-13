#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=fastqchip_peaks
#SBATCH -t 24:00:00
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --output=log/slurm-%j.out

set -euo pipefail

project_dir="${FASTQCHIP_PROJECT_DIR:-/gpfs/projects/b1042/LauberthLab/BenFolder}"
input_dir="${FASTQCHIP_BAM_DIR:-${project_dir}/fastqchip_bamfiles/human}"
dm6_input_dir="${FASTQCHIP_DM6_BAM_DIR:-${project_dir}/fastqchip_bamfiles/dm6}"
manifest="${FASTQCHIP_PEAK_MANIFEST:-${project_dir}/fastqchip_peak_manifest.tsv}"
outdir="${FASTQCHIP_PEAK_OUTDIR:-${project_dir}/fastqchip_macs3_results}"
annotation_script="${project_dir}/scripts/peak_annotation/chipseeker_fasterqchip.r"
threads="${SLURM_CPUS_PER_TASK:-16}"
track_threads=$(( threads < 8 ? threads : 8 ))

[[ -d "${project_dir}" ]] || { echo "Project directory not found: ${project_dir}" >&2; exit 1; }
[[ -d "${input_dir}" ]] || { echo "Aligned bam directory not found: ${input_dir}" >&2; exit 1; }
[[ -s "${manifest}" ]] || { echo "Peak manifest is missing or empty: ${manifest}" >&2; exit 1; }
[[ -s "${annotation_script}" ]] || { echo "ChIPseeker script is missing or empty: ${annotation_script}" >&2; exit 1; }
cd "${project_dir}"

module load anaconda3
source /software/anaconda3/2018.12/etc/profile.d/conda.sh
conda activate chipseeker
module load MACS3/3.0.2
module load samtools/1.14
module load deeptools/3.5.6

for program in samtools bamCoverage macs3 Rscript; do
  command -v "${program}" >/dev/null 2>&1 || { echo "Required program is unavailable: ${program}" >&2; exit 1; }
done

mkdir -p "${outdir}/merged_bams" "${outdir}/normalized_bams" "${outdir}/tracks" "${outdir}/peaks" "${outdir}/chipseeker"

declare -A merged_bams
declare -A roles
declare -A input_for
declare -A peak_modes
declare -A dm6_bams
chip_samples=()
target_samples=(input_0hr p53_0hr)
declare -A target_sample_set
for target_sample in "${target_samples[@]}"; do
  target_sample_set["${target_sample}"]=1
done

while IFS=$'\t' read -r sample condition target role input_sample peak_mode runs extra; do
  [[ "${sample}" == "sample" || -z "${sample}" || "${sample:0:1}" == "#" ]] && continue
  [[ -n "${target_sample_set[${sample}]:-}" ]] || continue
  [[ -z "${role}" || -z "${runs}" || -n "${extra:-}" ]] && {
    echo "Malformed manifest row for '${sample:-unknown}'; expected 7 tab-separated columns." >&2
    exit 1
  }
  [[ "${sample}" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo "Unsafe sample name: ${sample}" >&2; exit 1; }
  [[ "${role}" == "input" || "${role}" == "chip" ]] || { echo "role must be input or chip for ${sample}" >&2; exit 1; }
  if [[ "${role}" == "chip" ]]; then
    [[ "${input_sample}" != "." && ( "${peak_mode}" == "broad" || "${peak_mode}" == "narrow" ) ]] || {
      echo "ChIP sample ${sample} needs an input_sample and peak_mode (broad or narrow)." >&2
      exit 1
    }
  fi

  IFS=',' read -r -a run_ids <<< "${runs}"
  source_bams=()
  dm6_source_bams=()
  for run_id in "${run_ids[@]}"; do
    [[ "${run_id}" =~ ^SRR[0-9]+$ ]] || { echo "Invalid run accession '${run_id}' for ${sample}" >&2; exit 1; }
    bam="${input_dir}/${run_id}.sorted.bam"
    [[ -s "${bam}" ]] || { echo "Missing or empty bam for ${sample}: ${bam}" >&2; exit 1; }
    source_bams+=("${bam}")
    dm6_bam="${dm6_input_dir}/${run_id}.sorted.bam"
    [[ -s "${dm6_bam}" ]] || { echo "Missing dm6 bam for ${sample}: ${dm6_bam}" >&2; exit 1; }
    dm6_source_bams+=("${dm6_bam}")
  done

  merged_bam="${outdir}/merged_bams/${sample}.bam"
  samtools merge -f -@ "${track_threads}" "${merged_bam}" "${source_bams[@]}"
  samtools index -@ "${track_threads}" "${merged_bam}"
  dm6_merged_bam="${outdir}/merged_bams/${sample}.dm6.bam"
  samtools merge -f -@ "${track_threads}" "${dm6_merged_bam}" "${dm6_source_bams[@]}"
  samtools index -@ "${track_threads}" "${dm6_merged_bam}"

  merged_bams["${sample}"]="${merged_bam}"
  roles["${sample}"]="${role}"
  input_for["${sample}"]="${input_sample}"
  peak_modes["${sample}"]="${peak_mode}"
  dm6_bams["${sample}"]="${dm6_merged_bam}"
  [[ "${role}" == "chip" ]] && chip_samples+=("${sample}")
done < "${manifest}"

(( ${#chip_samples[@]} > 0 )) || { echo "No ChIP rows found in ${manifest}" >&2; exit 1; }

peak_files=()
dm6_counts=()
chip_dm6_counts=()
minimum_dm6_count=""
for sample in "${chip_samples[@]}"; do
  count="$(samtools view -@ "${track_threads}" -c -F 2308 "${dm6_bams[${sample}]}" )"
  (( count > 0 )) || { echo "No mapped dm6 reads for ${sample}" >&2; exit 1; }
  chip_dm6_counts+=("${count}")
  if [[ -z "${minimum_dm6_count}" ]] || (( count < minimum_dm6_count )); then minimum_dm6_count="${count}"; fi
done

factor_file="${outdir}/spikein_normalization_factors.tsv"
printf 'sample\tdm6_mapped_reads\ttarget_reads\tsubsampling_fraction\n' > "${factor_file}"
normalized_bams=()
normalize_one() {
  local sample="$1" source="$2" count="$3" output="$4" fraction
  fraction="$(awk -v target="${minimum_dm6_count}" -v observed="${count}" 'BEGIN { printf "%.12f", target/observed }')"
  printf '%s\t%s\t%s\t%s\n' "${sample}" "${count}" "${minimum_dm6_count}" "${fraction}" >> "${factor_file}"
  samtools view -@ "${track_threads}" -b -s "53${fraction#0}" "${source}" \
    | samtools sort -@ "${track_threads}" -o "${output}"
  samtools index -@ "${track_threads}" "${output}"
}

for sample in "${chip_samples[@]}"; do
  input_sample="${input_for[${sample}]}"
  count="${chip_dm6_counts[0]}"
  for i in "${!chip_samples[@]}"; do [[ "${chip_samples[i]}" == "${sample}" ]] && count="${chip_dm6_counts[i]}"; done
  chip_norm="${outdir}/normalized_bams/${sample}.spikein_normalized.bam"
  input_norm="${outdir}/normalized_bams/${input_sample}.spikein_normalized.bam"
  normalize_one "${sample}" "${merged_bams[${sample}]}" "${count}" "${chip_norm}"
  normalize_one "${input_sample}" "${merged_bams[${input_sample}]}" "${count}" "${input_norm}"
  normalized_bams+=("${chip_norm}" "${input_norm}")
done

for bam in "${normalized_bams[@]}"; do
  sample="$(basename "${bam%.bam}")"
  bamCoverage -b "${bam}" -o "${outdir}/tracks/${sample}.bw" \
    -of bigwig --normalizeUsing CPM --binSize 10 -p "${track_threads}"
done

for sample in "${chip_samples[@]}"; do
  input_sample="${input_for[${sample}]}"
  input_bam="${merged_bams[${input_sample}]:-}"
  [[ -n "${input_bam}" && "${roles[${input_sample}]:-}" == "input" ]] || {
    echo "Matched input '${input_sample}' for ${sample} is absent or is not role=input." >&2
    exit 1
  }
  peak_dir="${outdir}/peaks/${sample}"
  mkdir -p "${peak_dir}"
  if [[ "${peak_modes[${sample}]}" == "broad" ]]; then
  macs3 callpeak -t "${outdir}/normalized_bams/${sample}.spikein_normalized.bam" -c "${outdir}/normalized_bams/${input_sample}.spikein_normalized.bam" -f bam -g hs \
      -n "${sample}" --outdir "${peak_dir}" --broad --broad-cutoff "${BROAD_CUTOFF:-0.1}"
    peak_files+=("${peak_dir}/${sample}_peaks.broadPeak")
  else
    macs3 callpeak -t "${outdir}/normalized_bams/${sample}.spikein_normalized.bam" -c "${outdir}/normalized_bams/${input_sample}.spikein_normalized.bam" -f bam -g hs \
      -n "${sample}" --outdir "${peak_dir}" -q "${NARROW_QVALUE:-0.05}"
    peak_files+=("${peak_dir}/${sample}_peaks.narrowPeak")
  fi
done

for peak_file in "${peak_files[@]}"; do
  [[ -s "${peak_file}" ]] || { echo "MACS3 did not produce a peak file: ${peak_file}" >&2; exit 1; }
done
Rscript --vanilla "${annotation_script}" --output-dir "${outdir}/chipseeker" "${peak_files[@]}"
