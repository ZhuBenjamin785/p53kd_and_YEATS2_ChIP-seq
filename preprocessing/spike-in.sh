#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 2:00:00
#SBATCH --mem=10G
#SBATCH -n 2
#SBATCH -N 1

#SBATCH --cpus-per-task=8

set -euo pipefail

cd /projects/b1042/LauberthLab/BenFolder || exit 1

module load anaconda3
source /software/anaconda3/2018.12/etc/profile.d/conda.sh
conda activate rseqc_env
module load MACS3/3.0.2
module load deeptools/3.5.6
module load homer/5.1
module load subread/2.0.3
module load samtools

threads="${SLURM_CPUS_PER_TASK:-16}"
inputdir="${inputdir:-BAMfiles/dm6}"
outdir="${outdir:-spikein_results}"

mkdir -p \
  "${outdir}/filtered_bam" \
  "${outdir}/qc"

printf "Sample\tdm6_fragments\n" > "${outdir}/dm6_fragment_counts.tsv"

for bam in "${inputdir}"/*.bam; do
  [[ -e "${bam}" ]] || continue

  sample="$(basename "${bam%.bam}")"
  filtered_bam="${outdir}/filtered_bam/${sample}.filtered.bam"
  filtered_flagstat="${outdir}/qc/${sample}.filtered.flagstat.txt"
  filtered_idxstats="${outdir}/qc/${sample}.filtered.idxstats.txt"

  echo "Processing ${sample}"

  samtools view -@ "${threads}" -b -F 0x404 "${bam}" \
    | samtools sort -@ "${threads}" -o "${filtered_bam}"
  samtools index "${filtered_bam}"

  samtools flagstat "${filtered_bam}" > "${filtered_flagstat}"
  samtools idxstats "${filtered_bam}" > "${filtered_idxstats}"

  count="$(
    samtools view \
      -c \
      -f 0x42 \
      -F 0xF0C \
      "${filtered_bam}"
  )"

  printf "%s\t%s\n" "${sample}" "${count}" \
    >> "${outdir}/dm6_fragment_counts.tsv"
done

column -t -s $'\t' "${outdir}/dm6_fragment_counts.tsv"
