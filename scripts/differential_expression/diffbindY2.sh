#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=yeats2_diffbind
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

trap 'status=$?; echo "YEATS2 DiffBind exited with status ${status}" >&2; exit ${status}' EXIT

cd /projects/b1042/LauberthLab/BenFolder || exit 1

source /software/anaconda3/2022.05/etc/profile.d/conda.sh
module load anaconda3/2022.05
module load samtools

conda activate chipseeker

echo "Using R: $(command -v Rscript)"

index_threads="${SLURM_CPUS_PER_TASK:-1}"
required_bams=(
  "BAMfiles/human/H4K16ac_Scrameble_ChIP1_S155_L003_sorted.bam"
  "BAMfiles/human/H4K16ac_Scrameble_ChIP2_S156_L003_sorted.bam"
  "BAMfiles/human/H4K16ac_YEATS2_shRNA#3_ChIP1_S157_L003_sorted.bam"
  "BAMfiles/human/H4K16ac_YEATS2_shRNA#3_ChIP2_S158_L003_sorted.bam"
  "BAMfiles/human/scramble_Input_rep1_S151_L003_sorted.bam"
  "BAMfiles/human/scramble_Input_rep2_S152_L003_sorted.bam"
  "BAMfiles/human/YEATS2_shRNA#3_Input1_S153_L003_sorted.bam"
  "BAMfiles/human/YEATS2_shRNA#3_Input2_S154_L003_sorted.bam"
  "BAMfiles/dm6/H4K16ac_Scrameble_ChIP1_S155_L003_sorted.bam"
  "BAMfiles/dm6/H4K16ac_Scrameble_ChIP2_S156_L003_sorted.bam"
  "BAMfiles/dm6/H4K16ac_YEATS2_shRNA#3_ChIP1_S157_L003_sorted.bam"
  "BAMfiles/dm6/H4K16ac_YEATS2_shRNA#3_ChIP2_S158_L003_sorted.bam"
)

for bam in "${required_bams[@]}"; do
  [[ -s "${bam}" ]] || {
    echo "Missing or empty bam: ${bam}" >&2
    exit 1
  }
  if [[ ! -s "${bam}.bai" ]]; then
    echo "Creating missing bam index: ${bam}.bai"
    samtools index -@ "${index_threads}" "${bam}"
  fi
  samtools quickcheck -v "${bam}"
done

Rscript --vanilla -e 'if (!requireNamespace("DiffBind", quietly=TRUE)) stop("DiffBind is not installed in the active R environment")'
Rscript --vanilla scripts/differential_expression/diffbindY2.r
