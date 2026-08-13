#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=diffbind_p53_0hr
#SBATCH --output=log/slurm-%j.out
#SBATCH --error=log/slurm-%j.err

set -euo pipefail

cd /gpfs/projects/b1042/LauberthLab/BenFolder

module load anaconda3/2022.05
source /software/anaconda3/2022.05/etc/profile.d/conda.sh
module load samtools
conda activate chipseeker

echo "Using R: $(command -v Rscript)"
Rscript --vanilla -e 'if (!requireNamespace("DiffBind", quietly=TRUE)) stop("DiffBind is not installed in the active R environment")'

required_bams=(
  fastqchip_bamfiles/human/SRR5944063.sorted.bam
  fastqchip_bamfiles/human/SRR5944064.sorted.bam
  fastqchip_bamfiles/human/SRR5944081.sorted.bam
  fastqchip_bamfiles/human/SRR5944082.sorted.bam
  fastqchip_bamfiles/dm6/SRR5944063.sorted.bam
  fastqchip_bamfiles/dm6/SRR5944064.sorted.bam
)

for bam in "${required_bams[@]}"; do
  [[ -s "${bam}" ]] || { echo "Missing or empty bam: ${bam}" >&2; exit 1; }
  [[ -s "${bam}.bai" ]] || {
    echo "Creating missing bam index: ${bam}.bai"
    samtools index -@ "${SLURM_CPUS_PER_TASK:-1}" "${bam}"
  }
  samtools quickcheck -v "${bam}"
done

exec Rscript --vanilla scripts/differential_expression/diffbind_fastqchip_p53_0hr.r
