#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 3:00:00
#SBATCH --mem=32G
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=chipseq_summary_plots
#SBATCH --output=log/slurm-%j.out
#SBATCH --error=log/slurm-%j.err

set -euo pipefail
cd /gpfs/projects/b1042/LauberthLab/BenFolder

module load anaconda3/2022.05
source /software/anaconda3/2022.05/etc/profile.d/conda.sh
conda activate chipseeker

echo "Using R: $(command -v Rscript)"
Rscript --vanilla -e 'pkgs <- c("ChIPseeker","TxDb.Hsapiens.UCSC.hg38.knownGene","org.Hs.eg.db"); m <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if(length(m)) stop("Missing packages: ", paste(m, collapse=", "))'
Rscript --vanilla scripts/chipseq_summary/chipseq_summary_plots.r

conda activate pybw
python -E scripts/chipseq_summary/chipseq_summary_plots.py
