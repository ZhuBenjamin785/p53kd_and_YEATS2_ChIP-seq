#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 05:00:00
#SBATCH --mem=25G
#SBATCH -n 2
#SBATCH -N 1

#SBATCH --cpus-per-task=16

set -euo pipefail
cd /gpfs/projects/b1042/LauberthLab/BenFolder

# Legacy entry point for the shP53 single-end RNA-seq counting workflow.
bash shP53_RNAseq/featurecounts/postalignment_rnaseq.sh
