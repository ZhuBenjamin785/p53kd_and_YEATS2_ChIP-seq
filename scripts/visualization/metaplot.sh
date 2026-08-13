#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 05:00:00
#SBATCH --mem=25G
#SBATCH -N 1
#SBATCH --cpus-per-task=16

set -euo pipefail

cd /projects/b1042/LauberthLab/BenFolder || exit 1

module load deeptools/3.5.6

threads="${SLURM_CPUS_PER_TASK:-16}"
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-${SLURM_JOB_ID:-$$}"
mkdir -p "$MPLCONFIGDIR"

gtf="gencode.v50.basic.annotation.gtf"
genes_bed="protein_coding_genes.bed"
matrix="H4K16ac_protein_coding_gene_body_matrix.gz"
signal_files=(
    "macs3_results_p53kd/tracks/Scr_H4K16ac_2_S0_L001.spikein_normalized.bw"
    "macs3_results_p53kd/tracks/Scr_H4K16ac_1_S0_L001.spikein_normalized.bw"
    "macs3_results_p53kd/tracks/P53_H4K16ac_2_S0_L001.spikein_normalized.bw"
    "macs3_results_p53kd/tracks/P53_H4K16ac_1_S0_L001.spikein_normalized.bw"
)

for input_file in "$gtf" "${signal_files[@]}"; do
    if [[ ! -s "$input_file" ]]; then
        echo "Missing or empty input file: $input_file" >&2
        exit 1
    fi
done

for program in computeMatrix plotProfile plotHeatmap; do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "Required deepTools command not found: $program" >&2
        exit 1
    fi
done

awk 'BEGIN { OFS="\t" }
$3 == "gene" && /gene_type "protein_coding"/ {
    if (match($0, /gene_name "[^"]+"/)) {
        name = substr($0, RSTART + 11, RLENGTH - 12)
        print $1, $4 - 1, $5, name, 0, $7
    }
}' "$gtf" > "$genes_bed"

if [[ ! -s "$genes_bed" ]]; then
    echo "No protein-coding genes were written to $genes_bed" >&2
    exit 1
fi

computeMatrix scale-regions \
    -R "$genes_bed" \
    -S "${signal_files[@]}" \
    --beforeRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --afterRegionStartLength 3000 \
    --binSize 50 \
    --missingDataAsZero \
    --skipZeros \
    --numberOfProcessors "$threads" \
    --samplesLabel \
        "Scramble rep1" \
        "Scramble rep2" \
        "p53KD rep1" \
        "p53KD rep2" \
    -o "$matrix"


plotProfile \
    -m "$matrix" \
    -out H4K16ac_metaplot.pdf

plotHeatmap \
    -m "$matrix" \
    -out H4K16ac_heatmap.pdf

for output_file in "$matrix" H4K16ac_metaplot.pdf H4K16ac_heatmap.pdf; do
    if [[ ! -s "$output_file" ]]; then
        echo "Expected output was not created: $output_file" >&2
        exit 1
    fi
done

echo "Wrote matrix, profile, and heatmap outputs successfully."
