set -euo pipefail

project_root="${1:-.}"
dataset="${2:-all}"
package_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${project_root}"

conda run -n chipseq_chipseeker Rscript --vanilla "${package_dir}/chipseq_summary_plots.r" "${dataset}"
conda run -n chipseq_pybw python "${package_dir}/chipseq_summary_plots.py" "${dataset}"

if [[ "${dataset}" == "all" ]]; then
  conda run -n chipseq_pybw python "${package_dir}/chipseq_compare_p53_YEATS2.py"
fi
