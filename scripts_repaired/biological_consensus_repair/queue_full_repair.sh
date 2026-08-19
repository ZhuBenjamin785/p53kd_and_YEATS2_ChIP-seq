#!/bin/bash
# Submit the complete repaired biological-consensus workflow as one dependency
# graph. This script only queues jobs; SLURM performs the work overnight.
set -euo pipefail
ROOT="/gpfs/projects/b1042/LauberthLab/BenFolder"
S="$ROOT/shared/scripts/biological_consensus_repair"
cd "$ROOT"
mkdir -p shared/biological_consensus_repaired/logs

submit() { sbatch --parsable "$@"; }

# Primary RNA/ChIP integration must finish before ORA/FEA/GSEA.
j_integration=$(submit "$S/primary_integration.slurm")
j_enrichment=$(submit --dependency="afterok:$j_integration" "$S/run_corrected_enrichment.slurm")

# Independent ChIP QC and reproducibility jobs can run concurrently.
j_qc=$(submit "$S/chip_qc_repair.slurm")
j_depth=$(submit "$S/chip_complexity_downsampled.slurm")
j_insert=$(submit "$S/chip_insert_size_qc.slurm")
j_hash=$(submit "$S/checksum_chip_inputs.slurm")
j_consensus=$(submit "$ROOT/p53kdH4K16ac/preprocessing/consensus.sh")
j_broad_remove=$(submit "$S/diffbind_broad_sensitivity.slurm")
j_broad_keep=$(submit "$S/diffbind_broad_keep_duplicates.slurm")

# The corrected MAPQ30/deduplicated sensitivity waits for the complete
# eight-library species screen. It explicitly removes retained shared qnames.
j_chip_final=$(submit --dependency="afterok:$j_qc" "$S/chip_mapq30_reanalysis.slurm")

# Compact final integration, consensus sensitivity, comparisons, and validation
# wait for every result consumed by their checks.
deps="$j_chip_final:$j_depth:$j_insert:$j_hash:$j_consensus:$j_broad_remove:$j_broad_keep:$j_enrichment"
j_post=$(submit --dependency="afterok:$deps" "$S/postprocess_final_chip.slurm")

printf 'stage\tjob_id\n'
printf 'primary_integration\t%s\n' "$j_integration"
printf 'ORA_FEA_GSEA_plots\t%s\n' "$j_enrichment"
printf 'chip_full_QC\t%s\n' "$j_qc"
printf 'chip_depth_matched_QC\t%s\n' "$j_depth"
printf 'chip_insert_QC\t%s\n' "$j_insert"
printf 'chip_input_hashes\t%s\n' "$j_hash"
printf 'legacy_peak_consensus_repair\t%s\n' "$j_consensus"
printf 'broad_remove_flag_control\t%s\n' "$j_broad_remove"
printf 'broad_keep_flag_control\t%s\n' "$j_broad_keep"
printf 'MAPQ30_coordinate_deduplicated_sensitivity\t%s\n' "$j_chip_final"
printf 'postprocess_and_validate\t%s\n' "$j_post"
printf '\nFinal completion check: sacct -j %s --format=JobID,State,Elapsed,ExitCode\n' "$j_post"
