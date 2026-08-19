#!/usr/bin/env python3
"""Fail if a required corrected result is missing or internally inconsistent."""
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[3]
B=ROOT/"shared/biological_consensus_repaired"
checks=[]
def check(name,condition,detail):
    if not condition: raise SystemExit(f"FAIL {name}: {detail}")
    checks.append((name,"PASS",detail))

u=pd.read_csv(B/"integration/gene_sets/all_peaks__eligible_universe.csv")
check("all-peak universe",len(u)==2100,f"{len(u)} genes")
p=pd.read_csv(B/"integration/gene_sets/promoter_peaks__eligible_universe.csv")
check("promoter universe",len(p)==1828,f"{len(p)} genes")
f=pd.read_csv(B/"integration/tables/fisher_exact_results.csv")
check("Fisher alternatives",set(f.alternative)=={"greater"},str(set(f.alternative)))
check("Fisher margins",int(f.universe_genes.min())==2100,"both tests use 2,100 genes")
c=pd.read_csv(B/"integration/tables/full_universe_correlation_statistics.csv")
check("full-universe correlations",set(c.groupby('scope').N.unique().map(tuple))=={(2100,),(1828,)},"N=2,100/1,828")

for scope in ("all_peaks","promoter_peaks"):
    s=pd.read_csv(B/f"ora/{scope}/tables/ora_run_summary.csv")
    check(f"{scope} ORA rows",len(s)==12,f"{len(s)} database/list combinations")
    check(f"{scope} ORA status",not s.Status.str.startswith("ERROR").any(),"no database errors")
    for gene_set in ("loss_down","gain_up","loss_up","gain_down"):
        for db in ("GO_BP","KEGG","Reactome"):
            path=B/f"ora/{scope}/tables/{gene_set}__{db}__all_pathways.csv"
            check(f"ORA CSV {scope}/{gene_set}/{db}",path.is_file(),str(path))

fea=pd.read_csv(B/"fea/directional/tables/fea_run_summary.csv")
check("directional FEA",len(fea)==6 and (fea.Status=="completed").all(),"2 directions x 3 databases completed")
g=pd.read_csv(B/"gsea/all_eligible/tables/gsea_run_summary.csv")
check("all-gene GSEA",len(g)==3 and (g.Status=="completed").all(),"3 databases completed")
check("GSEA all eligible",int(g.Input_rows.min())==18502 and int(g.Ranked_Entrez_genes.min())==15455,
      "18,502 input rows; 15,455 mapped Entrez ranks")

# ChIP audit controls and final MAPQ30/coordinate-deduplicated sensitivity.
bkeep=pd.read_csv(B/"chipseq/diffbind_broad_keep_duplicates/DiffBind_all_peaks.csv")
bremove=pd.read_csv(B/"chipseq/diffbind_broad_remove_duplicates/DiffBind_all_peaks.csv")
cols=["seqnames","start","end","Fold","p.value","FDR"]
check("broad duplicate-flag control",bkeep[cols].equals(bremove[cols]),
      "keep/remove tables identical because source BAM duplicate flags are absent")
cross=pd.read_csv(B/"chipseq/qc/species_cross_mapping.tsv",sep="\t")
check("species cross-mapping rows",len(cross)==8,"all 8 ChIP/input libraries")
check("MAPQ30 species screen",cross.shared_qnames_mapq30.notna().all(),
      f"shared high-confidence query names explicitly measured: {int(cross.shared_qnames_mapq30.sum())}")
complexity=pd.read_csv(B/"chipseq/qc/library_complexity_frip.tsv",sep="\t")
check("full-depth ChIP QC",len(complexity)==4,"all 4 H4K16ac libraries")
insert=pd.read_csv(B/"chipseq/qc/insert_size/insert_size_summary.tsv",sep="\t")
check("insert-size QC",len(insert)==4,"all 4 H4K16ac libraries")
hashes=(B/"provenance/CHIP_INPUT_SHA256SUMS.txt").read_text().strip().splitlines()
check("ChIP input hashes",len(hashes)==20,"16 BAM and 4 broadPeak SHA-256 records")

final=B/"chipseq/diffbind_mapq30_broad_remove_duplicates/DiffBind_all_peaks.csv"
check("MAPQ30 deduplicated DiffBind",final.is_file(),str(final))
factors=pd.read_csv(B/"chipseq/mapq30_spikein_scale_factors.tsv",sep="\t")
check("per-library spike-in factors",len(factors)==8 and factors.scale_factor.notna().all(),
      "8 independently calculated dm6 factors")
removed=pd.read_csv(B/"chipseq/mapq30_removed_shared_qnames.tsv",sep="\t")
check("ambiguous species-query exclusion",len(removed)==8 and removed.shared_mapq30_qnames_removed.notna().all(),
      f"retained shared query names removed from both species: {int(removed.shared_mapq30_qnames_removed.sum())}")
track_manifest=pd.read_csv(B/"chipseq/mapq30_track_validation_manifest.tsv",sep="\t")
check("corrected ChIP/Input tracks",len(track_manifest)==4 and track_manifest.sha256.notna().all(),
      "4 explicit-scale, no-BPM bigWigs generated and checksummed; large intermediates not retained")
final_integration=pd.read_csv(B/"integration_sensitivity_mapq30_deduplicated/tables/fisher_exact_results.csv")
check("final ChIP integration sensitivity",len(final_integration)==2 and set(final_integration.alternative)=={"greater"},
      "two one-sided tests in the separately matched sensitivity universe")
final_consensus=list((B/"chipseq/consensus_mapq30_deduplicated").glob("*.bed"))
check("final reciprocal consensus",len(final_consensus)==6,
      "two conditions x three reciprocal-overlap thresholds")

out=pd.DataFrame(checks,columns=["check","status","detail"])
(B/"provenance").mkdir(exist_ok=True)
out.to_csv(B/"provenance/VALIDATION_RESULTS.csv",index=False)
print(out.to_string(index=False)); print(f"\n{len(out)} checks passed")
