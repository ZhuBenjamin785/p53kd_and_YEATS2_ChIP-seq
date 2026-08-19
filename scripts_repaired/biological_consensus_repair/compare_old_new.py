#!/usr/bin/env python3
"""Create a compact numerical comparison of superseded and repaired outputs."""
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[3]
NEW=ROOT/"shared/biological_consensus_repaired"
OLD=ROOT/"shared/rna_chip_integration/p53KD_H4K16ac_vs_RNAseq"
rows=[]
def add(section,metric,old,new,note=""):
    rows.append(dict(section=section,metric=metric,old_value=old,new_value=new,note=note))

of=pd.read_csv(OLD/"fisher_exact_results.csv"); nf=pd.read_csv(NEW/"integration/tables/fisher_exact_results.csv")
for test in nf.test:
    a=of[of.test.eq(test)].iloc[0]; b=nf[nf.test.eq(test)].iloc[0]
    for metric,oc,nc in (("universe_genes","universe_genes","universe_genes"),("overlap","overlap","overlap"),
                         ("odds_ratio","odds_ratio","odds_ratio"),("one_sided_p","fisher_p_value_greater","fisher_p_value")):
        add("Fisher",f"{test}: {metric}",a[oc],b[nc],"alternative=greater")

oc=pd.read_csv(OLD/"summary/all_correlation_statistics.csv")
nc=pd.read_csv(NEW/"integration/tables/full_universe_correlation_statistics.csv")
for scope in ("all_peaks","promoter_peaks"):
    for method in ("Pearson","Spearman"):
        a=oc[(oc.Analysis.eq(scope))&(oc.Method.eq(method))].iloc[0]
        b=nc[(nc.scope.eq(scope))&(nc.chip_gene_summary.eq("median"))&(nc.method.eq(method))].iloc[0]
        add("Correlation",f"{scope} {method}: N",a.N,b.N,"old: significant-only; new: full eligible universe")
        add("Correlation",f"{scope} {method}: coefficient",a.Correlation,b.correlation,
            "new primary ChIP summary=median Fold across assigned tested peaks")
        add("Correlation",f"{scope} {method}: p_value",a.P_value,b.p_value)

oc=pd.read_csv(OLD/"summary/all_category_counts.csv")
nc=pd.read_csv(NEW/"integration/tables/significant_category_counts.csv")
label={"H4K16ac loss + RNA down":"loss_down","H4K16ac gain + RNA up":"gain_up",
       "H4K16ac loss + RNA up":"loss_up","H4K16ac gain + RNA down":"gain_down"}
for _,a in oc.iterrows():
    b=nc[(nc.scope.eq(a.Analysis))&(nc.category.eq(label[a.Biological_category]))].iloc[0]
    add("Integrated categories",f"{a.Analysis}: {label[a.Biological_category]}",a.Gene_count,b['count'])

old_ora=OLD/"ora_matched_universes"
for scope in ("all_peaks","promoter_peaks"):
    op=old_ora/scope/"tables/ora_run_summary.csv"; np=NEW/"ora"/scope/"tables/ora_run_summary.csv"
    if op.exists() and np.exists():
        a=pd.read_csv(op); b=pd.read_csv(np)
        for _,x in b.iterrows():
            hit=a[(a.Gene_set.eq(x.Gene_set))&(a.Database.eq(x.Database))]
            old=hit.iloc[0].Significant_pathways if len(hit) else "NA"
            add("Integrated ORA",f"{scope} {x.Gene_set} {x.Database}: significant pathways",old,x.Significant_pathways)

old_g=pd.read_csv(ROOT/"shared/rna_seq_dea/shp53_vs_shLacZ_0hr/gsea_out/gseapy.gene_set.prerank.report.csv")
new_g=pd.read_csv(NEW/"gsea/all_eligible/tables/gsea_run_summary.csv")
add("GSEA","ranked input genes",1910,new_g.Ranked_Entrez_genes.iloc[0],
    "old input was significance-filtered symbols; new input is all finite Wald-statistic genes mapped to Entrez")
add("GSEA","KEGG significant pathways",int((old_g['FDR q-val']<.05).sum()),
    int(new_g.loc[new_g.Database.eq('KEGG'),'Significant_pathways'].iloc[0]),
    "collections/versions differ, so pathway counts are descriptive rather than a controlled method-only comparison")

# ChIP region-definition and duplicate-policy controls.  The source BAMs have
# no duplicate flags, so DiffBind's bRemoveDuplicates switch alone cannot
# remove coordinate duplicates.  Identical keep/remove results document that
# the observed change is caused by broad-region handling, not that duplicates
# have been resolved.
chip_files={
    "fixed_summit_keep_duplicates": ROOT/"shared/summary/tables/diffbind/DiffBind_all_peaks.csv",
    "broad_keep_duplicates": NEW/"chipseq/diffbind_broad_keep_duplicates/DiffBind_all_peaks.csv",
    "broad_remove_flagged_duplicates": NEW/"chipseq/diffbind_broad_remove_duplicates/DiffBind_all_peaks.csv",
    "mapq30_coordinate_deduplicated_broad": NEW/"chipseq/diffbind_mapq30_broad_remove_duplicates/DiffBind_all_peaks.csv",
}
for label,path in chip_files.items():
    if not path.exists():
        continue
    d=pd.read_csv(path); sig=d[d.FDR.lt(.05)]
    width=d['width'].median() if 'width' in d else (d.End-d.Start+1).median()
    add("ChIP sensitivity",f"{label}: tested regions","NA",len(d))
    add("ChIP sensitivity",f"{label}: significant loss regions","NA",int(sig.Fold.lt(0).sum()))
    add("ChIP sensitivity",f"{label}: significant gain regions","NA",int(sig.Fold.gt(0).sum()))
    add("ChIP sensitivity",f"{label}: median region width","NA",width,
        "bp; broad keep/remove equality reflects absent duplicate flags")

final_fisher=NEW/"integration_sensitivity_mapq30_deduplicated/tables/fisher_exact_results.csv"
if final_fisher.exists():
    primary=pd.read_csv(NEW/"integration/tables/fisher_exact_results.csv")
    sensitivity=pd.read_csv(final_fisher)
    for test in sensitivity.test:
        a=primary[primary.test.eq(test)].iloc[0]; b=sensitivity[sensitivity.test.eq(test)].iloc[0]
        for metric in ("universe_genes","overlap","odds_ratio","fisher_p_value"):
            old_col="fisher_p_value" if metric=="fisher_p_value" else metric
            add("Final ChIP sensitivity",f"{test}: {metric}",a[old_col],b[metric],
                "primary required 2,100-gene universe vs separately matched MAPQ30 coordinate-deduplicated universe")

out=pd.DataFrame(rows)
(NEW/"summary").mkdir(parents=True,exist_ok=True)
out.to_csv(NEW/"summary/old_vs_new_numeric_comparison.csv",index=False)
print(out.to_string(index=False))
