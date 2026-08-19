#!/usr/bin/env python3
"""Create compact summary plots and tables from corrected FEA and GSEA."""
from __future__ import annotations
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_consensus_repair")
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/"shared/biological_consensus_repaired"
OUT=BASE/"summary"
OUT.mkdir(parents=True,exist_ok=True)

def short(s,n=62):
    return s if len(s)<=n else s[:n-1]+"…"

fea_rows=[]
for p in sorted((BASE/"fea/directional/tables").glob("RNA_*__*__all_pathways.csv")):
    d=pd.read_csv(p); sig=d[pd.to_numeric(d.Adjusted_P_value_FDR,errors="coerce")<.05].copy()
    fea_rows.append(sig)
fea=pd.concat(fea_rows,ignore_index=True) if fea_rows else pd.DataFrame()
fea.to_csv(OUT/"directional_fea_significant_pathways.csv",index=False)
for db in ("GO_BP","KEGG","Reactome"):
    x=fea[fea.Database.eq(db)].sort_values(["Direction","Adjusted_P_value_FDR"]).groupby("Direction").head(10).copy()
    if x.empty: continue
    x=x.sort_values(["Direction","Adjusted_P_value_FDR"],ascending=[True,False])
    labels=[short(v) for v in x.Pathway]
    vals=-np.log10(x.Adjusted_P_value_FDR.clip(lower=np.finfo(float).tiny))
    colors=x.Direction.map({"RNA_up":"#B2182B","RNA_down":"#2166AC"})
    fig,h=plt.subplots(figsize=(9,max(4,.32*len(x)+1.3)))
    h.barh(range(len(x)),vals,color=colors)
    h.set_yticks(range(len(x)),labels); h.set_xlabel("−log10(BH FDR)"); h.set_title(f"Directional RNA FEA — {db}")
    h.axvline(-np.log10(.05),color="black",ls="--",lw=.8)
    h.spines[["top","right"]].set_visible(False); fig.tight_layout()
    fig.savefig(OUT/f"directional_FEA_{db}.pdf"); fig.savefig(OUT/f"directional_FEA_{db}.png",dpi=300); plt.close(fig)

gsea_rows=[]
for p in sorted((BASE/"gsea/all_eligible/tables").glob("*__all_pathways.csv")):
    d=pd.read_csv(p); gsea_rows.append(d[pd.to_numeric(d.Adjusted_P_value_FDR,errors="coerce")<.05].copy())
gsea=pd.concat(gsea_rows,ignore_index=True) if gsea_rows else pd.DataFrame()
gsea.to_csv(OUT/"all_gene_gsea_significant_pathways.csv",index=False)
for db in ("GO_BP","KEGG","Reactome"):
    d=gsea[gsea.Database.eq(db)].copy()
    pos=d[d.NES>0].nsmallest(8,"Adjusted_P_value_FDR")
    neg=d[d.NES<0].nsmallest(8,"Adjusted_P_value_FDR")
    x=pd.concat([neg.sort_values("NES",ascending=False),pos.sort_values("NES")])
    if x.empty: continue
    labels=[short(v) for v in x.Pathway]
    colors=np.where(x.NES>0,"#B2182B","#2166AC")
    fig,h=plt.subplots(figsize=(9,max(4,.34*len(x)+1.3)))
    h.barh(range(len(x)),x.NES,color=colors); h.axvline(0,color="black",lw=.8)
    h.set_yticks(range(len(x)),labels); h.set_xlabel("Normalized enrichment score")
    h.set_title(f"All-eligible-gene GSEA — {db}\nred: higher after p53 knockdown; blue: lower")
    h.spines[["top","right"]].set_visible(False); fig.tight_layout()
    fig.savefig(OUT/f"all_gene_GSEA_{db}.pdf"); fig.savefig(OUT/f"all_gene_GSEA_{db}.png",dpi=300); plt.close(fig)

print(f"FEA significant rows: {len(fea)}; GSEA significant rows: {len(gsea)}; output: {OUT}")
