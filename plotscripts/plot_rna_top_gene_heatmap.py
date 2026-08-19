"""Heatmap of top RNA genes; requires a sample-by-gene count matrix."""
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from plot_common import BASE, save, note
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--counts",default=""); ap.add_argument("--n",type=int,default=30); a=ap.parse_args()
    p=a.counts or str(BASE.parent/"rna_seq_dea/shp53_vs_shLacZ_0hr/counts_matrix.csv")
    try: x=pd.read_csv(p,index_col=0).select_dtypes("number")
    except FileNotFoundError: note("RNA heatmap",f"skipped; count matrix not found: {p}"); return
    res=pd.read_csv(BASE.parent/"rna_seq_dea/shp53_vs_shLacZ_0hr/results.csv").set_index("gene_id"); genes=res.reindex(x.index).sort_values("padj").head(a.n).index; z=np.log1p(x.loc[genes]).sub(np.log1p(x.loc[genes]).mean(1),axis=0).div(np.log1p(x.loc[genes]).std(1).replace(0,1),axis=0)
    fig,ax=plt.subplots(figsize=(10,10)); sns.heatmap(z,cmap="vlag",center=0,yticklabels=res.reindex(genes).gene_name.fillna(genes),ax=ax); ax.set_title("Top RNA differential genes"); save(fig,"rna_top_gene_heatmap")
if __name__=="__main__": main()
