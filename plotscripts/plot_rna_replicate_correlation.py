"""RNA replicate-correlation heatmap from a sample-by-gene count matrix."""
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from plot_common import BASE, save, note
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--counts", default=""); a=ap.parse_args()
    path=a.counts or str(BASE.parent/"rna_seq_dea/shp53_vs_shLacZ_0hr/counts_matrix.csv")
    try: df=pd.read_csv(path,index_col=0).select_dtypes("number").T
    except FileNotFoundError: note("RNA correlation",f"skipped; count matrix not found: {path}"); return
    if len(df)<2: note("RNA correlation","skipped; need at least two samples"); return
    fig,ax=plt.subplots(figsize=(8,7)); sns.heatmap(df.corr(method="spearman"),vmin=-1,vmax=1,cmap="vlag",annot=True,fmt=".2f",ax=ax)
    ax.set_title("RNA replicate correlation (Spearman)"); save(fig,"rna_replicate_correlation")
if __name__ == "__main__": main()
