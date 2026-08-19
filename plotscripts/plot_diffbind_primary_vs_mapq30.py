"""Compare primary and MAPQ30 DiffBind tables when both are available."""
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from plot_common import BASE,save,note
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--mapq30",default=""); a=ap.parse_args(); p=BASE/"chipseq/diffbind_broad_remove_duplicates/DiffBind_all_peaks.csv"; q=a.mapq30
    if not q: note("DiffBind comparison","skipped; supply --mapq30 path to corrected table"); return
    x=pd.read_csv(p); y=pd.read_csv(q); n=min(len(x),len(y)); fig,ax=plt.subplots(1,2,figsize=(13,5)); sns.scatterplot(x=x.Fold.iloc[:n],y=y.Fold.iloc[:n],ax=ax[0]); ax[0].axline((0,0),slope=1,color="black",ls="--"); ax[0].set(xlabel="Primary Fold",ylabel="MAPQ30 Fold",title="Peak Fold comparison"); d=pd.DataFrame({"primary":(x.FDR.iloc[:n]<.05).sum(),"MAPQ30":(y.FDR.iloc[:n]<.05).sum()},index=["significant"]); sns.barplot(data=d,ax=ax[1]); ax[1].set_title("Significant peak count"); save(fig,"primary_vs_mapq30_diffbind")
if __name__=="__main__": main()
