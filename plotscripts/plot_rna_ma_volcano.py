"""RNA MA and volcano plots from the canonical DE results."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plot_common import BASE, save
def main():
    f=BASE.parent/"rna_seq_dea/shp53_vs_shLacZ_0hr/results.csv"; df=pd.read_csv(f)
    df["sig"]=(df.padj<=.05)&(df.log2FoldChange.abs()>=1); df["neglog10padj"]=-np.log10(df.padj.clip(lower=1e-300))
    fig,ax=plt.subplots(1,2,figsize=(14,6)); ax[0].scatter(np.log10(df.baseMean.clip(lower=1)),df.log2FoldChange,c=df.sig.map({True:"#D55E00",False:"#BDBDBD"}),s=7,alpha=.6)
    ax[0].axhline(0,color="black",lw=.8); ax[0].set(xlabel="log10(baseMean)",ylabel="log2 fold change",title="RNA MA plot")
    ax[1].scatter(df.log2FoldChange,df.neglog10padj,c=df.sig.map({True:"#D55E00",False:"#BDBDBD"}),s=7,alpha=.6); ax[1].axvline(1,color="black",ls="--"); ax[1].axvline(-1,color="black",ls="--"); ax[1].set(xlabel="log2 fold change",ylabel="-log10(BH padj)",title="RNA volcano plot")
    save(fig,"rna_ma_volcano")
if __name__=="__main__": main()
