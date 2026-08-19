"""Supplementary ChIP FRiP, complexity, and insert-size plots."""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from plot_common import BASE, save
def main():
    qc=BASE/"chipseq/qc"; fr=pd.read_csv(qc/"library_complexity_frip.tsv",sep="\t"); dep=pd.read_csv(qc/"library_complexity_10pct.tsv",sep="\t"); ins=pd.read_csv(qc/"insert_size/insert_size_summary.tsv",sep="\t"); d=fr.merge(dep,on="sample",suffixes=("","_depth10" )).merge(ins,on="sample"); d["label"]=d["sample"].str.replace("_S0_L001","",regex=False)
    fig,ax=plt.subplots(2,2,figsize=(14,10)); sns.barplot(data=d,x="label",y="FRiP_read1_fraction",ax=ax[0,0]); ax[0,0].set_title("Full-depth FRiP"); sns.barplot(data=d,x="label",y="NRF",ax=ax[0,1]); sns.barplot(data=d,x="label",y="PBC1",ax=ax[1,0]); sns.barplot(data=d,x="label",y="average_insert_size",ax=ax[1,1]); ax[1,1].set_title("Mean insert size");
    for a in ax.flat: a.tick_params(axis="x",rotation=35)
    fig.suptitle("ChIP-seq supplementary QC"); save(fig,"chip_frip_complexity_insert_size")
if __name__=="__main__": main()
