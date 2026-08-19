"""ChIP replicate-correlation heatmap from the repaired DiffBind table."""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from plot_common import BASE, save
def main():
    f=BASE/"chipseq/diffbind_broad_remove_duplicates/replicate_correlations.csv"; d=pd.read_csv(f); names=sorted(set(d.sample_1)|set(d.sample_2)); m=pd.DataFrame(index=names,columns=names,dtype=float); m.loc[:,:]=1
    for _,r in d.iterrows(): m.loc[r.sample_1,r.sample_2]=m.loc[r.sample_2,r.sample_1]=r.pearson_log2_score
    fig,ax=plt.subplots(figsize=(8,7)); sns.heatmap(m.astype(float),vmin=-1,vmax=1,cmap="vlag",annot=True,fmt=".2f",ax=ax); ax.set_title("ChIP replicate Pearson correlation"); save(fig,"chip_replicate_correlation")
if __name__=="__main__": main()
