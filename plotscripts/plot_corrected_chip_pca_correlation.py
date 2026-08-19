"""Corrected ChIP PCA and correlation heatmap from a region-by-sample matrix."""
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from plot_common import save,note
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--matrix",default=""); a=ap.parse_args()
    if not a.matrix: note("ChIP PCA","skipped; supply --matrix region-by-sample signal matrix"); return
    x=pd.read_csv(a.matrix,index_col=0).select_dtypes("number"); p=PCA(2).fit_transform(np.log1p(x).T); fig,ax=plt.subplots(1,2,figsize=(14,6)); ax[0].scatter(p[:,0],p[:,1],s=100); [ax[0].text(p[i,0],p[i,1],str(s)) for i,s in enumerate(x.columns)]; ax[0].set_title("Corrected ChIP PCA"); sns.heatmap(x.corr(),vmin=-1,vmax=1,cmap="vlag",annot=True,fmt=".2f",ax=ax[1]); ax[1].set_title("Corrected ChIP correlation"); save(fig,"corrected_chip_pca_correlation")
if __name__=="__main__": main()
