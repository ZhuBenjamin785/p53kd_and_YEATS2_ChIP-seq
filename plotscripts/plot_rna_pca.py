"""RNA PCA from a sample-by-gene count matrix, if one is available."""
import argparse
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from plot_common import BASE, save, note

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--counts", default="")
    a = ap.parse_args(); path = a.counts or str(BASE.parent / "rna_seq_dea/shp53_vs_shLacZ_0hr/counts_matrix.csv")
    try: df = pd.read_csv(path, index_col=0)
    except FileNotFoundError:
        note("RNA PCA", f"skipped; count matrix not found: {path}"); return
    x = df.select_dtypes("number").T
    if x.shape[0] < 2: note("RNA PCA", "skipped; need at least two samples"); return
    z = StandardScaler().fit_transform(np.log1p(x))
    p = PCA(n_components=2).fit_transform(z)
    fig, ax = plt.subplots(figsize=(7, 6)); ax.scatter(p[:,0], p[:,1], s=100)
    for i, label in enumerate(x.index): ax.text(p[i,0], p[i,1], str(label), fontsize=9)
    ax.set(xlabel="PC1", ylabel="PC2", title="RNA-seq PCA (log1p standardized counts)")
    save(fig, "rna_pca")
if __name__ == "__main__": main()
