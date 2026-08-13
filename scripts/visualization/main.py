
from pathlib import Path
import pickle
import re
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.decomposition import pca
from sklearn.preprocessing import StandardScaler
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats



projectdir = Path(__file__).resolve().parent
datafile = projectdir / "data" / "dataset.csv"
outputdir = projectdir / "output_files"
outputdir.mkdir(exist_ok=True)

conditions = {
    "21161R-62-01": "TOP1KD",
    "21161R-62-02": "TOP1KD",
    "21161R-62-03": "dTAG_neg",
    "21161R-62-04": "dTAG_neg",
}


def load_counts_and_metadata(datafile: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load gene-by-sample counts and build PyDESeq2 sample metadata."""
    gene_table = pd.read_csv(datafile, sep="\t", index_col="Geneid")

    counts_r = gene_table.drop(columns="Length")

                 
                                                                  

                                     
                                            
                                       
                                 
                                                                 
    
    
    
    minCPM = 0.3
    library_sizes = counts_r.sum(axis=0)
    cpm_data = counts_r.div(library_sizes, axis=1) * 1_000_000
    count_check = cpm_data > minCPM
    keep = count_check.sum(axis=1) >= 2
    print(keep.sum())
    countData = counts_r.loc[keep].copy()
    countData.columns = [re.sub(r"_S.*", "", col) for col in countData.columns]
    
    
    counts_df = countData.T.astype(int)

    missing_samples = set(counts_df.index) - set(conditions)
    if missing_samples:
        raise ValueError(f"No condition assigned to: {sorted(missing_samples)}")

    metadata = pd.DataFrame(
        {"condition": [conditions[sample] for sample in counts_df.index]},
        index=counts_df.index,
    )
    metadata.index.name = "sample"
    return counts_df, metadata

def makepcaplot(counts_df: pd.DataFrame, metadata: pd.DataFrame, outpath: Path) -> None:
    

    

    
    

    
    
    """Make a PCA plot from the sample count matrix."""
                                                       
    log_counts = np.log2(counts_df + 1)

                                                                 
    scaled = StandardScaler().fit_transform(log_counts)

    pca = pca(n_components=2)
    pcs = pca.fit_transform(scaled)

    pca_df = pd.DataFrame(
        pcs,
        index=counts_df.index,
        columns=["PC1", "PC2"],
    ).join(metadata)

    plt.figure(figsize=(9, 6))
    sns.scatterplot(
        data=pca_df,
        x="PC1",
        y="PC2",
        hue="condition",
        s=100,
    )
    for sample, row in pca_df.iterrows():
        plt.text(row["PC1"] + 0.02, row["PC2"] + 0.02, sample, fontsize=8)

    plt.title(
        f"PCA of samples "
        f"(PC1 {pca.explained_variance_ratio_[0] * 100:.1f}%, "
        f"PC2 {pca.explained_variance_ratio_[1] * 100:.1f}%)"
    )
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


def run_analysis() -> None:
    """Create metadata, fit DESeq2, and save the analysis outputs."""
    counts_df, metadata = load_counts_and_metadata(datafile)
    metadata.to_csv(outputdir / "metadata.csv")

    counts_df = counts_df.loc[:, counts_df.sum(axis=0) >= 10]
    
    makepcaplot(counts_df, metadata, outputdir / "pca_plot.png")

    print(f"Loaded {counts_df.shape[0]} samples and {counts_df.shape[1]} genes.")
    print("Metadata written to", outputdir / "metadata.csv")
    print(metadata)

    inference = DefaultInference(n_cpus=1)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design="~condition",
        refit_cooks=True,
        inference=inference,
    )
    
    dds.deseq2()

    with open(outputdir / "result_adata.pkl", "wb") as result_file:
        pickle.dump(dds.to_picklable_anndata(), result_file)

    stats = DeseqStats(
        dds,
        contrast=["condition", "TOP1KD", "dTAG_neg"],
        inference=inference,
    )
    stats.summary()
    stats.results_df.to_csv(outputdir / "results.csv")


if __name__ == "__main__":
    run_analysis()
