"""Plot MAPQ30 human/dm6 cross-mapping and spike-in factors."""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from plot_common import BASE, save
def main():
    x=pd.read_csv(BASE/"chipseq/qc/species_cross_mapping.tsv",sep="\t"); x["label"]=x["sample"].str.replace("_S0_L001","",regex=False); f=pd.read_csv(BASE.parent/"summary/metadata/p53KD/spikein_normalization_factors.tsv",sep="\t"); f["label"]=f.peak_sample.str.replace("_H4K16ac_rep"," rep",regex=False)
    fig,ax=plt.subplots(1,2,figsize=(14,6)); z=x.melt(id_vars="label",value_vars=["human_primary_proper_fragments_mapq30","dm6_primary_proper_fragments_mapq30"],var_name="species",value_name="fragments"); sns.barplot(data=z,x="label",y="fragments",hue="species",ax=ax[0]); ax[0].tick_params(axis="x",rotation=35); ax[0].set_yscale("log"); ax[0].set_title("MAPQ30 species counts"); sns.barplot(data=f,x="label",y="subsampling_fraction",ax=ax[1],color="#4C78A8"); ax[1].axhline(1,color="black",ls="--"); ax[1].set_title("Spike-in subsampling fraction"); ax[1].tick_params(axis="x",rotation=35); save(fig,"mapq30_crossmapping_spikein")
if __name__=="__main__": main()
