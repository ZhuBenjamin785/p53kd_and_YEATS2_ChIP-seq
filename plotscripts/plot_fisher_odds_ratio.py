"""Forest-style plot of Fisher odds ratios and confidence intervals."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plot_common import BASE,save
def main():
    d=pd.read_csv(BASE/"integration/tables/fisher_exact_results.csv"); y=np.arange(len(d)); lo=[]; hi=[]
    for _,r in d.iterrows():
        a,b,c,e=r.chip_and_rna,r.chip_and_not_rna,r.no_chip_and_rna,r.no_chip_and_not_rna; lo.append(np.exp(np.log(max(r.odds_ratio,1e-6))-1.96*np.sqrt(1/max(a,1)+1/max(b,1)+1/max(c,1)+1/max(e,1)))); hi.append(np.exp(np.log(max(r.odds_ratio,1e-6))+1.96*np.sqrt(1/max(a,1)+1/max(b,1)+1/max(c,1)+1/max(e,1))))
    fig,ax=plt.subplots(figsize=(10,4)); ax.errorbar(d.odds_ratio,y,xerr=[d.odds_ratio-lo, np.array(hi)-d.odds_ratio],fmt="o",capsize=4); ax.axvline(1,color="black",ls="--"); ax.set_yticks(y,d.test); ax.set_xscale("log"); ax.set_xlabel("Odds ratio (log scale)"); ax.set_title("Corrected Fisher odds ratios"); save(fig,"fisher_odds_ratio")
if __name__=="__main__": main()
