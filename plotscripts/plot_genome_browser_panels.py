"""Make lightweight genome-browser-style panels from BED intervals and optional bedGraph tracks.

Use --locus name:chrom:start:end repeatedly. BigWig rendering is intentionally not implicit;
convert selected tracks to bedGraph beforehand or provide --bedgraph files.
"""
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from plot_common import BASE,save,note
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--locus",action="append",default=[]); ap.add_argument("--bed",action="append",default=[]); ap.add_argument("--bedgraph",action="append",default=[]); a=ap.parse_args()
    if not a.locus: note("browser panels","skipped; provide --locus name:chr:start:end"); return
    for spec in a.locus:
        name,chrom,start,end=spec.split(":"); start=int(start); end=int(end); fig,ax=plt.subplots(figsize=(14,3)); ax.set_xlim(start,end); ax.set_yticks([]); ax.set_title(f"{name} ({chrom}:{start:,}-{end:,})")
        for j,b in enumerate(a.bed):
            try: d=pd.read_csv(b,sep="\t",header=None,comment="#")
            except Exception: continue
            d=d[(d[0]==chrom)&(d[1]<end)&(d[2]>start)]; ax.broken_barh([(max(start,int(r[1])),min(end,int(r[2]))-max(start,int(r[1]))) for _,r in d.iterrows()],(j-.35,.7),label=Path(b).stem)
        if a.bed: ax.legend(fontsize=8,loc="upper right")
        save(fig,f"genome_browser_{name}")
if __name__=="__main__": main()
