#!/usr/bin/env python3
"""Build merged consensus regions supported by reciprocal peak overlap."""

from __future__ import annotations
import argparse
from pathlib import Path


def read_bed(path: Path):
    rows=[]
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"): continue
        f=line.split("\t"); rows.append((f[0],int(f[1]),int(f[2])))
    return sorted(rows)


def supported(a,b,fraction):
    by={}
    for row in b: by.setdefault(row[0],[]).append(row)
    out=[]
    left_by_chrom={}
    for chrom,start,end in a:
        candidates=by.get(chrom,[])
        left=left_by_chrom.get(chrom,0)
        while left<len(candidates) and candidates[left][2]<=start:
            left+=1
        left_by_chrom[chrom]=left
        j=left
        while j<len(candidates):
            _,bs,be=candidates[j]
            if bs>=end: break
            ov=min(end,be)-max(start,bs)
            if ov>0 and ov/(end-start)>=fraction and ov/(be-bs)>=fraction:
                out.append((chrom,start,end)); break
            j+=1
    return out


def merge(rows):
    out=[]
    for chrom,start,end in sorted(set(rows)):
        if out and out[-1][0]==chrom and start<=out[-1][2]:
            out[-1]=(chrom,out[-1][1],max(end,out[-1][2]))
        else: out.append((chrom,start,end))
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--rep1",type=Path,required=True)
    p.add_argument("--rep2",type=Path,required=True); p.add_argument("--output",type=Path,required=True)
    p.add_argument("--reciprocal-fraction",type=float,default=0.5); a=p.parse_args()
    one,two=read_bed(a.rep1),read_bed(a.rep2)
    consensus=merge(supported(one,two,a.reciprocal_fraction)+supported(two,one,a.reciprocal_fraction))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text("".join(f"{c}\t{s}\t{e}\n" for c,s,e in consensus))
    print(f"rep1={len(one)} rep2={len(two)} reciprocal_merged_consensus={len(consensus)} output={a.output}")


if __name__=="__main__": main()
