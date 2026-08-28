#!/usr/bin/env python3
# =============================================================================
# mosaic_stage6_pathways.py
# Differential proteins on BATCH-CORRECTED data + pathway enrichment.
# Output: corrected differential table, up/down protein lists, enriched pathways.
# Requires: pip install gseapy  (queries Enrichr; needs internet)
# =============================================================================
import os, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from scipy import stats

PROJ=os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
INT=os.path.join(PROJ,"data/processed/integrated"); OUT=os.path.join(PROJ,"results/mosaic")

prot=pd.read_csv(os.path.join(INT,"proteomics_olink.csv"),index_col=0)
meta=pd.read_csv(os.path.join(INT,"metabolomics.csv"),index_col=0)
outc=pd.read_csv(os.path.join(INT,"outcomes.csv")); outc.index=outc["COCOMO_ID"].astype(str)
endo=pd.read_csv(os.path.join(OUT,"endotypes_batch_corrected.csv"))
endo.index=endo["COCOMO_ID"].astype(str)
ids=[str(i) for i in prot.index]; outc=outc.loc[ids]; endo=endo.loc[ids]
g=endo["endotype_corrected"].values

# batch-corrected proteomics (same median-centering)
P=prot.values.astype(float); P=P-np.nanmedian(P,1)[:,None]

# define high-risk group by METS rate
mets=outc["METS"].values.astype(float)
rate={k: np.nanmean(mets[g==k]) for k in np.unique(g)}
high=max(rate,key=rate.get); low=min(rate,key=rate.get)
gi=(g==high).astype(int)
print(f"High-risk = E{high} (METS {rate[high]:.2f}, n={(g==high).sum()}) | "
      f"low-risk = E{low} (METS {rate[low]:.2f}, n={(g==low).sum()})")

def bh(pv):
    pv=np.asarray(pv,float);n=len(pv);o=np.argsort(pv);r=pv[o]
    q=r*n/(np.arange(n)+1);q=np.minimum.accumulate(q[::-1])[::-1]
    out=np.empty(n);out[o]=np.clip(q,0,1);return out

def diff(X, names, layer):
    a=X[gi==1]; b=X[gi==0]
    p=np.array([stats.mannwhitneyu(a[:,j],b[:,j],alternative="two-sided").pvalue
                for j in range(X.shape[1])])
    d=np.nanmean(a,0)-np.nanmean(b,0); q=bh(p)
    r=pd.DataFrame({"feature":names,"layer":layer,"diff_high_minus_low":d,"p":p,"q":q})
    print(f"  {layer:14s}: {len(names):5d} tested | {(q<0.05).sum():4d} sig (FDR<0.05)")
    return r.sort_values("q")

print("\n=== Differential features (batch-corrected) ===")
res=pd.concat([diff(P,list(prot.columns),"proteomics"),
               diff(meta.values.astype(float),list(meta.columns),"metabolomics")],
              ignore_index=True)
res.to_csv(os.path.join(OUT,"corrected_differential_features.csv"),index=False)

sig=res[(res.q<0.05)&(res.layer=="proteomics")]
up=sig[sig.diff_high_minus_low>0]["feature"].tolist()
dn=sig[sig.diff_high_minus_low<0]["feature"].tolist()
print(f"\nProteins UP in high-risk: {len(up)} | DOWN: {len(dn)}")
pd.Series(up).to_csv(os.path.join(OUT,"proteins_up_highrisk.txt"),index=False,header=False)
pd.Series(dn).to_csv(os.path.join(OUT,"proteins_down_highrisk.txt"),index=False,header=False)
print("Top 10 UP:  ", ", ".join(up[:10]))
print("Top 10 DOWN:", ", ".join(dn[:10]))

# ---- pathway enrichment -----------------------------------------------------
print("\n=== Pathway enrichment (Enrichr) ===")
try:
    import gseapy as gp
    libs=["KEGG_2021_Human","GO_Biological_Process_2023","Reactome_2022"]
    for label,genes in [("UP_in_high_risk",up),("DOWN_in_high_risk",dn)]:
        if len(genes)<5:
            print(f"  {label}: too few genes ({len(genes)}), skipped"); continue
        try:
            e=gp.enrichr(gene_list=genes, gene_sets=libs, organism="Human",
                         outdir=None, no_plot=True)
            r=e.results
            r=r[r["Adjusted P-value"]<0.05].sort_values("Adjusted P-value")
            r.to_csv(os.path.join(OUT,f"pathways_{label}.csv"),index=False)
            print(f"\n  {label}: {len(r)} enriched terms (FDR<0.05). Top 8:")
            for _,x in r.head(8).iterrows():
                print(f"    [{x['Gene_set'][:18]:18s}] {x['Term'][:55]:55s} q={x['Adjusted P-value']:.2e}")
        except Exception as ex:
            print(f"  {label}: enrichr failed ({ex})")
except ImportError:
    print("  gseapy not installed. Run:  pip install gseapy")
    print("  (protein lists saved; you can paste them into https://maayanlab.cloud/Enrichr/)")

print("\nStage 6 complete. Saved corrected_differential_features.csv + protein lists.")
