#!/usr/bin/env python3
# =============================================================================
# mosaic_stage4b_artifact_check.py
# Is the proteomic endotype signature real biology or a technical artifact?
# Runs four honest diagnostics:
#   (1) Global-shift test: do "significant" proteins reflect a uniform per-sample
#       offset (classic batch effect) rather than specific biology?
#   (2) Survival test: does the endotype persist WITHOUT proteomics
#       (RNA + metabolomics only)? If not, proteomics drove everything.
#   (3) Per-layer endotype agreement (which layer really defines the split).
#   (4) Direction balance: batch shifts push most features one way; biology mixes.
# =============================================================================
import os, numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from scipy import stats

np.random.seed(0)
PROJ=os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
INT=os.path.join(PROJ,"data/processed/integrated"); OUT=os.path.join(PROJ,"results/mosaic")

def load(n): return pd.read_csv(os.path.join(INT,n),index_col=0)
rna,prot,meta=load("RNA_top_variable.csv"),load("proteomics_olink.csv"),load("metabolomics.csv")
endo=pd.read_csv(os.path.join(OUT,"definitive_endotypes.csv")); endo.index=endo["COCOMO_ID"].astype(str)
ids=[str(i) for i in rna.index]; endo=endo.loc[ids]
g=endo["endotype"].values

def prep(df,cap=1000):
    X=df.values.astype(float); mu,sd=np.nanmean(X,0),np.nanstd(X,0)+1e-8
    X=np.nan_to_num((X-mu)/sd,nan=0.0)
    if cap and X.shape[1]>cap: X=X[:,np.argsort(X.var(0))[::-1][:cap]]
    return X
Pn,Mn,Rn=prep(prot),prep(meta),prep(rna)
prot_raw=prot.values.astype(float)  # for median test use raw (pre-zscore) values

print("="*64)
print("ARTIFACT DIAGNOSTICS for the proteomic endotype signature")
print("="*64)

# (1) global-shift / batch signature: per-sample median protein level
med=np.nanmedian(prot_raw,1)
p_med=stats.mannwhitneyu(med[g==0],med[g==1],alternative="two-sided").pvalue
print(f"\n(1) Per-sample MEDIAN protein level by endotype:")
print(f"    E0 median-of-medians {np.median(med[g==0]):.3f} | E1 {np.median(med[g==1]):.3f} | p={p_med:.2e}")
print("    -> very small p suggests a GLOBAL SHIFT (batch/normalization), not specific biology")

# (4) direction balance among significant proteins
def bh(pv):
    pv=np.asarray(pv,float);n=len(pv);o=np.argsort(pv);r=pv[o]
    q=r*n/(np.arange(n)+1);q=np.minimum.accumulate(q[::-1])[::-1];out=np.empty(n);out[o]=np.clip(q,0,1);return out
X=prot.values.astype(float)
p=np.array([stats.mannwhitneyu(X[g==0,j],X[g==1,j],alternative="two-sided").pvalue for j in range(X.shape[1])])
lfc=np.nanmean(X[g==0],0)-np.nanmean(X[g==1],0)
q=bh(p); sigmask=q<0.05
up=(lfc[sigmask]>0).sum(); dn=(lfc[sigmask]<0).sum()
print(f"\n(4) Direction of significant proteins: {up} UP vs {dn} DOWN in high-risk")
print("    -> heavily one-sided (e.g. >90% one direction) points to a global shift/batch;")
print("       a real biological contrast is usually more balanced.")

# (2) survival test: endotype WITHOUT proteomics
def cluster(Xcat,k=2):
    return KMeans(k,n_init=20,random_state=0).fit_predict(PCA(min(10,Xcat.shape[1]),random_state=0).fit_transform(Xcat))
endo_full   = cluster(np.concatenate([Pn,Mn,Rn],1))
endo_no_prot= cluster(np.concatenate([Mn,Rn],1))
endo_prot   = cluster(Pn)
ari_np  = adjusted_rand_score(g, endo_no_prot)
ari_p   = adjusted_rand_score(g, endo_prot)
print(f"\n(2/3) Which layer defines the endotype? (agreement with final endotype, ARI)")
print(f"    proteomics-only ARI:        {ari_p:.3f}   (near 1 => proteomics defines it)")
print(f"    RNA+metabolomics-only ARI:  {ari_np:.3f}   (near 0 => split vanishes without proteomics)")

# verdict
print("\n"+"="*64)
print("HONEST VERDICT")
print("="*64)
flags=0
if p_med<1e-3: flags+=1
if max(up,dn)/(up+dn+1e-9)>0.9: flags+=1
if ari_np<0.2: flags+=1
if ari_p>0.7: flags+=1
print(f"warning flags raised: {flags}/4")
if flags>=3:
    print("STRONG artifact suspicion: the proteomic split looks largely technical.")
    print("Recommended: batch-correct proteomics (or drop it) and re-derive endotypes")
    print("from RNA+metabolomics, OR find and regress out the batch variable.")
elif flags>=1:
    print("PARTIAL concern: some technical contribution likely; proceed with caution")
    print("and report the check. Consider batch correction as a sensitivity analysis.")
else:
    print("No strong artifact signature detected: the proteomic differences appear")
    print("biologically driven. Safe to proceed to pathways / drug repurposing.")
print("\nStage 4b complete.")
