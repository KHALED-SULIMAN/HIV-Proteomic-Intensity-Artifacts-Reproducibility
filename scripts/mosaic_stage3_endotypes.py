#!/usr/bin/env python3
# =============================================================================
# mosaic_stage3_endotypes.py
# Definitive integration for the biology/drug story: compare PCA vs MOFA-like,
# pick endotype number rigorously, test endotypes vs held-out clinical outcomes,
# and save the winning endotype assignment for downstream analysis.
# =============================================================================
import os, numpy as np, pandas as pd
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from scipy import stats

np.random.seed(0)
PROJ=os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
INT=os.path.join(PROJ,"data/processed/integrated"); OUT=os.path.join(PROJ,"results/mosaic")
os.makedirs(OUT,exist_ok=True)

def load(n): return pd.read_csv(os.path.join(INT,n),index_col=0)
rna,prot,meta=load("RNA_top_variable.csv"),load("proteomics_olink.csv"),load("metabolomics.csv")
outc=pd.read_csv(os.path.join(INT,"outcomes.csv")); outc.index=outc["COCOMO_ID"].astype(str)
ids=[str(i) for i in rna.index]; outc=outc.loc[ids]
assert [str(i) for i in prot.index]==ids==[str(i) for i in meta.index], "alignment!"
print(f"Aligned: {len(ids)} patients")

def prep(df,cap=1000):
    X=df.values.astype(np.float32); mu,sd=np.nanmean(X,0),np.nanstd(X,0)+1e-8
    X=np.nan_to_num((X-mu)/sd,nan=0.0)
    if cap and X.shape[1]>cap: X=X[:,np.argsort(X.var(0))[::-1][:cap]]
    return X
Xn={"rna":prep(rna),"prot":prep(prot),"meta":prep(meta)}
Xcat=np.concatenate([Xn[k] for k in Xn],1)

ym=outc["METS"].values.astype(float)
ytgl=(outc["Tgl"].values.astype(float)>=np.nanquantile(outc["Tgl"].values.astype(float),0.5)).astype(float)

def pca_emb(k=10): return PCA(k,random_state=0).fit_transform(Xcat)
def mofa_emb(k=4): return np.concatenate([TruncatedSVD(k,random_state=0).fit_transform(Xn[m]) for m in Xn],1)

def pick_k(emb):
    best=(2,-1)
    for k in range(2,7):
        lab=KMeans(k,n_init=10,random_state=0).fit_predict(emb)
        s=silhouette_score(emb,lab)
        print(f"    k={k}: silhouette {s:.3f}")
        if s>best[1]: best=(k,s)
    return best

def auroc_dist(emb,y,seeds=range(50)):
    m=~np.isnan(y);e=emb[m];yy=y[m].astype(int)
    if len(np.unique(yy))<2: return np.array([np.nan])
    out=[]
    for s in seeds:
        pr=np.zeros(len(yy))
        for tr,te in StratifiedKFold(5,shuffle=True,random_state=s).split(e,yy):
            pr[te]=LogisticRegression(max_iter=2000).fit(e[tr],yy[tr]).predict_proba(e[te])[:,1]
        out.append(roc_auc_score(yy,pr))
    return np.array(out)

print("\n=== PCA integration: choosing endotype number ===")
pe=pca_emb(); pk,ps=pick_k(pe)
print(f"  -> PCA best k = {pk} (silhouette {ps:.3f})")
print("\n=== MOFA-like integration: choosing endotype number ===")
me=mofa_emb(); mk,ms=pick_k(me)
print(f"  -> MOFA best k = {mk} (silhouette {ms:.3f})")

print("\n=== Held-out clinical stratification (AUROC, 50 splits) ===")
for name,emb in [("PCA",pe),("MOFA-like",me)]:
    dm=auroc_dist(emb,ym); dt=auroc_dist(emb,ytgl)
    print(f"  {name:10s}  METS {dm.mean():.3f} [{np.percentile(dm,2.5):.3f},{np.percentile(dm,97.5):.3f}]"
          f"   highTgl {dt.mean():.3f} [{np.percentile(dt,2.5):.3f},{np.percentile(dt,97.5):.3f}]")
dp,dmo=auroc_dist(pe,ym),auroc_dist(me,ym); diff=dp-dmo
lo,hi=np.percentile(diff,2.5),np.percentile(diff,97.5)
winner = "PCA" if (lo>0) else ("MOFA" if hi<0 else "TIE (pick PCA=simpler)")
print(f"\n  PCA - MOFA on METS: {diff.mean():+.3f} CI[{lo:+.3f},{hi:+.3f}] -> {winner}")

# ---- lock in the winning representation & endotypes ------------------------
use_emb, use_k, use_name = (pe,pk,"PCA") if "PCA" in winner or "TIE" in winner else (me,mk,"MOFA")
endotype = KMeans(use_k,n_init=20,random_state=0).fit_predict(use_emb)
res=outc.copy(); res["endotype"]=endotype
res.to_csv(os.path.join(OUT,"definitive_endotypes.csv"))
pd.DataFrame(use_emb,index=ids,
             columns=[f"f{i+1}" for i in range(use_emb.shape[1])]
            ).to_csv(os.path.join(OUT,"definitive_embedding.csv"))
print(f"\nLocked in: {use_name} integration, {use_k} endotypes.")
print("Endotype sizes:", dict(zip(*np.unique(endotype,return_counts=True))))

# quick clinical characterization of endotypes
print("\n=== Endotype clinical profile (mean) ===")
for c in ["METS","BMI","Tgl","Hdl","Ldl","CD4","AGE"]:
    if c in res: 
        gm=res.groupby("endotype")[c].mean()
        print(f"  {c:5s}", " ".join(f"E{int(g)}={v:.2f}" for g,v in gm.items()))
print("\nStage 3 complete. Saved definitive_endotypes.csv + definitive_embedding.csv")
