#!/usr/bin/env python3
# =============================================================================
# mosaic_stage5_fix_and_verify.py   -- ALL-IN-ONE
# 1) inspect raw Olink file for batch/plate info
# 2) remove the global per-sample shift (median-centering, a standard Olink fix)
# 3) re-derive endotypes on CORRECTED data
# 4) re-run the full artifact check
# 5) re-test clinical stratification (METS) vs baselines
# Verdict at the end tells you if real biology survives.
# =============================================================================
import os, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score, roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from scipy import stats

np.random.seed(0)
PROJ=os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
INT=os.path.join(PROJ,"data/processed/integrated"); OUT=os.path.join(PROJ,"results/mosaic")
REPO=os.path.join(PROJ,"external/Immunometabolism_HIV")

print("="*66); print("STEP 1: look for batch/plate info in raw Olink files"); print("="*66)
for f in ["data/Olink/results_master_file_future_olink.xlsx",
          "data/Olink/Preliminary Olink Explore 3072 Assay List_2021-06-16.xlsx"]:
    p=os.path.join(REPO,f)
    if os.path.exists(p):
        try:
            xl=pd.ExcelFile(p); print(f"\n{os.path.basename(f)} sheets: {xl.sheet_names[:5]}")
            d=pd.read_excel(p, sheet_name=xl.sheet_names[0], nrows=3)
            cols=[c for c in d.columns]
            hits=[c for c in cols if any(k in str(c).lower() for k in
                  ["plate","batch","run","panel","block","chip","date","qc"])]
            print(f"  columns: {len(cols)} | batch-like columns found: {hits if hits else 'NONE'}")
        except Exception as e: print(f"  could not read: {e}")
    else: print(f"  missing: {f}")

print("\n"+"="*66); print("STEP 2: remove global per-sample shift (median-centering)"); print("="*66)
def load(n): return pd.read_csv(os.path.join(INT,n),index_col=0)
rna,prot,meta=load("RNA_top_variable.csv"),load("proteomics_olink.csv"),load("metabolomics.csv")
outc=pd.read_csv(os.path.join(INT,"outcomes.csv")); outc.index=outc["COCOMO_ID"].astype(str)
ids=[str(i) for i in rna.index]; outc=outc.loc[ids]

P=prot.values.astype(float)
med_before=np.nanmedian(P,1)
Pc = P - med_before[:,None]          # remove per-sample global offset
print(f"  per-sample median spread BEFORE: {med_before.std():.4f}")
print(f"  per-sample median spread AFTER : {np.nanmedian(Pc,1).std():.4f}")

def z(X,cap=1000):
    mu,sd=np.nanmean(X,0),np.nanstd(X,0)+1e-8
    X=np.nan_to_num((X-mu)/sd,nan=0.0)
    if cap and X.shape[1]>cap: X=X[:,np.argsort(X.var(0))[::-1][:cap]]
    return X
Pz,Mz,Rz = z(Pc), z(meta.values.astype(float)), z(rna.values.astype(float))

print("\n"+"="*66); print("STEP 3: re-derive endotypes on CORRECTED data"); print("="*66)
Xcat=np.concatenate([Pz,Mz,Rz],1)
emb=PCA(10,random_state=0).fit_transform(Xcat)
best=(2,-1)
for k in range(2,6):
    lab=KMeans(k,n_init=20,random_state=0).fit_predict(emb); s=silhouette_score(emb,lab)
    print(f"  k={k}: silhouette {s:.3f}")
    if s>best[1]: best=(k,s)
K=best[0]; endo=KMeans(K,n_init=20,random_state=0).fit_predict(emb)
print(f"  -> k={K} (silhouette {best[1]:.3f}), sizes {dict(zip(*np.unique(endo,return_counts=True)))}")

print("\n"+"="*66); print("STEP 4: re-run ARTIFACT CHECK on corrected endotypes"); print("="*66)
def bh(pv):
    pv=np.asarray(pv,float);n=len(pv);o=np.argsort(pv);r=pv[o]
    q=r*n/(np.arange(n)+1);q=np.minimum.accumulate(q[::-1])[::-1];out=np.empty(n);out[o]=np.clip(q,0,1);return out
g=endo
p=np.array([stats.mannwhitneyu(Pc[g==0,j],Pc[g==1,j],alternative="two-sided").pvalue
            for j in range(Pc.shape[1])])
lfc=np.nanmean(Pc[g==0],0)-np.nanmean(Pc[g==1],0); q=bh(p); sig=q<0.05
up=(lfc[sig]>0).sum(); dn=(lfc[sig]<0).sum()
medc=np.nanmedian(Pc,1); p_med=stats.mannwhitneyu(medc[g==0],medc[g==1]).pvalue
ari_np=adjusted_rand_score(g, KMeans(K,n_init=20,random_state=0).fit_predict(
        PCA(10,random_state=0).fit_transform(np.concatenate([Mz,Rz],1))))
ari_p =adjusted_rand_score(g, KMeans(K,n_init=20,random_state=0).fit_predict(
        PCA(10,random_state=0).fit_transform(Pz)))
print(f"  significant proteins: {sig.sum()} ({up} up / {dn} down)")
print(f"  per-sample median by endotype p={p_med:.2e}")
print(f"  ARI proteomics-only {ari_p:.3f} | RNA+metab-only {ari_np:.3f}")
flags=sum([p_med<1e-3, (max(up,dn)/(up+dn+1e-9))>0.9 if (up+dn)>0 else False,
           ari_np<0.2, ari_p>0.7])
print(f"  warning flags now: {flags}/4  (was 4/4 before correction)")

print("\n"+"="*66); print("STEP 5: clinical stratification after correction"); print("="*66)
ym=outc["METS"].values.astype(float)
def auroc(E,y,seeds=range(30)):
    m=~np.isnan(y);e=E[m];yy=y[m].astype(int)
    if len(np.unique(yy))<2: return np.array([np.nan])
    o=[]
    for s in seeds:
        pr=np.zeros(len(yy))
        for tr,te in StratifiedKFold(5,shuffle=True,random_state=s).split(e,yy):
            pr[te]=LogisticRegression(max_iter=2000).fit(e[tr],yy[tr]).predict_proba(e[te])[:,1]
        o.append(roc_auc_score(yy,pr))
    return np.array(o)
d_corr=auroc(emb,ym)
emb_nop=PCA(10,random_state=0).fit_transform(np.concatenate([Mz,Rz],1))
d_nop=auroc(emb_nop,ym)
print(f"  corrected all-layers  METS AUROC {d_corr.mean():.3f} [{np.percentile(d_corr,2.5):.3f},{np.percentile(d_corr,97.5):.3f}]")
print(f"  RNA+metab only        METS AUROC {d_nop.mean():.3f} [{np.percentile(d_nop,2.5):.3f},{np.percentile(d_nop,97.5):.3f}]")

pd.DataFrame({"COCOMO_ID":ids,"endotype_corrected":endo}).to_csv(
    os.path.join(OUT,"endotypes_batch_corrected.csv"),index=False)

print("\n"+"="*66); print("VERDICT"); print("="*66)
if flags<=1 and d_corr.mean()>0.65:
    print("GOOD: artifact largely removed AND endotypes still stratify METS.")
    print("-> Proceed: this is a defensible, publishable result.")
elif flags<=1:
    print("MIXED: artifact removed, but clinical stratification is weak.")
    print("-> Honest framing: no strong multi-omics endotype in this cohort.")
else:
    print("STILL ARTIFACTUAL: correction insufficient.")
    print("-> Recommend dropping proteomics; use RNA+metabolomics only, or report")
    print("   the artifact itself as the paper's methodological contribution.")
print("\nSaved endotypes_batch_corrected.csv")
