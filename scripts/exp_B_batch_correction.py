#!/usr/bin/env python3
# =============================================================================
# exp_B_batch_correction.py   -> FIGURE 7
# Which batch-correction strategy best removes the Olink global-shift artifact
# while preserving clinical signal?
#   none | median-centering | quantile normalisation | rank-inverse-normal
# Evaluated on: artifact flags (0-4), endotype stability, METS AUROC,
#               up/down direction balance.
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

def load(n): return pd.read_csv(os.path.join(INT,n),index_col=0)
rna,prot,meta=load("RNA_top_variable.csv"),load("proteomics_olink.csv"),load("metabolomics.csv")
outc=pd.read_csv(os.path.join(INT,"outcomes.csv")); outc.index=outc["COCOMO_ID"].astype(str)
ids=[str(i) for i in rna.index]; outc=outc.loc[ids]
y=outc["METS"].values.astype(float)
P0=prot.values.astype(float)

def none_(X): return X.copy()
def median_(X): return X-np.nanmedian(X,1)[:,None]
def quantile_(X):
    R=np.argsort(np.argsort(X,axis=1),axis=1)
    ref=np.sort(X,axis=1).mean(0)
    return ref[R]
def rankinv_(X):
    out=np.empty_like(X)
    for i in range(X.shape[0]):
        r=stats.rankdata(X[i])/(X.shape[1]+1); out[i]=stats.norm.ppf(r)
    return out
METHODS={"none":none_,"median_center":median_,"quantile":quantile_,"rank_inverse_normal":rankinv_}

def z(X,cap=1000):
    mu,sd=np.nanmean(X,0),np.nanstd(X,0)+1e-8
    X=np.nan_to_num((X-mu)/sd,nan=0.0)
    if X.shape[1]>cap: X=X[:,np.argsort(X.var(0))[::-1][:cap]]
    return X
Mz,Rz=z(meta.values.astype(float)),z(rna.values.astype(float))

def bh(pv):
    pv=np.asarray(pv,float);n=len(pv);o=np.argsort(pv);r=pv[o]
    q=r*n/(np.arange(n)+1);q=np.minimum.accumulate(q[::-1])[::-1]
    out=np.empty(n);out[o]=np.clip(q,0,1);return out

def auroc(E,yy,seeds=range(30)):
    m=~np.isnan(yy);E=E[m];v=yy[m].astype(int)
    if len(np.unique(v))<2: return np.nan
    o=[]
    for s in seeds:
        pr=np.zeros(len(v))
        for tr,te in StratifiedKFold(5,shuffle=True,random_state=s).split(E,v):
            pr[te]=LogisticRegression(max_iter=2000).fit(E[tr],v[tr]).predict_proba(E[te])[:,1]
        o.append(roc_auc_score(v,pr))
    return float(np.mean(o))

rows=[]
print(f"{'method':22s} {'flags':>5} {'up/down':>12} {'silh':>6} {'AUROC':>7} {'stability':>9}")
for name,fn in METHODS.items():
    Pc=fn(P0); Pz=z(Pc)
    emb=PCA(10,random_state=0).fit_transform(np.concatenate([Pz,Mz,Rz],1))
    g=KMeans(2,n_init=20,random_state=0).fit_predict(emb)
    sil=silhouette_score(emb,g)
    # artifact diagnostics
    med=np.nanmedian(Pc,1); p_med=stats.mannwhitneyu(med[g==0],med[g==1]).pvalue
    p=np.array([stats.mannwhitneyu(Pc[g==0,j],Pc[g==1,j]).pvalue for j in range(Pc.shape[1])])
    lfc=np.nanmean(Pc[g==0],0)-np.nanmean(Pc[g==1],0); q=bh(p); s=q<0.05
    up=int((lfc[s]>0).sum()); dn=int((lfc[s]<0).sum())
    ari_np=adjusted_rand_score(g,KMeans(2,n_init=20,random_state=0).fit_predict(
            PCA(10,random_state=0).fit_transform(np.concatenate([Mz,Rz],1))))
    ari_p=adjusted_rand_score(g,KMeans(2,n_init=20,random_state=0).fit_predict(
            PCA(10,random_state=0).fit_transform(Pz)))
    flags=int(p_med<1e-3)+int((max(up,dn)/(up+dn+1e-9))>0.9 if up+dn>0 else 0)+ \
          int(ari_np<0.2)+int(ari_p>0.7)
    au=auroc(emb,y)
    # stability: ARI across bootstrap resamples
    aris=[]
    for b in range(10):
        idx=np.random.RandomState(b).choice(len(ids),len(ids),replace=True)
        gb=KMeans(2,n_init=10,random_state=0).fit_predict(emb[idx])
        aris.append(adjusted_rand_score(g[idx],gb))
    stab=float(np.mean(aris))
    rows.append({"method":name,"artifact_flags":flags,"n_up":up,"n_down":dn,
                 "silhouette":sil,"METS_AUROC":au,"stability_ARI":stab,
                 "p_median_shift":p_med})
    print(f"{name:22s} {flags:5d} {f'{up}/{dn}':>12} {sil:6.3f} {au:7.3f} {stab:9.3f}")

df=pd.DataFrame(rows); df.to_csv(os.path.join(OUT,"expB_batch_correction.csv"),index=False)

try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,3,figsize=(14,4))
    ax[0].bar(df["method"],df["artifact_flags"],color="crimson"); ax[0].set_ylabel("artifact flags (0-4)")
    ax[0].set_title("Artifact removal"); ax[0].tick_params(axis='x',rotation=30)
    ax[1].bar(df["method"],df["METS_AUROC"],color="steelblue"); ax[1].set_ylabel("METS AUROC")
    ax[1].set_ylim(0.5,1); ax[1].set_title("Clinical signal retained"); ax[1].tick_params(axis='x',rotation=30)
    ax[2].bar(df["method"],df["stability_ARI"],color="seagreen"); ax[2].set_ylabel("bootstrap ARI")
    ax[2].set_title("Endotype stability"); ax[2].tick_params(axis='x',rotation=30)
    plt.tight_layout(); plt.savefig(os.path.join(OUT,"FIG7_batch_correction.png"),dpi=150)
    print("\nSaved FIG7_batch_correction.png")
except Exception as e: print("plot skipped:",e)

best=df[(df.artifact_flags<=1)].sort_values("METS_AUROC",ascending=False)
print("\n=== RECOMMENDATION ===")
if len(best):
    b=best.iloc[0]
    print(f"Best: {b['method']} -> flags {int(b['artifact_flags'])}/4, "
          f"AUROC {b['METS_AUROC']:.3f}, stability {b['stability_ARI']:.3f}")
else:
    print("No method reduced artifact flags to <=1. Consider dropping proteomics.")
