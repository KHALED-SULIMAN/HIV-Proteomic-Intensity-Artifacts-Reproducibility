#!/usr/bin/env python3
# =============================================================================
# exp_A_missing_modality.py   -> FIGURE 6
# Does MOSAIC's native missing-block handling beat imputation when whole
# omics layers are missing?  (PCA won on complete data - but PCA cannot
# handle missing blocks without imputation.)
#
# Compares at 0/10/20/30/40/50% of patients missing one whole random layer:
#   MOSAIC (PoE masking, no imputation) | PCA+mean | PCA+KNN | complete-case
# Metric: cross-validated AUROC for held-out metabolic syndrome.
# =============================================================================
import os, numpy as np, pandas as pd, torch, torch.nn as nn, warnings
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer, SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

np.random.seed(0); torch.manual_seed(0)
PROJ=os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
INT=os.path.join(PROJ,"data/processed/integrated"); OUT=os.path.join(PROJ,"results/mosaic")
LAT=8; HID=64; EPOCHS=250; FRACS=[0.0,0.1,0.2,0.3,0.4,0.5]; REPS=5

def load(n): return pd.read_csv(os.path.join(INT,n),index_col=0)
rna,prot,meta=load("RNA_top_variable.csv"),load("proteomics_olink.csv"),load("metabolomics.csv")
outc=pd.read_csv(os.path.join(INT,"outcomes.csv")); outc.index=outc["COCOMO_ID"].astype(str)
ids=[str(i) for i in rna.index]; outc=outc.loc[ids]
P=prot.values.astype(float); P=P-np.nanmedian(P,1)[:,None]   # batch-corrected

def z(X,cap=800):
    mu,sd=np.nanmean(X,0),np.nanstd(X,0)+1e-8
    X=np.nan_to_num((X-mu)/sd,nan=0.0)
    if X.shape[1]>cap: X=X[:,np.argsort(X.var(0))[::-1][:cap]]
    return X.astype(np.float32)
Xn={"rna":z(rna.values.astype(float)),"prot":z(P),"meta":z(meta.values.astype(float))}
keys=list(Xn); N=len(ids); L=len(keys)
y=outc["METS"].values.astype(float); ok=~np.isnan(y); y=y[ok].astype(int)
print(f"n={N} patients, layers={ {k:Xn[k].shape[1] for k in Xn} }, METS available n={ok.sum()}")

def make_mask(frac,seed):
    rng=np.random.RandomState(seed); M=np.ones((N,L),bool)
    for p in rng.choice(N,int(frac*N),replace=False): M[p,rng.randint(0,L)]=False
    return M

class Enc(nn.Module):
    def __init__(s,di):
        super().__init__(); s.n=nn.Sequential(nn.Linear(di,HID),nn.ReLU(),nn.Dropout(0.2))
        s.mu=nn.Linear(HID,LAT); s.lv=nn.Linear(HID,LAT)
    def forward(s,x): h=s.n(x); return s.mu(h),s.lv(h)

def mosaic(M):
    T_=[torch.tensor(Xn[k]) for k in keys]
    enc=nn.ModuleList([Enc(Xn[k].shape[1]) for k in keys])
    dec=nn.ModuleList([nn.Sequential(nn.Linear(LAT,HID),nn.ReLU(),
                        nn.Linear(HID,Xn[k].shape[1])) for k in keys])
    opt=torch.optim.Adam(list(enc.parameters())+list(dec.parameters()),
                         lr=1e-3,weight_decay=1e-4)
    Mt=torch.tensor(M.astype(np.float32))
    for e in range(EPOCHS):
        opt.zero_grad(); mus=[];lvs=[]
        for i in range(L):
            m,l=enc[i](T_[i]); mus.append(m); lvs.append(l)
        mus=torch.stack(mus,0); lvs=torch.stack(lvs,0)
        prec=torch.exp(-lvs)*Mt.T.unsqueeze(-1)
        Ts=1.0+prec.sum(0); mu=(mus*prec).sum(0)/Ts; lv=torch.log(1.0/Ts)
        zz=mu+torch.exp(0.5*lv)*torch.randn_like(lv)
        rec=0
        for i in range(L):
            w=Mt[:,i:i+1]
            rec=rec+(((dec[i](zz)-T_[i])**2)*w).sum()/(w.sum()*T_[i].shape[1]+1e-8)
        kl=torch.clamp((-0.5*(1+lv-mu.pow(2)-lv.exp())).mean(0),min=0.5).sum()
        (rec+min(1,e/60)*kl).backward(); opt.step()
    with torch.no_grad():
        mus=[];lvs=[]
        for i in range(L):
            m,l=enc[i](T_[i]); mus.append(m); lvs.append(l)
        mus=torch.stack(mus,0); lvs=torch.stack(lvs,0)
        prec=torch.exp(-lvs)*Mt.T.unsqueeze(-1)
        Ts=1.0+prec.sum(0); mu=(mus*prec).sum(0)/Ts
    return mu.numpy()

def cat_with_missing(M,impute):
    parts=[]
    for i,k in enumerate(keys):
        A=Xn[k].astype(float).copy(); A[~M[:,i]]=np.nan; parts.append(A)
    C=np.concatenate(parts,1)
    if impute=="mean": C=SimpleImputer().fit_transform(C)
    elif impute=="knn": C=KNNImputer(n_neighbors=5).fit_transform(C)
    return C

def auroc(E,yy,seeds=range(20)):
    if len(np.unique(yy))<2 or len(yy)<25: return np.nan
    o=[]
    for s in seeds:
        pr=np.zeros(len(yy))
        for tr,te in StratifiedKFold(5,shuffle=True,random_state=s).split(E,yy):
            pr[te]=LogisticRegression(max_iter=2000).fit(E[tr],yy[tr]).predict_proba(E[te])[:,1]
        o.append(roc_auc_score(yy,pr))
    return float(np.mean(o))

rows=[]
print("\nfrac | MOSAIC | PCA+mean | PCA+KNN | complete-case (n)")
for frac in FRACS:
    accs={"MOSAIC":[],"PCA_mean":[],"PCA_knn":[],"complete_case":[]}
    ncc=[]
    for rep in range(REPS):
        M=make_mask(frac,seed=100+rep)
        accs["MOSAIC"].append(auroc(mosaic(M)[ok],y))
        accs["PCA_mean"].append(auroc(PCA(LAT,random_state=0).fit_transform(
                                cat_with_missing(M,"mean"))[ok],y))
        accs["PCA_knn"].append(auroc(PCA(LAT,random_state=0).fit_transform(
                                cat_with_missing(M,"knn"))[ok],y))
        keep=M.all(1)&ok
        ncc.append(keep.sum())
        if keep.sum()>25:
            E=PCA(LAT,random_state=0).fit_transform(cat_with_missing(M,"mean"))[keep]
            accs["complete_case"].append(auroc(E,outc["METS"].values.astype(float)[keep].astype(int)))
        else: accs["complete_case"].append(np.nan)
    r={"missing_frac":frac,"n_complete":int(np.mean(ncc))}
    for k in accs: r[k]=float(np.nanmean(accs[k]))
    rows.append(r)
    print(f"{frac:4.1f} | {r['MOSAIC']:.3f}  | {r['PCA_mean']:.3f}    | "
          f"{r['PCA_knn']:.3f}   | {r['complete_case']:.3f} ({r['n_complete']})")

df=pd.DataFrame(rows); df.to_csv(os.path.join(OUT,"expA_missing_modality.csv"),index=False)

try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.figure(figsize=(8,5))
    for c,mk_,lb in [("MOSAIC","o-","MOSAIC (native masking)"),
                     ("PCA_mean","s--","PCA + mean imputation"),
                     ("PCA_knn","^--","PCA + KNN imputation"),
                     ("complete_case","d:","complete-case only")]:
        plt.plot(df["missing_frac"]*100, df[c], mk_, label=lb)
    plt.xlabel("% patients missing one whole omics layer")
    plt.ylabel("cross-validated AUROC (metabolic syndrome)")
    plt.title("Robustness to missing omics blocks")
    plt.legend(); plt.grid(alpha=.3); plt.tight_layout()
    plt.savefig(os.path.join(OUT,"FIG6_missing_modality.png"),dpi=150)
    print("\nSaved FIG6_missing_modality.png")
except Exception as e: print("plot skipped:",e)

d0,d5=df.iloc[0],df.iloc[-1]
print("\n=== VERDICT ===")
print(f"MOSAIC   {d0['MOSAIC']:.3f} -> {d5['MOSAIC']:.3f}  (drop {d0['MOSAIC']-d5['MOSAIC']:+.3f})")
print(f"PCA+mean {d0['PCA_mean']:.3f} -> {d5['PCA_mean']:.3f}  (drop {d0['PCA_mean']-d5['PCA_mean']:+.3f})")
if d5['MOSAIC']>d5['PCA_mean'] and d5['MOSAIC']>d5['PCA_knn']:
    print("-> MOSAIC degrades more gracefully: a REAL niche for the deep model.")
else:
    print("-> No advantage for MOSAIC even with missing blocks. Report honestly.")
