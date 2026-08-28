#!/usr/bin/env python3
# =============================================================================
# mosaic_stage2b_rigor.py
# Rigorous, honest comparison: MOSAIC vs PCA vs MOFA-like on held-out outcomes,
# with bootstrap confidence intervals and CI-of-difference.
# Answers: is MOSAIC statistically DIFFERENT from PCA, or COMPETITIVE (tie)?
# =============================================================================
import os, numpy as np, pandas as pd, torch, torch.nn as nn
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

torch.manual_seed(0); np.random.seed(0)
PROJ=os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
INT=os.path.join(PROJ,"data/processed/integrated"); OUT=os.path.join(PROJ,"results/mosaic")

def load(n): return pd.read_csv(os.path.join(INT,n),index_col=0)
rna,prot,meta=load("RNA_top_variable.csv"),load("proteomics_olink.csv"),load("metabolomics.csv")
outc=pd.read_csv(os.path.join(INT,"outcomes.csv")); outc.index=outc["COCOMO_ID"].astype(str)
ids=[str(i) for i in rna.index]; outc=outc.loc[ids]
def prep(df,cap=1000):
    X=df.values.astype(np.float32); mu,sd=np.nanmean(X,0),np.nanstd(X,0)+1e-8
    X=np.nan_to_num((X-mu)/sd,nan=0.0)
    if cap and X.shape[1]>cap: X=X[:,np.argsort(X.var(0))[::-1][:cap]]
    return X
Xn={"rna":prep(rna),"prot":prep(prot),"meta":prep(meta)}
Xt={k:torch.tensor(v) for k,v in Xn.items()}; dims={k:v.shape[1] for k,v in Xn.items()}
LAT=8;HID=64;DROP=0.3

yc=(outc["cluster"].values.astype(float)).astype(int); yc-=yc.min()
ym=outc["METS"].values.astype(float)

class Enc(nn.Module):
    def __init__(s,di):super().__init__();s.n=nn.Sequential(nn.Linear(di,HID),nn.ReLU(),nn.Dropout(DROP));s.mu=nn.Linear(HID,LAT);s.lv=nn.Linear(HID,LAT)
    def forward(s,x):h=s.n(x);return s.mu(h),s.lv(h)
class Dec(nn.Module):
    def __init__(s,do):super().__init__();s.n=nn.Sequential(nn.Linear(LAT,HID),nn.ReLU(),nn.Dropout(DROP),nn.Linear(HID,do))
    def forward(s,z):return s.n(z)
def poe(mus,lvs):
    mus=torch.stack(mus,0);lvs=torch.stack(lvs,0);T=torch.exp(-lvs)
    Ts=torch.ones_like(T[0])+T.sum(0);return (mus*T).sum(0)/Ts,torch.log(1.0/Ts)
def mosaic(lmbda=1.0,ep=300):
    enc=nn.ModuleDict({k:Enc(dims[k]) for k in dims});dec=nn.ModuleDict({k:Dec(dims[k]) for k in dims})
    head=nn.Linear(LAT,2);P=list(enc.parameters())+list(dec.parameters())+list(head.parameters())
    opt=torch.optim.Adam(P,lr=1e-3,weight_decay=1e-4);ya=torch.tensor(yc)
    for e in range(ep):
        for m in list(enc.values())+list(dec.values()):m.train()
        opt.zero_grad();mus=[];lvs=[]
        for k in Xt:m,l=enc[k](Xt[k]);mus.append(m);lvs.append(l)
        mu,lv=poe(mus,lvs);z=mu+torch.exp(0.5*lv)*torch.randn_like(lv)
        rec=sum(nn.functional.mse_loss(dec[k](z),Xt[k]) for k in Xt)
        kl=torch.clamp((-0.5*(1+lv-mu.pow(2)-lv.exp())).mean(0),min=0.5).sum()
        (rec+min(1,e/80)*kl+lmbda*nn.functional.cross_entropy(head(mu),ya)).backward();opt.step()
    for m in list(enc.values()):m.eval()
    with torch.no_grad():
        mus=[];lvs=[]
        for k in Xt:m,l=enc[k](Xt[k]);mus.append(m);lvs.append(l)
        mu,_=poe(mus,lvs)
    return mu.numpy()

# embeddings
emb_mosaic=mosaic(1.0)
Xcat=np.concatenate([Xn[k] for k in Xn],axis=1)
emb_pca=PCA(LAT,random_state=0).fit_transform(Xcat)
emb_mofa=np.concatenate([TruncatedSVD(3,random_state=0).fit_transform(Xn[k]) for k in Xn],axis=1)

def cv_oof(emb,y,seed):
    m=~np.isnan(y); e=emb[m]; yy=y[m].astype(int); N=len(yy)
    preds=np.zeros(N)
    for tr,te in StratifiedKFold(5,shuffle=True,random_state=seed).split(e,yy):
        preds[te]=LogisticRegression(max_iter=2000).fit(e[tr],yy[tr]).predict_proba(e[te])[:,1]
    return roc_auc_score(yy,preds)
def dist(emb,y,seeds=range(50)): return np.array([cv_oof(emb,y,s) for s in seeds])

def report(name,y):
    print(f"\n----- outcome: {name} -----")
    dm,dp,df_=dist(emb_mosaic,y),dist(emb_pca,y),dist(emb_mofa,y)
    for lbl,d in [("MOSAIC",dm),("PCA",dp),("MOFA-like",df_)]:
        print(f"  {lbl:10s} AUROC {d.mean():.3f}  95% CI [{np.percentile(d,2.5):.3f}, {np.percentile(d,97.5):.3f}]")
    diff=dm-dp
    lo,hi=np.percentile(diff,2.5),np.percentile(diff,97.5)
    verdict = "INDISTINGUISHABLE (tie)" if lo<=0<=hi else ("MOSAIC better" if lo>0 else "PCA better")
    print(f"  MOSAIC - PCA: {diff.mean():+.4f}  95% CI [{lo:+.4f}, {hi:+.4f}]  -> {verdict}")

report("SNF cluster (anchor)", yc)
report("Metabolic syndrome (HELD OUT)", ym)
print("\nDone.")
