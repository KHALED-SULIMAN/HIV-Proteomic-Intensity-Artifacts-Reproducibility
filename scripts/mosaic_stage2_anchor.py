#!/usr/bin/env python3
# =============================================================================
# mosaic_stage2_anchor.py
# MOSAIC-HIV Stage 2: the ANCHOR and the lambda-sweep (heart of the novelty).
#
#   anchor:      SNF cluster (published benchmark), weighted by lambda
#   validation:  METS and lipids are HELD OUT (never trained on)
#   novelty:     sweep lambda -> stability-vs-clinical-relevance frontier
#   benchmark:   compare MOSAIC latent vs PCA and per-layer-SVD (MOFA-like)
#                on held-out outcome separation (cross-validated AUROC)
#
# Alignment guardrail: asserts identical patient order across all layers.
# =============================================================================
import os, numpy as np, pandas as pd, torch, torch.nn as nn
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score

torch.manual_seed(0); np.random.seed(0)

PROJ = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
INT  = os.path.join(PROJ, "data/processed/integrated")
OUT  = os.path.join(PROJ, "results/mosaic"); os.makedirs(OUT, exist_ok=True)

LATENT_DIM=8; HIDDEN=64; DROPOUT=0.3; WD=1e-4
FREE_BITS=0.5; ANNEAL_EP=80; EPOCHS=300; LR=1e-3
LAMBDAS=[0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]   # the anchor sweep

# ---- load aligned layers + outcomes ----------------------------------------
def load(n): return pd.read_csv(os.path.join(INT,n), index_col=0)
rna, prot, meta = load("RNA_top_variable.csv"), load("proteomics_olink.csv"), load("metabolomics.csv")
outc = pd.read_csv(os.path.join(INT,"outcomes.csv"))
outc.index = outc["COCOMO_ID"].astype(str)

# GUARDRAIL: identical patient order everywhere
ids = [str(i) for i in rna.index]
assert [str(i) for i in prot.index]==ids and [str(i) for i in meta.index]==ids, "omics order mismatch!"
outc = outc.loc[ids]
assert list(outc.index)==ids, "outcomes not aligned to omics order!"
print(f"Alignment OK: {len(ids)} patients, identical order across all layers.")

def prep(df, cap=1000):
    X=df.values.astype(np.float32); mu,sd=np.nanmean(X,0),np.nanstd(X,0)+1e-8
    X=np.nan_to_num((X-mu)/sd,nan=0.0)
    if cap and X.shape[1]>cap:
        v=X.var(0); X=X[:,np.argsort(v)[::-1][:cap]]
    return torch.tensor(X)
X={"rna":prep(rna),"prot":prep(prot),"meta":prep(meta)}
dims={k:v.shape[1] for k,v in X.items()}

# targets
y_cluster = (outc["cluster"].values.astype(float)).astype(int)      # anchor
y_cluster = y_cluster - y_cluster.min()                              # -> 0/1
y_mets    = outc["METS"].values.astype(float)                       # held-out
# a lipid-based held-out label: high triglycerides (top tertile) as a clean binary
tgl = outc["Tgl"].values.astype(float)
y_tgl = (tgl >= np.nanquantile(tgl, 0.5)).astype(int)               # held-out

def clean_xy(emb, y):
    m=~np.isnan(y); return emb[m], y[m].astype(int)

# ---- model ------------------------------------------------------------------
class Enc(nn.Module):
    def __init__(s,di): super().__init__(); s.n=nn.Sequential(nn.Linear(di,HIDDEN),nn.ReLU(),nn.Dropout(DROPOUT)); s.mu=nn.Linear(HIDDEN,LATENT_DIM); s.lv=nn.Linear(HIDDEN,LATENT_DIM)
    def forward(s,x): h=s.n(x); return s.mu(h),s.lv(h)
class Dec(nn.Module):
    def __init__(s,do): super().__init__(); s.n=nn.Sequential(nn.Linear(LATENT_DIM,HIDDEN),nn.ReLU(),nn.Dropout(DROPOUT),nn.Linear(HIDDEN,do))
    def forward(s,z): return s.n(z)
def poe(mus,lvs):
    mus=torch.stack(mus,0);lvs=torch.stack(lvs,0);T=torch.exp(-lvs)
    Ts=torch.ones_like(T[0])+T.sum(0); return (mus*T).sum(0)/Ts, torch.log(1.0/Ts)

def train_mosaic(lmbda):
    enc=nn.ModuleDict({k:Enc(dims[k]) for k in dims}); dec=nn.ModuleDict({k:Dec(dims[k]) for k in dims})
    head=nn.Linear(LATENT_DIM,2)
    P=list(enc.parameters())+list(dec.parameters())+list(head.parameters())
    opt=torch.optim.Adam(P,lr=LR,weight_decay=WD)
    ya=torch.tensor(y_cluster)
    for ep in range(EPOCHS):
        for m in list(enc.values())+list(dec.values()): m.train()
        opt.zero_grad(); mus=[];lvs=[]
        for k in X: m,l=enc[k](X[k]); mus.append(m);lvs.append(l)
        mu,lv=poe(mus,lvs); z=mu+torch.exp(0.5*lv)*torch.randn_like(lv)
        rec=sum(nn.functional.mse_loss(dec[k](z),X[k]) for k in X)
        kl=torch.clamp((-0.5*(1+lv-mu.pow(2)-lv.exp())).mean(0),min=FREE_BITS).sum()
        beta=min(1.0,ep/ANNEAL_EP)
        anchor=nn.functional.cross_entropy(head(mu),ya)
        (rec+beta*kl+lmbda*anchor).backward(); opt.step()
    for m in list(enc.values()): m.eval()
    with torch.no_grad():
        mus=[];lvs=[]
        for k in X: m,l=enc[k](X[k]); mus.append(m);lvs.append(l)
        mu,lv=poe(mus,lvs)
    return mu.numpy(), torch.exp(lv).mean(1).numpy()

def auroc(emb,y):
    e,yy=clean_xy(emb,y)
    if len(np.unique(yy))<2: return np.nan
    return cross_val_score(LogisticRegression(max_iter=2000), e, yy, cv=5, scoring="roc_auc").mean()

# ---- the lambda sweep -------------------------------------------------------
print("\n===== ANCHOR SWEEP (cross-validated AUROC) =====")
print(f"{'lambda':>7} | {'cluster(anchor)':>15} | {'METS(heldout)':>13} | {'highTgl(heldout)':>16}")
rows=[]
best=None
for lmb in LAMBDAS:
    emb,unc=train_mosaic(lmb)
    a=auroc(emb,y_cluster); m=auroc(emb,y_mets); t=auroc(emb,y_tgl)
    print(f"{lmb:7.2f} | {a:15.3f} | {m:13.3f} | {t:16.3f}")
    rows.append({"lambda":lmb,"AUROC_cluster":a,"AUROC_METS":m,"AUROC_highTgl":t})
    # keep embedding at a moderate lambda for downstream use
    if lmb==1.0:
        best=(emb,unc)
pd.DataFrame(rows).to_csv(os.path.join(OUT,"stage2_lambda_sweep.csv"),index=False)

# ---- baselines on the same held-out outcomes -------------------------------
Xcat=np.concatenate([X[k].numpy() for k in X],axis=1)
pca=PCA(n_components=LATENT_DIM,random_state=0).fit_transform(Xcat)
# MOFA-like: per-layer SVD factors concatenated
svd_parts=[TruncatedSVD(n_components=3,random_state=0).fit_transform(X[k].numpy()) for k in X]
mofa=np.concatenate(svd_parts,axis=1)

print("\n===== BASELINES (same held-out outcomes) =====")
print(f"{'method':>12} | {'cluster':>8} | {'METS':>6} | {'highTgl':>8}")
for name,emb in [("PCA",pca),("MOFA-like",mofa)]:
    print(f"{name:>12} | {auroc(emb,y_cluster):8.3f} | {auroc(emb,y_mets):6.3f} | {auroc(emb,y_tgl):8.3f}")

# ---- save the lambda=1 embedding + uncertainty for downstream stages -------
if best is not None:
    emb,unc=best
    pd.DataFrame(emb,index=ids,columns=[f"z{i+1}" for i in range(LATENT_DIM)]).to_csv(
        os.path.join(OUT,"mosaic_stage2_embedding.csv"))
    pd.DataFrame({"latent_var_mean":unc},index=ids).to_csv(
        os.path.join(OUT,"mosaic_stage2_uncertainty.csv"))

# ---- plot the frontier ------------------------------------------------------
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    d=pd.DataFrame(rows)
    plt.figure(figsize=(8,5))
    plt.plot(d["lambda"],d["AUROC_cluster"],"o-",label="cluster (anchor)")
    plt.plot(d["lambda"],d["AUROC_METS"],"s-",label="METS (held out)")
    plt.plot(d["lambda"],d["AUROC_highTgl"],"^-",label="high triglycerides (held out)")
    plt.xlabel("anchor weight  lambda"); plt.ylabel("cross-validated AUROC")
    plt.title("MOSAIC anchor-sweep: stability vs clinical relevance frontier")
    plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
    plt.savefig(os.path.join(OUT,"stage2_lambda_frontier.png"),dpi=150)
    print("\nSaved frontier plot: results/mosaic/stage2_lambda_frontier.png")
except Exception as e:
    print("plot skipped:",e)

print("\nStage 2 complete.")
