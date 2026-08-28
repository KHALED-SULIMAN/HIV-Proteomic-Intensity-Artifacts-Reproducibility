#!/usr/bin/env python3
# =============================================================================
# mosaic_stage1_core.py  (v2: KL annealing + free-bits, train/val split)
# MOSAIC-HIV, Stage 1: core multi-view VAE with product-of-experts fusion.
# Fixes posterior collapse (KL->0) seen in v1.
# =============================================================================
import os, numpy as np, pandas as pd, torch, torch.nn as nn

torch.manual_seed(0); np.random.seed(0)

PROJ = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
INT  = os.path.join(PROJ, "data/processed/integrated")
OUT  = os.path.join(PROJ, "results/mosaic"); os.makedirs(OUT, exist_ok=True)

LATENT_DIM = 16
HIDDEN     = 256
EPOCHS     = 400
LR         = 1e-3
FREE_BITS  = 0.5     # minimum nats each latent dim must carry (anti-collapse)
ANNEAL_EP  = 100     # ramp beta from 0 -> 1 over this many epochs
VAL_FRAC   = 0.2     # hold out 20% of patients to watch for overfitting

def load(name):
    return pd.read_csv(os.path.join(INT, name), index_col=0)
rna  = load("RNA_top_variable.csv")
prot = load("proteomics_olink.csv")
meta = load("metabolomics.csv")
assert list(rna.index) == list(prot.index) == list(meta.index), "patient order mismatch!"
patients = list(rna.index)
print(f"Loaded {len(patients)} patients | RNA {rna.shape[1]}, "
      f"PROT {prot.shape[1]}, META {meta.shape[1]} features")

def prep(df):
    X = df.values.astype(np.float32)
    mu, sd = np.nanmean(X,0), np.nanstd(X,0)+1e-8
    X = (X-mu)/sd
    return torch.tensor(np.nan_to_num(X, nan=0.0))
layers = {"rna": prep(rna), "prot": prep(prot), "meta": prep(meta)}
dims   = {k: v.shape[1] for k,v in layers.items()}

# train/val split (by patient index)
n = len(patients); idx = np.random.permutation(n)
n_val = int(round(VAL_FRAC*n)); val_idx = idx[:n_val]; tr_idx = idx[n_val:]
def subset(t, ix): return t[ix]

class Encoder(nn.Module):
    def __init__(self, d_in, d_lat):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(d_in, HIDDEN), nn.ReLU(),
                                 nn.Linear(HIDDEN, HIDDEN), nn.ReLU())
        self.mu     = nn.Linear(HIDDEN, d_lat)
        self.logvar = nn.Linear(HIDDEN, d_lat)
    def forward(self, x):
        h = self.net(x); return self.mu(h), self.logvar(h)

class Decoder(nn.Module):
    def __init__(self, d_lat, d_out):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(d_lat, HIDDEN), nn.ReLU(),
                                 nn.Linear(HIDDEN, HIDDEN), nn.ReLU(),
                                 nn.Linear(HIDDEN, d_out))
    def forward(self, z): return self.net(z)

def poe(mus, logvars):
    mus     = torch.stack(mus, dim=0)
    logvars = torch.stack(logvars, dim=0)
    T = torch.exp(-logvars)
    prior_T = torch.ones_like(T[0])
    T_sum  = prior_T + T.sum(dim=0)
    mu_sum = (mus*T).sum(dim=0)
    mu     = mu_sum / T_sum
    logvar = torch.log(1.0/T_sum)
    return mu, logvar

class MOSAIC(nn.Module):
    def __init__(self, dims, d_lat):
        super().__init__()
        self.enc = nn.ModuleDict({k: Encoder(d, d_lat) for k,d in dims.items()})
        self.dec = nn.ModuleDict({k: Decoder(d_lat, d) for k,d in dims.items()})
    def forward(self, batch, mask):
        mus, lvs = [], []
        for k in batch:
            if mask[k]:
                mu, lv = self.enc[k](batch[k]); mus.append(mu); lvs.append(lv)
        mu, logvar = poe(mus, lvs)
        std = torch.exp(0.5*logvar); z = mu + std*torch.randn_like(std)
        recon = {k: self.dec[k](z) for k in batch}
        return recon, mu, logvar, z

device = "cuda" if torch.cuda.is_available() else "cpu"
model  = MOSAIC(dims, LATENT_DIM).to(device)
opt    = torch.optim.Adam(model.parameters(), lr=LR)
Xfull  = {k: v.to(device) for k,v in layers.items()}
Xtr    = {k: subset(v, tr_idx) for k,v in Xfull.items()}
Xval   = {k: subset(v, val_idx) for k,v in Xfull.items()}
mask   = {k: True for k in Xfull}

def losses(Xb):
    recon, mu, logvar, z = model(Xb, mask)
    rec = sum(nn.functional.mse_loss(recon[k], Xb[k]) for k in Xb)
    kl_dim = (-0.5*(1+logvar-mu.pow(2)-logvar.exp())).mean(0)   # per-dim
    kl = torch.clamp(kl_dim, min=FREE_BITS).sum()               # free bits
    return rec, kl

print(f"Training on {device} | train {len(tr_idx)}  val {len(val_idx)}")
for ep in range(1, EPOCHS+1):
    model.train(); opt.zero_grad()
    rec, kl = losses(Xtr)
    beta = min(1.0, ep/ANNEAL_EP)
    (rec + beta*kl).backward(); opt.step()
    if ep % 50 == 0 or ep == 1:
        model.eval()
        with torch.no_grad(): vrec, vkl = losses(Xval)
        print(f"  ep {ep:4d} | train rec {rec.item():.3f} kl {kl.item():.3f} "
              f"| val rec {vrec.item():.3f} | beta {beta:.2f}")

# final embedding on ALL patients
model.eval()
with torch.no_grad():
    _, mu, logvar, _ = model(Xfull, mask)
emb = pd.DataFrame(mu.cpu().numpy(), index=patients,
                   columns=[f"z{i+1}" for i in range(LATENT_DIM)])
emb.to_csv(os.path.join(OUT, "mosaic_stage1_embedding.csv"))
# also save per-patient latent uncertainty (mean posterior variance) for later
unc = pd.DataFrame({"latent_var_mean": torch.exp(logvar).mean(1).cpu().numpy()},
                   index=patients)
unc.to_csv(os.path.join(OUT, "mosaic_stage1_uncertainty.csv"))
print(f"\nSaved embedding: {emb.shape[0]} x {emb.shape[1]}  (+ uncertainty)")

try:
    import umap, matplotlib
    matplotlib.use("Agg"); import matplotlib.pyplot as plt
    u = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=0).fit_transform(mu.cpu().numpy())
    plt.figure(figsize=(7,6)); plt.scatter(u[:,0], u[:,1], s=40)
    plt.title("MOSAIC Stage 1 - patient latent embedding (UMAP)")
    plt.xlabel("UMAP-1"); plt.ylabel("UMAP-2"); plt.tight_layout()
    p = os.path.join(OUT, "mosaic_stage1_umap.png"); plt.savefig(p, dpi=150)
    print(f"Saved UMAP: {p}")
except Exception as e:
    print("UMAP step skipped:", e)

print("\nStage 1 (v2) complete.")
