#!/usr/bin/env python3
# =============================================================================
# mosaic_stage4_characterize.py
# Characterize the high-risk (E0) vs low-risk (E1) endotypes:
#   - clinical differences with proper significance tests
#   - differential molecular features (proteins, metabolites, genes)
#     via Mann-Whitney U + Benjamini-Hochberg FDR + effect sizes
#   - saves ranked feature tables (the drug-signature input) and a volcano plot
# Honest: reports how many survive multiple-testing correction, no inflation.
# =============================================================================
import os, numpy as np, pandas as pd
from scipy import stats

np.random.seed(0)
PROJ=os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
INT=os.path.join(PROJ,"data/processed/integrated"); OUT=os.path.join(PROJ,"results/mosaic")
os.makedirs(OUT,exist_ok=True)

def load(n): return pd.read_csv(os.path.join(INT,n),index_col=0)
rna,prot,meta=load("RNA_top_variable.csv"),load("proteomics_olink.csv"),load("metabolomics.csv")
endo=pd.read_csv(os.path.join(OUT,"definitive_endotypes.csv"))
endo.index=endo["COCOMO_ID"].astype(str)
ids=[str(i) for i in rna.index]; endo=endo.loc[ids]
g=endo["endotype"].values
# define E0 as the HIGHER-risk group by METS rate
mets_by=pd.Series(endo["METS"].values).groupby(g).mean()
high=int(mets_by.idxmax()); low=int(mets_by.idxmin())
gi=np.where(g==high,1,0)  # 1 = high-risk
print(f"High-risk endotype = E{high} (METS {mets_by[high]:.2f}, n={ (g==high).sum() }); "
      f"low-risk = E{low} (METS {mets_by[low]:.2f}, n={ (g==low).sum() })")

def bh_fdr(pv):
    pv=np.asarray(pv,float); n=len(pv); order=np.argsort(pv); ranked=pv[order]
    q=ranked*n/(np.arange(n)+1); q=np.minimum.accumulate(q[::-1])[::-1]
    out=np.empty(n); out[order]=np.clip(q,0,1); return out

def differential(df, layer):
    X=df.values.astype(float); feats=list(df.columns)
    a=X[gi==1]; b=X[gi==0]     # high vs low
    p=np.array([stats.mannwhitneyu(a[:,j],b[:,j],alternative="two-sided").pvalue
                for j in range(X.shape[1])])
    lfc=np.nanmean(a,0)-np.nanmean(b,0)   # mean diff (high - low)
    q=bh_fdr(p)
    out=pd.DataFrame({"feature":feats,"layer":layer,"mean_diff_high_minus_low":lfc,
                      "p":p,"q_FDR":q}).sort_values("q_FDR")
    nsig=(q<0.05).sum()
    print(f"  {layer:12s}: {X.shape[1]:5d} features | {nsig:4d} significant at FDR<0.05")
    return out

print("\n=== Differential molecular features (high-risk vs low-risk) ===")
res=pd.concat([differential(prot,"proteomics"),
               differential(meta,"metabolomics"),
               differential(rna,"transcriptomics")],ignore_index=True)
res.to_csv(os.path.join(OUT,"endotype_differential_features.csv"),index=False)

sig=res[res["q_FDR"]<0.05].sort_values("q_FDR")
print(f"\nTotal significant features (FDR<0.05): {len(sig)}")
print("Top 15 discriminating features:")
for _,r in sig.head(15).iterrows():
    direction="UP in high-risk" if r["mean_diff_high_minus_low"]>0 else "DOWN in high-risk"
    print(f"  {r['feature'][:32]:32s} {r['layer']:14s} q={r['q_FDR']:.2e}  {direction}")

# ---- clinical significance tests -------------------------------------------
print("\n=== Clinical differences (high vs low), tested ===")
clin_num=["BMI","Tgl","Hdl","Ldl","CD4","CD8","AGE","VAT","SAT"]
for c in clin_num:
    if c in endo:
        v=endo[c].values.astype(float)
        a=v[gi==1]; b=v[gi==0]; a=a[~np.isnan(a)]; b=b[~np.isnan(b)]
        if len(a)>3 and len(b)>3:
            p=stats.mannwhitneyu(a,b,alternative="two-sided").pvalue
            flag="*" if p<0.05 else " "
            print(f"  {c:5s} high {np.mean(a):7.2f}  low {np.mean(b):7.2f}  p={p:.3g} {flag}")
# METS as categorical
mets=endo["METS"].values.astype(float); m=~np.isnan(mets)
tab=pd.crosstab(gi[m], mets[m].astype(int))
chi=stats.chi2_contingency(tab)[1]
print(f"  METS  chi-square p={chi:.3g}  {'*' if chi<0.05 else ''}")

# ---- volcano plot -----------------------------------------------------------
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.figure(figsize=(8,6))
    x=res["mean_diff_high_minus_low"]; y=-np.log10(res["q_FDR"]+1e-300)
    colors=np.where(res["q_FDR"]<0.05,"crimson","lightgray")
    plt.scatter(x,y,c=colors,s=12,alpha=0.6)
    plt.axhline(-np.log10(0.05),ls="--",c="black",lw=0.8)
    plt.xlabel("mean difference (high-risk minus low-risk)")
    plt.ylabel("-log10 FDR q-value")
    plt.title("Molecular features distinguishing high- vs low-risk HIV endotype")
    plt.tight_layout(); plt.savefig(os.path.join(OUT,"endotype_volcano.png"),dpi=150)
    print("\nSaved volcano plot: results/mosaic/endotype_volcano.png")
except Exception as e:
    print("plot skipped:",e)

print("\nStage 4 complete. Saved endotype_differential_features.csv")
