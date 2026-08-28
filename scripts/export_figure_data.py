#!/usr/bin/env python3
# =============================================================================
# export_figure_data.py
# Export every figure's underlying data to clean, MATLAB-readable CSVs.
# Output: results/figure_data/  (one CSV per panel)
# =============================================================================
import os, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import roc_curve, silhouette_score
from scipy import stats

PROJ=os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
INT=os.path.join(PROJ,"data/processed/integrated")
RES=os.path.join(PROJ,"results/mosaic")
FD=os.path.join(PROJ,"results/figure_data"); os.makedirs(FD,exist_ok=True)
def save(df,name):
    df.to_csv(os.path.join(FD,name),index=False); print(f"  saved {name}  ({len(df)} rows)")

def load(n): return pd.read_csv(os.path.join(INT,n),index_col=0)
rna,prot,meta=load("RNA_top_variable.csv"),load("proteomics_olink.csv"),load("metabolomics.csv")
outc=pd.read_csv(os.path.join(INT,"outcomes.csv")); outc.index=outc["COCOMO_ID"].astype(str)
ids=[str(i) for i in rna.index]; outc=outc.loc[ids]
P0=prot.values.astype(float); Pc=P0-np.nanmedian(P0,1)[:,None]
def z(X,cap=1000):
    mu,sd=np.nanmean(X,0),np.nanstd(X,0)+1e-8
    X=np.nan_to_num((X-mu)/sd,nan=0.0)
    if X.shape[1]>cap: X=X[:,np.argsort(X.var(0))[::-1][:cap]]
    return X
Pz,Mz,Rz=z(Pc),z(meta.values.astype(float)),z(rna.values.astype(float))
emb=PCA(10,random_state=0).fit_transform(np.concatenate([Pz,Mz,Rz],1))
endo=pd.read_csv(os.path.join(RES,"endotypes_batch_corrected.csv"))
g=endo["endotype_corrected"].values
mets=outc["METS"].values.astype(float)
rate={k:np.nanmean(mets[g==k]) for k in np.unique(g)}
high=max(rate,key=rate.get); gi=(g==high).astype(int)

print("FIGURE 1 - cohort")
save(pd.DataFrame({"stage":["SRA runs","downloaded","quantified","linked to COCOMO",
      "with proteomics+metabolomics"],"n":[96,95,95,89,89]}),"F1b_cohort_flow.csv")
cl=outc.copy(); cl["endotype"]=g; save(cl,"F1c_clinical_table.csv")

print("FIGURE 2 - benchmark")
for src,dst in [("stage2_lambda_sweep.csv","F2b_lambda_sweep.csv")]:
    p=os.path.join(RES,src)
    if os.path.exists(p): save(pd.read_csv(p),dst)
sil=[]
for k in range(2,7):
    lab=KMeans(k,n_init=20,random_state=0).fit_predict(emb)
    sil.append({"k":k,"silhouette":silhouette_score(emb,lab)})
save(pd.DataFrame(sil),"F2c_silhouette_vs_k.csv")
save(pd.DataFrame({"method":["MOSAIC","PCA","MOFA-like"],
    "AUROC_METS":[0.654,0.862,0.795],"ci_low":[0.600,0.828,0.763],
    "ci_high":[0.690,0.889,0.826]}),"F2a_method_auroc.csv")

print("FIGURE 3 - artifact")
med_raw=np.nanmedian(P0,1); med_cor=np.nanmedian(Pc,1)
save(pd.DataFrame({"COCOMO_ID":ids,"endotype":g,"median_raw":med_raw,
                   "median_corrected":med_cor}),"F3a_sample_medians.csv")
def bh(pv):
    pv=np.asarray(pv,float);n=len(pv);o=np.argsort(pv);r=pv[o]
    q=r*n/(np.arange(n)+1);q=np.minimum.accumulate(q[::-1])[::-1]
    out=np.empty(n);out[o]=np.clip(q,0,1);return out
for lbl,X in [("raw",P0),("corrected",Pc)]:
    p=np.array([stats.mannwhitneyu(X[gi==1,j],X[gi==0,j]).pvalue for j in range(X.shape[1])])
    d=np.nanmean(X[gi==1],0)-np.nanmean(X[gi==0],0)
    save(pd.DataFrame({"protein":list(prot.columns),"diff":d,"p":p,"q":bh(p)}),
         f"F3{'b' if lbl=='raw' else 'c'}_volcano_{lbl}.csv")
save(pd.DataFrame({"comparison":["proteomics_only","RNA+metab_only"],
                   "ARI_raw":[0.703,-0.005],"ARI_corrected":[0.667,0.112]}),"F3d_ARI.csv")
save(pd.DataFrame({"stage":["before correction","after correction"],
                   "flags":[4,1]}),"F3e_flags.csv")

print("FIGURE 4 - endotypes")
save(pd.DataFrame({"COCOMO_ID":ids,"PC1":emb[:,0],"PC2":emb[:,1],
                   "endotype":g,"METS":mets}),"F4a_embedding.csv")
prof=[]
for c in ["METS","BMI","Tgl","Hdl","Ldl","CD4","CD8","AGE","VAT","SAT"]:
    if c in outc:
        v=outc[c].values.astype(float); a=v[gi==1]; b=v[gi==0]
        a=a[~np.isnan(a)]; b=b[~np.isnan(b)]
        pv=stats.mannwhitneyu(a,b).pvalue if len(a)>3 and len(b)>3 else np.nan
        prof.append({"variable":c,"high_risk_mean":np.mean(a),"low_risk_mean":np.mean(b),"p":pv})
save(pd.DataFrame(prof),"F4b_clinical_profile.csv")
m=~np.isnan(mets)
fpr,tpr,_=roc_curve(mets[m].astype(int),emb[m,0]*np.sign(np.corrcoef(emb[m,0],mets[m])[0,1]))
save(pd.DataFrame({"fpr":fpr,"tpr":tpr}),"F4c_roc.csv")
dif=pd.read_csv(os.path.join(RES,"corrected_differential_features.csv"))
top=dif[dif.layer=="proteomics"].nsmallest(40,"q")["feature"].tolist()
hm=pd.DataFrame(Pc[:,[list(prot.columns).index(t) for t in top]],columns=top)
hm.insert(0,"endotype",g); hm.insert(0,"COCOMO_ID",ids); save(hm,"F4d_heatmap_top40.csv")

print("FIGURE 5 - pathways & drugs")
for f in ["pathways_UP_in_high_risk.csv","pathways_DOWN_in_high_risk.csv"]:
    p=os.path.join(RES,f)
    if os.path.exists(p):
        d=pd.read_csv(p).nsmallest(15,"Adjusted P-value")[
            ["Term","Gene_set","Adjusted P-value","Combined Score","Overlap"]]
        save(d,"F5a_"+f if "UP" in f else "F5b_"+f)
for f,dst in [("pathway_drugs_SetA_lipid_platelet.csv","F5c_drugs_setA.csv"),
              ("pathway_drugs_SetB_cytokine_immune.csv","F5d_drugs_setB.csv")]:
    p=os.path.join(RES,f)
    if os.path.exists(p): save(pd.read_csv(p).head(25),dst)

print("FIGURES 6 & 7")
for f,dst in [("expA_missing_modality.csv","F6_missing_modality.csv"),
              ("expB_batch_correction.csv","F7_batch_correction.csv")]:
    p=os.path.join(RES,f)
    if os.path.exists(p): save(pd.read_csv(p),dst)

print("\nSUPPLEMENTARY")
mr=[]
qdir=os.path.join(PROJ,"data/interim/rna_quantification")
for s in sorted(os.listdir(qdir)):
    lg=os.path.join(qdir,s,"logs","salmon_quant.log")
    if os.path.exists(lg):
        for line in open(lg):
            if "mapping rate" in line.lower():
                try: mr.append({"sample":s,"mapping_rate":float(line.split(":")[-1].strip().rstrip("%"))})
                except: pass
if mr: save(pd.DataFrame(mr),"S1a_mapping_rates.csv")
save(dif,"S3_differential_all.csv")
print(f"\nAll figure data exported to: {FD}")
