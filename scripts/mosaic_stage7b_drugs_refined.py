#!/usr/bin/env python3
# =============================================================================
# mosaic_stage7b_drugs_refined.py
# ART-aware drug repurposing - REFINED.
# Fixes over v1: tight signature (top-N by q), strict FDR, clean drug-name
# parsing (cell-line tokens removed), and a cardiometabolic/anti-inflammatory
# mechanism-plausibility filter driven by the observed pathways
# (lipid & atherosclerosis, platelet activation, cytokine signalling).
#
# *** COMPUTATIONAL HYPOTHESES ONLY - NOT CLINICAL ADVICE ***
# *** ART flags are a screening aid; verify with Liverpool HIV Drug        ***
# *** Interactions + clinical pharmacologist before any manuscript claim.  ***
# =============================================================================
import os, re, pandas as pd, numpy as np, warnings
warnings.filterwarnings("ignore")

PROJ=os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT=os.path.join(PROJ,"results/mosaic")
TOP_N=150          # tighten signature
FDR=0.05           # strict

res=pd.read_csv(os.path.join(OUT,"corrected_differential_features.csv"))
p=res[(res.layer=="proteomics")&(res.q<0.05)].sort_values("q")
up=p[p.diff_high_minus_low>0]["feature"].head(TOP_N).tolist()
dn=p[p.diff_high_minus_low<0]["feature"].head(TOP_N).tolist()
print(f"Tightened signature: {len(up)} UP, {len(dn)} DOWN (top {TOP_N} by q)")

import gseapy as gp
CELL={"mcf10a","ha1e","a375","mcf7","npc","mdamb231","pc3","vcap","hepg2","ht29",
      "a549","haie","hcc515","hek293","hs578t","skbr3","jurkat","thp1","huvec",
      "neu","fibrntpc","nkdba","ashp","24h","6h","3h","um","nm"}
PREFIX=re.compile(r"^(ljp|cpc|cpd|dos|bra|erg|hog|pcl|ctpc|rad|mucl|kdc|nmh|xpr)\d*$",re.I)

def clean_name(term):
    toks=re.split(r"[\s\-_:,\.]+",str(term))
    out=[]
    for t in toks:
        tl=t.lower()
        if not tl.isalpha(): continue
        if len(tl)<4: continue
        if tl in CELL or PREFIX.match(tl): continue
        out.append(tl)
    return out[-1] if out else None

def enrich(genes,lib):
    try:
        r=gp.enrichr(gene_list=genes,gene_sets=[lib],organism="human",
                     outdir=None,no_plot=True).results
        return r[r["Adjusted P-value"]<FDR].copy()
    except Exception as e:
        print(f"  {lib} failed: {e}"); return pd.DataFrame()

print("\nQuerying LINCS L1000 (reversal direction, FDR<0.05)...")
r1=enrich(up,"LINCS_L1000_Chem_Pert_down")   # drug DOWNs what's UP
r2=enrich(dn,"LINCS_L1000_Chem_Pert_up")     # drug UPs what's DOWN
print(f"  raw hits: {len(r1)} / {len(r2)}")

def tally(df,col):
    if df.empty: return pd.DataFrame(columns=["drug",col])
    d=df.copy(); d["drug"]=d["Term"].map(clean_name); d=d.dropna(subset=["drug"])
    return d.groupby("drug")["Adjusted P-value"].min().reset_index().rename(columns={"Adjusted P-value":col})

m=pd.merge(tally(r1,"q_downUP"),tally(r2,"q_upDOWN"),on="drug",how="outer")
m["n_dir"]=m[["q_downUP","q_upDOWN"]].notna().sum(axis=1)
m["best_q"]=m[["q_downUP","q_upDOWN"]].min(axis=1)
m=m.sort_values(["n_dir","best_q"],ascending=[False,True])
print(f"  parsed unique drug names: {len(m)}")

# ---- mechanism plausibility for cardiometabolic / inflammatory biology -----
PLAUSIBLE={
 # lipid
 "atorvastatin","simvastatin","rosuvastatin","pravastatin","lovastatin","fluvastatin",
 "mevastatin","cerivastatin","fenofibrate","gemfibrozil","bezafibrate","clofibrate",
 "ezetimibe","niacin","nicotinic",
 # glucose / insulin sensitising
 "metformin","phenformin","pioglitazone","rosiglitazone","troglitazone","ciglitazone",
 "acarbose","sitagliptin","glipizide","glyburide","glimepiride",
 # anti-inflammatory
 "aspirin","salicylate","celecoxib","indomethacin","ibuprofen","diclofenac","naproxen",
 "sulindac","meloxicam","piroxicam","dexamethasone","prednisolone","hydrocortisone",
 "budesonide","colchicine","methotrexate","hydroxychloroquine","chloroquine",
 "sulfasalazine","mesalamine","auranofin",
 # antiplatelet / vascular
 "clopidogrel","ticlopidine","dipyridamole","cilostazol","pentoxifylline",
 "telmisartan","losartan","valsartan","irbesartan","candesartan","captopril",
 "enalapril","ramipril","amlodipine","nifedipine","propranolol","carvedilol",
 "spironolactone","eplerenone",
 # antioxidant / other repurposing
 "acetylcysteine","resveratrol","curcumin","quercetin","allopurinol","probucol",
 "doxycycline","minocycline","azithromycin","rapamycin","sirolimus","everolimus",
 "metronidazole","thalidomide","lenalidomide",
}
m["mechanistically_plausible"]=m["drug"].isin(PLAUSIBLE)

# ---- ART interaction screen (VERIFY) ---------------------------------------
ART={
 "simvastatin":"CONTRAINDICATED with boosted PI/cobicistat (myopathy)",
 "lovastatin":"CONTRAINDICATED with boosted PI/cobicistat (myopathy)",
 "atorvastatin":"dose-limit with boosted PI/cobicistat",
 "rosuvastatin":"dose adjust with some PI (e.g. atazanavir)",
 "rifampicin":"strong CYP3A4 inducer - avoid with PI/INSTI",
 "rifampin":"strong CYP3A4 inducer - avoid with PI/INSTI",
 "ketoconazole":"strong CYP3A4 inhibitor","itraconazole":"strong CYP3A4 inhibitor",
 "clarithromycin":"CYP3A4 inhibitor + QT","azithromycin":"lower CYP risk; monitor QT",
 "carbamazepine":"CYP3A4 inducer","phenytoin":"CYP3A4 inducer",
 "clopidogrel":"CYP2C19 activation may be reduced by some ART",
 "rapamycin":"CYP3A4 substrate - levels raised by boosters",
 "sirolimus":"CYP3A4 substrate - levels raised by boosters",
 "everolimus":"CYP3A4 substrate - levels raised by boosters",
 "colchicine":"CYP3A4/P-gp substrate - toxicity risk with boosters",
 "amlodipine":"CYP3A4 substrate - monitor",
 "spironolactone":"monitor K+ with TDF/other nephrotoxicity",
}
m["ART_flag"]=m["drug"].map(lambda d: ART.get(d,""))
m["ART_screen"]=np.where(m["ART_flag"]=="","no flag in this screen","REVIEW REQUIRED")
m.to_csv(os.path.join(OUT,"drug_candidates_refined.csv"),index=False)

pl=m[m["mechanistically_plausible"]].sort_values(["n_dir","best_q"],ascending=[False,True])
print("\n"+"="*70)
print("MECHANISTICALLY PLAUSIBLE CANDIDATES (cardiometabolic / inflammatory)")
print("="*70)
if len(pl)==0:
    print("  NONE. Signature reversal did not recover plausible cardiometabolic")
    print("  agents. This is an honest negative result - report it as such.")
else:
    for _,r in pl.head(25).iterrows():
        print(f"  {r['drug'][:22]:22s} dirs={int(r['n_dir'])} q={r['best_q']:.2e}  "
              f"ART: {r['ART_flag'] or 'no flag in screen'}")

print(f"\nTotal parsed drugs: {len(m)} | plausible: {len(pl)} | "
      f"plausible with ART flag: {int((pl['ART_flag']!='').sum()) if len(pl) else 0}")
print("\n"+"!"*70)
print("COMPUTATIONAL HYPOTHESES ONLY - NOT CLINICAL RECOMMENDATIONS.")
print("Verify every ART flag against the Liverpool HIV Drug Interactions")
print("database and have a clinical pharmacologist review before publication.")
print("!"*70)
print("\nSaved: drug_candidates_refined.csv")
