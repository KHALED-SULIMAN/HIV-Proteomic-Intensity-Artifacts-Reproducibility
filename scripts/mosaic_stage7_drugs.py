#!/usr/bin/env python3
# =============================================================================
# mosaic_stage7_drugs.py
# ART-aware drug repurposing by LINCS L1000 signature reversal.
#
# Logic: a drug that REVERSES the high-risk endotype signature should
#        DOWN-regulate the proteins UP in high-risk, and
#        UP-regulate the proteins DOWN in high-risk.
#        Drugs hit in BOTH directions = strongest reversal candidates.
#
# Then flags antiretroviral (ART) interaction risk.
#
# *** ALL OUTPUT IS COMPUTATIONAL HYPOTHESIS, NOT CLINICAL ADVICE. ***
# *** The ART interaction flags MUST be verified by a pharmacist/       ***
# *** clinical pharmacologist before any manuscript claim.              ***
# =============================================================================
import os, re, pandas as pd, warnings
warnings.filterwarnings("ignore")

PROJ=os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT=os.path.join(PROJ,"results/mosaic")

up=[l.strip() for l in open(os.path.join(OUT,"proteins_up_highrisk.txt")) if l.strip()]
dn=[l.strip() for l in open(os.path.join(OUT,"proteins_down_highrisk.txt")) if l.strip()]
print(f"Signature: {len(up)} UP, {len(dn)} DOWN in high-risk endotype")

import gseapy as gp
LIB_UP="LINCS_L1000_Chem_Pert_up"; LIB_DN="LINCS_L1000_Chem_Pert_down"

def enrich(genes, lib):
    try:
        r=gp.enrichr(gene_list=genes, gene_sets=[lib], organism="human",
                     outdir=None, no_plot=True).results
        return r[r["Adjusted P-value"]<0.25].copy()
    except Exception as e:
        print(f"  query failed for {lib}: {e}"); return pd.DataFrame()

def drug_name(term):
    # LINCS terms look like: "CPC006_HA1E_6H-BRD-K12345-001-01-1-10" or contain a drug name
    t=re.split(r"[-_]", str(term))
    cands=[x for x in t if x.isalpha() and len(x)>3]
    return cands[-1].lower() if cands else str(term)

print("\nQuerying LINCS L1000 via Enrichr (reversal direction)...")
rev1=enrich(up, LIB_DN)   # drugs that DOWN what's UP
rev2=enrich(dn, LIB_UP)   # drugs that UP what's DOWN
print(f"  hits: {len(rev1)} (down-regulate UP set), {len(rev2)} (up-regulate DOWN set)")

def tally(df,label):
    if df.empty: return pd.DataFrame(columns=["drug",label])
    d=df.copy(); d["drug"]=d["Term"].map(drug_name)
    g=d.groupby("drug")["Adjusted P-value"].min().reset_index()
    g.columns=["drug",label]; return g

t1,t2=tally(rev1,"q_down_of_UP"), tally(rev2,"q_up_of_DOWN")
merged=pd.merge(t1,t2,on="drug",how="outer")
merged["n_directions"]=merged[["q_down_of_UP","q_up_of_DOWN"]].notna().sum(axis=1)
merged["best_q"]=merged[["q_down_of_UP","q_up_of_DOWN"]].min(axis=1)
merged=merged.sort_values(["n_directions","best_q"],ascending=[False,True])

# ---- ART interaction flags (VERIFY BEFORE USE) -----------------------------
# Common mechanisms: CYP3A4 inhibition/induction, P-gp, QT prolongation,
# nephrotoxicity (with TDF), myopathy (with statins + boosters).
ART_RISK = {
 "rifampicin":"strong CYP3A4 inducer - contraindicated with PI/INSTI",
 "rifampin":"strong CYP3A4 inducer - contraindicated with PI/INSTI",
 "ketoconazole":"strong CYP3A4 inhibitor",
 "itraconazole":"strong CYP3A4 inhibitor",
 "clarithromycin":"CYP3A4 inhibitor + QT risk",
 "simvastatin":"contraindicated with boosted PI (myopathy)",
 "lovastatin":"contraindicated with boosted PI (myopathy)",
 "carbamazepine":"CYP3A4 inducer",
 "phenytoin":"CYP3A4 inducer",
 "amiodarone":"CYP3A4 + QT risk",
 "quinidine":"CYP3A4 + QT risk",
 "midazolam":"CYP3A4 substrate - accumulation risk",
 "cisapride":"QT risk - contraindicated",
 "st-john":"CYP3A4 inducer - contraindicated",
}
def art_flag(d):
    for k,v in ART_RISK.items():
        if k in d: return v
    return ""
merged["ART_interaction_flag"]=merged["drug"].map(art_flag)
merged["ART_safe_candidate"]=merged["ART_interaction_flag"]==""

merged.to_csv(os.path.join(OUT,"drug_repurposing_candidates.csv"),index=False)

print("\n=== TOP REVERSAL CANDIDATES (both directions first) ===")
top=merged[merged["n_directions"]==2].head(20)
if len(top)==0: top=merged.head(20)
for _,r in top.iterrows():
    flag=r["ART_interaction_flag"] or "no known ART interaction in this list"
    print(f"  {r['drug'][:28]:28s} dirs={int(r['n_directions'])} q={r['best_q']:.2e}  [{flag}]")

n_safe=int(merged["ART_safe_candidate"].sum())
print(f"\nTotal candidates: {len(merged)} | flagged for ART interaction: {len(merged)-n_safe}")
print("\n" + "!"*66)
print("COMPUTATIONAL HYPOTHESES ONLY - NOT CLINICAL RECOMMENDATIONS.")
print("ART interaction list is a starting screen compiled from general")
print("pharmacology and MUST be verified against a current interaction")
print("database (e.g. HIV Drug Interactions, Liverpool) and reviewed by a")
print("clinical pharmacologist before any manuscript claim.")
print("!"*66)
print("\nSaved: drug_repurposing_candidates.csv")
