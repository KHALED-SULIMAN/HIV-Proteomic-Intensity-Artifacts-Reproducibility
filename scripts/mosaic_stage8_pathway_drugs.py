#!/usr/bin/env python3
# =============================================================================
# mosaic_stage8_pathway_drugs.py
# PATHWAY-TARGETED drug discovery (replaces noisy signature reversal).
#
# Set A: lipid & atherosclerosis + platelet activation / hemostasis
# Set B: cytokine-cytokine receptor interaction + immunoregulatory signalling
#
# For each set: extract the differential proteins driving that pathway,
# query DGIdb for drugs targeting them, rank drugs by number of targets hit,
# then apply the ART interaction screen.
#
# *** COMPUTATIONAL HYPOTHESES ONLY - NOT CLINICAL ADVICE ***
# =============================================================================
import os, json, time, pandas as pd, numpy as np, warnings
warnings.filterwarnings("ignore")
import urllib.request, urllib.parse

PROJ=os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT=os.path.join(PROJ,"results/mosaic")

SET_A=["Lipid and atherosclerosis","Platelet Activation","Hemostasis",
       "Platelet activation","Lipid And Atherosclerosis"]
SET_B=["Cytokine-cytokine receptor interaction","Immunoregulatory Interactions",
       "Cytokine Signaling","Cytokine-Cytokine Receptor Interaction"]

def genes_from_pathways(fn, keys):
    p=os.path.join(OUT,fn)
    if not os.path.exists(p): return []
    df=pd.read_csv(p)
    g=set()
    for _,r in df.iterrows():
        term=str(r.get("Term",""))
        if any(k.lower() in term.lower() for k in keys):
            for x in str(r.get("Genes","")).split(";"):
                if x.strip(): g.add(x.strip().upper())
    return sorted(g)

A = sorted(set(genes_from_pathways("pathways_UP_in_high_risk.csv",SET_A) +
               genes_from_pathways("pathways_DOWN_in_high_risk.csv",SET_A)))
B = sorted(set(genes_from_pathways("pathways_UP_in_high_risk.csv",SET_B) +
               genes_from_pathways("pathways_DOWN_in_high_risk.csv",SET_B)))
print(f"Set A (lipid/atherosclerosis + platelet): {len(A)} target proteins")
print(f"  {', '.join(A[:20])}{' ...' if len(A)>20 else ''}")
print(f"Set B (cytokine/immune signalling):       {len(B)} target proteins")
print(f"  {', '.join(B[:20])}{' ...' if len(B)>20 else ''}")
if not A and not B:
    print("\nNo pathway gene lists found. Re-run stage 6 first."); raise SystemExit

def dgidb(genes, chunk=40):
    """Query DGIdb GraphQL for drug-gene interactions."""
    rows=[]
    url="https://dgidb.org/api/graphql"
    for i in range(0,len(genes),chunk):
        sub=genes[i:i+chunk]
        q={"query":"{genes(names:[%s]){nodes{name interactions{drug{name approved} "
                   "interactionScore interactionTypes{type}}}}}"
                   % ",".join(f'"{g}"' for g in sub)}
        try:
            req=urllib.request.Request(url,data=json.dumps(q).encode(),
                 headers={"Content-Type":"application/json"})
            d=json.loads(urllib.request.urlopen(req,timeout=60).read())
            for node in d.get("data",{}).get("genes",{}).get("nodes",[]):
                gname=node.get("name")
                for it in node.get("interactions",[]):
                    dr=it.get("drug",{}) or {}
                    rows.append({"gene":gname,"drug":(dr.get("name") or "").lower(),
                                 "approved":dr.get("approved"),
                                 "score":it.get("interactionScore")})
        except Exception as e:
            print(f"  DGIdb chunk {i//chunk+1} failed: {e}")
        time.sleep(1)
    return pd.DataFrame(rows)

ART={
 "simvastatin":"CONTRAINDICATED with boosted PI/cobicistat",
 "lovastatin":"CONTRAINDICATED with boosted PI/cobicistat",
 "atorvastatin":"dose-limit with boosted PI/cobicistat",
 "rosuvastatin":"dose adjust with atazanavir/others",
 "clopidogrel":"CYP2C19 activation may be affected",
 "colchicine":"CYP3A4/P-gp - toxicity risk with boosters",
 "sirolimus":"CYP3A4 substrate - raised by boosters",
 "everolimus":"CYP3A4 substrate - raised by boosters",
 "rifampicin":"strong CYP3A4 inducer - avoid",
 "ketoconazole":"strong CYP3A4 inhibitor",
 "clarithromycin":"CYP3A4 inhibitor + QT",
 "amiodarone":"CYP3A4 + QT",
 "carbamazepine":"CYP3A4 inducer",
 "warfarin":"narrow therapeutic index - monitor INR",
 "tacrolimus":"CYP3A4 substrate - raised by boosters",
 "cyclosporine":"CYP3A4 substrate - raised by boosters",
}

def summarise(df,label,targets):
    if df.empty:
        print(f"\n{label}: no DGIdb interactions retrieved."); return pd.DataFrame()
    d=df[df["drug"]!=""].copy()
    agg=d.groupby("drug").agg(n_targets=("gene","nunique"),
                              targets=("gene",lambda s:";".join(sorted(set(s)))),
                              approved=("approved","max"),
                              mean_score=("score","mean")).reset_index()
    agg["ART_flag"]=agg["drug"].map(lambda x: ART.get(x,""))
    agg=agg.sort_values(["n_targets","mean_score"],ascending=[False,False])
    agg.to_csv(os.path.join(OUT,f"pathway_drugs_{label}.csv"),index=False)
    print(f"\n{'='*70}\n{label}: {len(agg)} candidate drugs "
          f"(from {len(targets)} pathway proteins)\n{'='*70}")
    show=agg[agg["approved"]==True] if (agg["approved"]==True).any() else agg
    for _,r in show.head(20).iterrows():
        flag=r["ART_flag"] or "no flag in screen"
        print(f"  {r['drug'][:26]:26s} targets={int(r['n_targets']):2d}  "
              f"approved={r['approved']}  [{flag}]")
        print(f"      -> {r['targets'][:70]}")
    return agg

print("\nQuerying DGIdb (this takes a minute)...")
aA=summarise(dgidb(A),"SetA_lipid_platelet",A) if A else pd.DataFrame()
aB=summarise(dgidb(B),"SetB_cytokine_immune",B) if B else pd.DataFrame()

print("\n"+"!"*70)
print("COMPUTATIONAL HYPOTHESES ONLY - NOT CLINICAL RECOMMENDATIONS.")
print("Verify ART flags (Liverpool HIV Drug Interactions) and have a clinical")
print("pharmacologist review before any manuscript claim.")
print("!"*70)
print("\nSaved: pathway_drugs_SetA_*.csv / pathway_drugs_SetB_*.csv")
