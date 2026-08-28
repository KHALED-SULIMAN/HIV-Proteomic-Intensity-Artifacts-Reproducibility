from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

LINK = ROOT / "data/processed/master_patient_linkage.csv"

PROT = (
    ROOT
    / "external/Immunometabolism_HIV/processing"
    / "COCOMO_proteomics_olink_filt.csv"
)

MET = (
    ROOT
    / "external/Immunometabolism_HIV/processing"
    / "COCOMO_metabolomics_log2_norm_filt.csv"
)

CLIN = (
    ROOT
    / "external/Immunometabolism_HIV/processing"
    / "clinical_data_clean_with_clusters_and_categories.csv"
)

OUT = ROOT / "data/intermediate_public"
OUT.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# 1. Validated 89-participant linkage
# ------------------------------------------------------------
link = pd.read_csv(LINK)

link["COCOMO_ID"] = (
    pd.to_numeric(link["COCOMO_ID"], errors="raise")
    .astype("int64")
    .astype(str)
)

if len(link) != 89:
    raise RuntimeError(
        f"Expected 89 linked participants, found {len(link)}"
    )

if link["COCOMO_ID"].nunique() != 89:
    raise RuntimeError("COCOMO_ID values are not unique.")

ids = link["COCOMO_ID"].tolist()

print("Validated linkage:", len(ids), "participants")

# ------------------------------------------------------------
# 2. Proteomics
# ------------------------------------------------------------
prot = pd.read_csv(PROT, index_col=0)

prot.index = (
    pd.to_numeric(prot.index, errors="raise")
    .astype("int64")
    .astype(str)
)

missing_prot = [x for x in ids if x not in prot.index]

if missing_prot:
    raise RuntimeError(
        f"Proteomics missing {len(missing_prot)} participants: "
        f"{missing_prot}"
    )

prot89 = prot.loc[ids].copy()
prot89.index.name = "COCOMO_ID"

# ------------------------------------------------------------
# 3. Metabolomics
# ------------------------------------------------------------
met = pd.read_csv(MET, index_col=0)

met.index = (
    pd.to_numeric(met.index, errors="raise")
    .astype("int64")
    .astype(str)
)

missing_met = [x for x in ids if x not in met.index]

if missing_met:
    raise RuntimeError(
        f"Metabolomics missing {len(missing_met)} participants: "
        f"{missing_met}"
    )

met89 = met.loc[ids].copy()
met89.index.name = "COCOMO_ID"

# ------------------------------------------------------------
# 4. Clinical/outcome variables
# ------------------------------------------------------------
clinical = pd.read_csv(CLIN)

clinical["COCOMO_ID"] = (
    pd.to_numeric(clinical["COCOMO_ID"], errors="raise")
    .astype("int64")
    .astype(str)
)

if clinical["COCOMO_ID"].duplicated().any():
    raise RuntimeError("Clinical COCOMO_ID values are not unique.")

clinical = clinical.set_index("COCOMO_ID")

missing_clin = [x for x in ids if x not in clinical.index]

if missing_clin:
    raise RuntimeError(
        f"Clinical table missing {len(missing_clin)} participants: "
        f"{missing_clin}"
    )

wanted = [
    "METS",
    "AGE",
    "GENDER",
    "BMI",
    "Duration",
    "CD4",
    "CD8",
    "central_obesity",
    "hypertension",
    "diabetes",
    "VAT",
    "SAT",
    "Tgl",
    "Hdl",
    "Ldl",
    "cluster",
]

present = [x for x in wanted if x in clinical.columns]

outcomes = clinical.loc[ids, present].copy()
outcomes.index.name = "COCOMO_ID"

# ------------------------------------------------------------
# 5. Integrity checks
# ------------------------------------------------------------
assert prot89.index.tolist() == ids
assert met89.index.tolist() == ids
assert outcomes.index.tolist() == ids

# ------------------------------------------------------------
# 6. Save public intermediate matrices
# ------------------------------------------------------------
link.to_csv(
    OUT / "master_patient_linkage.csv",
    index=False
)

prot89.to_csv(
    OUT / "proteomics_olink.csv"
)

met89.to_csv(
    OUT / "metabolomics.csv"
)

outcomes.to_csv(
    OUT / "outcomes.csv"
)

# ------------------------------------------------------------
# 7. Summary
# ------------------------------------------------------------
print()
print("=" * 65)
print("SUCCESS")
print("=" * 65)
print("Participant linkage :", link.shape)
print("Proteomics          :", prot89.shape)
print("Metabolomics        :", met89.shape)
print("Outcomes            :", outcomes.shape)

print("\nFiles created:")

for name in [
    "master_patient_linkage.csv",
    "proteomics_olink.csv",
    "metabolomics.csv",
    "outcomes.csv",
]:
    path = OUT / name
    print(
        f"{name:30s} "
        f"{path.stat().st_size / (1024 * 1024):.2f} MB"
    )
