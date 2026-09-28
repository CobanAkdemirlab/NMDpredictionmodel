#1. Motif QC
# ==============================================================================
# MOTIF FEATURE QC
# ==============================================================================

MOTIF_PREFIXES = (
    "utr3_all.",
    "utr3_200.",
    "ptc_to_ejc.",
    "newutr_all.",
    "newutr_200.",
    "ptc_pm100.",
    "ejc_pm100.",
)

motif_cols = [
    col for col in X.columns
    if col.startswith(MOTIF_PREFIXES)
]

print("\n" + "=" * 80)
print("MOTIF FEATURE QC")
print("=" * 80)

print(f"Motif features found: {len(motif_cols)}")


# ------------------------------------------------------------------------------
# Calculate prevalence
# ------------------------------------------------------------------------------

motif_qc = []

for col in motif_cols:

    x = pd.to_numeric(
        X[col],
        errors="coerce"
    )

    n_available = x.notna().sum()
    n_missing = x.isna().sum()

    n_present = (x == 1).sum()
    n_absent = (x == 0).sum()

    prevalence = (
        n_present / n_available
        if n_available > 0
        else np.nan
    )

    motif_qc.append({
        "feature": col,
        "n_available": n_available,
        "n_missing": n_missing,
        "n_present": n_present,
        "n_absent": n_absent,
        "prevalence": prevalence
    })


motif_qc = pd.DataFrame(motif_qc)


# ------------------------------------------------------------------------------
# Flag uninformative motifs
# ------------------------------------------------------------------------------

motif_qc["all_zero"] = (
    motif_qc["n_present"] == 0
)

motif_qc["all_one"] = (
    motif_qc["n_absent"] == 0
)

motif_qc["constant"] = (
    motif_qc["all_zero"] |
    motif_qc["all_one"]
)


print(
    "All-zero motifs:",
    motif_qc["all_zero"].sum()
)

print(
    "All-one motifs:",
    motif_qc["all_one"].sum()
)

print(
    "Constant motifs:",
    motif_qc["constant"].sum()
)


# ------------------------------------------------------------------------------
# Save QC report
# ------------------------------------------------------------------------------

motif_qc.to_csv(
    "results/feature_qc/motif_feature_qc.csv",
    index=False
)

motif_constant_drop = (
    motif_qc.loc[
        motif_qc["constant"],
        "feature"
    ]
    .tolist()
)

X = X.drop(
    columns=motif_constant_drop
)

print(
    f"Dropped {len(motif_constant_drop)} "
    "constant motif features."
)
motif_qc["present_lt_5"] = (
    motif_qc["n_present"] < 5
)

motif_qc["present_lt_1pct"] = (
    motif_qc["prevalence"] < 0.01
)

