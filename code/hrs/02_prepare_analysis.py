from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = Path(os.getenv("BMI_BP_OUTPUT_DIR", ROOT / "outputs"))
SOURCE = OUTPUT_ROOT / "hrs_formal_cohort"
OUT = OUTPUT_ROOT / "hrs_analysis"
OUT.mkdir(parents=True, exist_ok=True)

ID_DTYPE = {"person_id": "string", "hhid": "string", "pn": "string"}
baseline = pd.read_csv(SOURCE / "hrs_o1_baseline_2014.csv", dtype=ID_DTYPE, low_memory=False)
long = pd.read_csv(SOURCE / "hrs_o1_person_wave_2014_2022.csv", dtype=ID_DTYPE, low_memory=False)

state_labels = {
    1: "Normal untreated",
    2: "Elevated untreated",
    3: "Untreated hypertension",
    4: "Treated controlled",
    5: "Treated uncontrolled",
    6: "Death",
}


def accaha_state(sbp: pd.Series, dbp: pd.Series, med: pd.Series) -> pd.Series:
    state = pd.Series(pd.NA, index=sbp.index, dtype="Int64")
    valid = sbp.notna() & dbp.notna() & med.isin([0, 1])
    normal = (sbp < 120) & (dbp < 80)
    elevated = (sbp >= 120) & (sbp < 130) & (dbp < 80)
    hypertension = (sbp >= 130) | (dbp >= 80)
    state.loc[valid & med.eq(0) & normal] = 1
    state.loc[valid & med.eq(0) & elevated] = 2
    state.loc[valid & med.eq(0) & hypertension] = 3
    state.loc[valid & med.eq(1) & ~hypertension] = 4
    state.loc[valid & med.eq(1) & hypertension] = 5
    return state


# Baseline covariates and exposure parameterization.
baseline["female"] = baseline["ragender"].eq(2).astype(int)
baseline["married_partnered"] = baseline["r12mstat"].isin([1, 2, 3]).astype(float)
baseline["current_smoker"] = baseline["r12smoken"]
baseline["current_drinker"] = baseline["r12drink"]
baseline["cvd_history"] = baseline[["r12hearte", "r12stroke"]].max(axis=1)
baseline["age10"] = baseline["r12agey_e"] / 10.0
baseline["bmi_change_pct"] = (baseline["r12pmbmi"] - baseline["r8pmbmi"]) / baseline["r8pmbmi"] * 100
baseline["bmi_change_cat"] = pd.cut(
    baseline["bmi_change_pct"],
    bins=[-np.inf, -5, 5, np.inf],
    labels=["Loss >5%", "Stable +/-5%", "Gain >5%"],
    right=False,
).astype("string")
baseline["bmi_change_cat_code"] = np.select(
    [baseline["bmi_change_pct"].lt(-5), baseline["bmi_change_pct"].ge(5)],
    [-1, 1],
    default=0,
)
baseline["obese_2014"] = baseline["r12pmbmi"].ge(30).astype(int)
baseline["age65"] = baseline["r12agey_e"].ge(65).astype(int)
baseline["rural"] = baseline["beale2023_14"].eq(3).astype(int)

exposure_vars = ["excess_bmi25", "bmi_auc", "bmi_twmean", "bmi_mean", "bmi_vim", "bmi_arv", "bmi_cv", "bmi_sd", "bmi_change_pct"]
scaling = []
for var in exposure_vars:
    mean = float(baseline[var].mean())
    sd = float(baseline[var].std(ddof=0))
    baseline[f"{var}_z"] = (baseline[var] - mean) / sd
    scaling.append({"variable": var, "mean": mean, "sd": sd})

baseline["burden_q"] = pd.qcut(baseline["excess_bmi25"], 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop").astype("string")
baseline["vim_q"] = pd.qcut(baseline["bmi_vim"], 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop").astype("string")

# ACC/AHA sensitivity states reconstructed from the same measured BP and medication fields.
long["state5_accaha"] = accaha_state(long["sbp"], long["dbp"], long["med_current"])

covariates = [
    "r12agey_e", "female", "ragender", "raracem", "raedyrs", "r12mstat", "married_partnered",
    "r12shlt", "r12smoken", "current_smoker", "r12drink", "current_drinker", "r12diabe",
    "r12hearte", "r12stroke", "cvd_history", "r12cesd", "region14", "beale2023_14",
    "obese_2014", "age65", "rural", "secu", "stratum", "base_weight_candidate",
]
exposure_columns = exposure_vars + [f"{v}_z" for v in exposure_vars] + ["burden_q", "vim_q", "bmi_change_cat", "bmi_change_cat_code"]
base_map = baseline.set_index("person_id")
long_index = long.set_index(["person_id", "year"])

rows = []
for start, end in [(2014, 2018), (2018, 2022)]:
    for person_id, base in base_map.iterrows():
        start_row = long_index.loc[(person_id, start)]
        end_row = long_index.loc[(person_id, end)]
        row = {
            "person_id": person_id,
            "interval": f"{start}-{end}",
            "start_year": start,
            "end_year": end,
            "age_start": float(base["r12agey_e"] + (start - 2014)),
            "age_start10": float((base["r12agey_e"] + (start - 2014)) / 10.0),
            "origin": start_row["state5"],
            "destination": end_row["state5"],
            "destination6": end_row["state6_death"],
            "origin_accaha": start_row["state5_accaha"],
            "destination_accaha": end_row["state5_accaha"],
            "died": int(end_row["died_by_wave"] == 1),
            "known_six_state": int(pd.notna(end_row["state5"]) or end_row["died_by_wave"] == 1),
            "model_weight": base[f"analysis_weight_{end}"],
            "base_weight": base["base_weight_candidate"],
            "secu": base["secu"],
            "stratum": base["stratum"],
        }
        for col in covariates + exposure_columns:
            if col not in row:
                row[col] = base[col]

        o, d = row["origin"], row["destination"]
        oa, da = row["origin_accaha"], row["destination_accaha"]
        row["hypertension_onset"] = int(d in [3, 4, 5]) if o in [1, 2] and pd.notna(d) else np.nan
        row["treatment_initiation"] = int(d in [4, 5]) if o == 3 and pd.notna(d) else np.nan
        row["untreated_bp_improvement"] = int(d in [1, 2]) if o == 3 and d in [1, 2, 3] else np.nan
        row["control_loss"] = int(d == 5) if o == 4 and d in [4, 5] else np.nan
        row["control_achievement"] = int(d == 4) if o == 5 and d in [4, 5] else np.nan
        row["hypertension_onset_accaha"] = int(da in [3, 4, 5]) if pd.notna(oa) and pd.notna(da) and oa in [1, 2] else np.nan
        row["control_loss_accaha"] = int(da == 5) if pd.notna(oa) and pd.notna(da) and oa == 4 and da in [4, 5] else np.nan
        row["control_achievement_accaha"] = int(da == 4) if pd.notna(oa) and pd.notna(da) and oa == 5 and da in [4, 5] else np.nan
        row["death_next_wave"] = int(row["died"]) if pd.notna(o) and row["known_six_state"] == 1 else np.nan
        rows.append(row)

intervals = pd.DataFrame(rows)

# Transition matrices and event gates.
transition_rows = []
for start, end in [(2014, 2018), (2018, 2022)]:
    subset = intervals.loc[intervals["interval"].eq(f"{start}-{end}") & intervals["origin"].notna() & intervals["destination"].notna()]
    matrix = pd.crosstab(subset["origin"].astype(int), subset["destination"].astype(int))
    for origin in range(1, 6):
        for destination in range(1, 6):
            n = int(matrix.loc[origin, destination]) if origin in matrix.index and destination in matrix.columns else 0
            transition_rows.append(
                {
                    "interval": f"{start}-{end}", "origin": origin, "destination": destination, "n": n,
                    "origin_label": state_labels[origin], "destination_label": state_labels[destination],
                }
            )
transition_counts = pd.DataFrame(transition_rows)

endpoint_names = [
    "hypertension_onset", "treatment_initiation", "untreated_bp_improvement",
    "control_loss", "control_achievement", "death_next_wave",
    "hypertension_onset_accaha", "control_loss_accaha", "control_achievement_accaha",
]
event_gate = []
for endpoint in endpoint_names:
    eligible = intervals[endpoint].notna()
    event_gate.append(
        {
            "endpoint": endpoint,
            "eligible_intervals": int(eligible.sum()),
            "events": int(intervals.loc[eligible, endpoint].sum()),
            "nonevents": int(eligible.sum() - intervals.loc[eligible, endpoint].sum()),
            "participants": int(intervals.loc[eligible, "person_id"].nunique()),
            "weighted_eligible": int((eligible & intervals["model_weight"].gt(0)).sum()) if endpoint != "death_next_wave" else int((eligible & intervals["base_weight"].gt(0)).sum()),
        }
    )
event_gate = pd.DataFrame(event_gate)

# Covariate audit and exposure diagnostics.
audit_vars = [
    "r12agey_e", "ragender", "raracem", "raedyrs", "r12mstat", "r12shlt", "r12smoken",
    "r12drink", "r12diabe", "r12hearte", "r12stroke", "r12cesd", "region14", "beale2023_14",
] + exposure_vars
covariate_audit = pd.DataFrame(
    {
        "variable": audit_vars,
        "nonmissing_n": [int(baseline[v].notna().sum()) for v in audit_vars],
        "missing_n": [int(baseline[v].isna().sum()) for v in audit_vars],
        "missing_percent": [float(100 * baseline[v].isna().mean()) for v in audit_vars],
        "unique_values": [int(baseline[v].nunique(dropna=True)) for v in audit_vars],
    }
)
exposure_summary = baseline[exposure_vars].describe(percentiles=[0.25, 0.5, 0.75]).T.reset_index(names="variable")
exposure_correlations = baseline[exposure_vars].corr().reset_index(names="variable")

# A compact descriptive table by cumulative excess BMI quartile.
def summarize_group(frame: pd.DataFrame, label: str) -> dict:
    return {
        "group": label,
        "n": len(frame),
        "age_mean": frame["r12agey_e"].mean(),
        "female_percent": 100 * frame["female"].mean(),
        "black_percent": 100 * frame["raracem"].eq(2).mean(),
        "education_mean": frame["raedyrs"].mean(),
        "smoker_percent": 100 * frame["current_smoker"].mean(),
        "diabetes_percent": 100 * frame["r12diabe"].mean(),
        "cvd_percent": 100 * frame["cvd_history"].mean(),
        "bmi_mean": frame["bmi_mean"].mean(),
        "excess_bmi25_mean": frame["excess_bmi25"].mean(),
        "vim_mean": frame["bmi_vim"].mean(),
        "baseline_htn_percent": 100 * frame["state5_2014"].isin([3, 4, 5]).mean(),
    }


table1 = [summarize_group(baseline, "Overall")]
for quartile in ["Q1", "Q2", "Q3", "Q4"]:
    table1.append(summarize_group(baseline.loc[baseline["burden_q"].eq(quartile)], quartile))
table1 = pd.DataFrame(table1)

# Write analysis-ready files.
baseline.to_csv(OUT / "analysis_baseline.csv", index=False, encoding="utf-8")
intervals.to_csv(OUT / "analysis_intervals.csv", index=False, encoding="utf-8")
pd.DataFrame(scaling).to_csv(OUT / "exposure_scaling.csv", index=False, encoding="utf-8-sig")
covariate_audit.to_csv(OUT / "covariate_missingness.csv", index=False, encoding="utf-8-sig")
exposure_summary.to_csv(OUT / "exposure_summary.csv", index=False, encoding="utf-8-sig")
exposure_correlations.to_csv(OUT / "exposure_correlations.csv", index=False, encoding="utf-8-sig")
transition_counts.to_csv(OUT / "transition_counts.csv", index=False, encoding="utf-8-sig")
event_gate.to_csv(OUT / "endpoint_event_gate.csv", index=False, encoding="utf-8-sig")
table1.to_csv(OUT / "table1_by_burden_quartile.csv", index=False, encoding="utf-8-sig")

summary = {
    "baseline_n": int(len(baseline)),
    "interval_rows": int(len(intervals)),
    "maximum_covariate_missing_percent": float(covariate_audit["missing_percent"].max()),
    "corr_bmi_auc_twmean": float(baseline[["bmi_auc", "bmi_twmean"]].corr().iloc[0, 1]),
    "corr_excess_burden_mean_bmi": float(baseline[["excess_bmi25", "bmi_mean"]].corr().iloc[0, 1]),
    "weight_change_categories": baseline["bmi_change_cat"].value_counts(dropna=False).to_dict(),
    "event_gate": event_gate.to_dict("records"),
}
(OUT / "analysis_preparation_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

print(json.dumps(summary, indent=2, ensure_ascii=False))
print("\nTRANSITION COUNTS")
print(transition_counts.pivot_table(index=["interval", "origin"], columns="destination", values="n", fill_value=0).to_string())
