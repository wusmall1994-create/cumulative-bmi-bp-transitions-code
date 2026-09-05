from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadstat
import scipy.linalg
import statsmodels.api as sm


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = Path(os.getenv("BMI_BP_OUTPUT_DIR", ROOT / "outputs"))
DATA = Path(os.getenv("CHNS_DATA_DIR", ROOT / "data" / "chns"))
OUT = OUTPUT_ROOT / "chns_analysis"
OUT.mkdir(parents=True, exist_ok=True)


def five_state(sbp: pd.Series, dbp: pd.Series, med: pd.Series) -> pd.Series:
    state = pd.Series(pd.NA, index=sbp.index, dtype="Int64")
    valid = sbp.notna() & dbp.notna() & med.isin([0, 1])
    normal = (sbp < 120) & (dbp < 80)
    below_htn = (sbp < 140) & (dbp < 90)
    htn = (sbp >= 140) | (dbp >= 90)
    state.loc[valid & med.eq(0) & normal] = 1
    state.loc[valid & med.eq(0) & below_htn & ~normal] = 2
    state.loc[valid & med.eq(0) & htn] = 3
    state.loc[valid & med.eq(1) & below_htn] = 4
    state.loc[valid & med.eq(1) & htn] = 5
    return state


def build_design(frame: pd.DataFrame, continuous: list[str], categorical: list[str]) -> pd.DataFrame:
    parts = []
    for column in continuous:
        values = pd.to_numeric(frame[column], errors="coerce")
        missing = values.isna().astype(float)
        median = values.median()
        if pd.isna(median):
            median = 0.0
        values = values.fillna(median).astype(float)
        sd = values.std(ddof=0)
        if not np.isfinite(sd) or sd == 0:
            sd = 1.0
        parts.append(pd.DataFrame({column: (values - values.mean()) / sd, f"{column}_missing": missing}, index=frame.index))
    for column in categorical:
        values = frame[column].astype("string").fillna("MISSING")
        parts.append(pd.get_dummies(values, prefix=column, drop_first=True, dtype=float))
    design = pd.concat(parts, axis=1) if parts else pd.DataFrame(index=frame.index)
    design = design.loc[:, design.nunique(dropna=False) > 1]
    design = sm.add_constant(design.astype(float), has_constant="add")
    rank = np.linalg.matrix_rank(design.to_numpy())
    if rank < design.shape[1]:
        _, _, pivot = scipy.linalg.qr(design.to_numpy(), mode="economic", pivoting=True)
        design = design.iloc[:, sorted(pivot[:rank])]
    return design


def effective_sample_size(weight: pd.Series) -> float:
    weight = pd.to_numeric(weight, errors="coerce").dropna()
    return float(weight.sum() ** 2 / (weight**2).sum()) if len(weight) and (weight**2).sum() else np.nan


def fit_sw(name: str, frame: pd.DataFrame, outcome: pd.Series, continuous: list[str], categorical: list[str]):
    y = pd.to_numeric(outcome, errors="coerce")
    eligible = y.isin([0, 1])
    data = frame.loc[eligible].copy()
    yy = y.loc[eligible].astype(float)
    x = build_design(data, continuous, categorical)
    try:
        fit = sm.GLM(yy, x, family=sm.families.Binomial()).fit(maxiter=200, disp=0)
        model_type = "GLM"
    except Exception:
        fit = sm.GLM(yy, x, family=sm.families.Binomial()).fit_regularized(alpha=1e-6, L1_wt=0, maxiter=500)
        model_type = "ridge_regularized_GLM"
    p = pd.Series(np.nan, index=frame.index, dtype=float)
    p.loc[eligible] = np.asarray(fit.predict(x), dtype=float)
    p = p.clip(0.02, 0.98)
    numerator = float(yy.mean())
    sw = pd.Series(np.nan, index=frame.index, dtype=float)
    selected = eligible & y.eq(1)
    sw.loc[selected] = numerator / p.loc[selected]
    low, high = sw.dropna().quantile([0.01, 0.99]).tolist()
    sw_t = sw.clip(low, high)
    diag = {
        "model": name, "model_type": model_type, "eligible_n": int(eligible.sum()),
        "selected_n": int(selected.sum()), "selected_rate": numerator,
        "probability_min": float(p.loc[eligible].min()), "probability_p01": float(p.loc[eligible].quantile(0.01)),
        "probability_median": float(p.loc[eligible].median()), "probability_p99": float(p.loc[eligible].quantile(0.99)),
        "probability_max": float(p.loc[eligible].max()), "sw_p99": float(sw.dropna().quantile(0.99)),
        "sw_max": float(sw.dropna().max()), "truncate_low": low, "truncate_high": high,
        "ess_truncated": effective_sample_size(sw_t.dropna()), "design_columns": int(x.shape[1]),
    }
    return p, sw_t, diag


def as_binary(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").where(pd.to_numeric(series, errors="coerce").isin([0, 1]))


# Load only required public-use variables.
pexam_cols = [
    "IDind", "WAVE", "HEIGHT", "WEIGHT", "SYSTOL1", "DIASTOL1", "SYSTOL2", "DIASTOL2",
    "SYSTOL3", "DIASTOL3", "U22", "U24", "U24A", "U24J", "U24L", "U25", "U27", "U40", "U48A",
    "COMMID", "T1", "T2",
]
pexam, _ = pyreadstat.read_sas7bdat(DATA / "pexam_00.sas7bdat", usecols=pexam_cols)
surveys, _ = pyreadstat.read_sas7bdat(DATA / "surveys_pub_12.sas7bdat", usecols=["Idind", "wave", "age", "urban", "stratum", "commid"])
education, _ = pyreadstat.read_sas7bdat(DATA / "educ_12.sas7bdat", usecols=["IDind", "WAVE", "A11"])
master, _ = pyreadstat.read_sas7bdat(DATA / "mast_pub_12.sas7bdat", usecols=["Idind", "GENDER", "DOD_RPT"])

for frame in [pexam, surveys, education, master]:
    frame.columns = frame.columns.str.lower()

pexam = pexam.loc[pexam["wave"].isin([2000, 2004, 2006, 2009, 2011])].copy()
surveys = surveys.loc[surveys["wave"].isin([2000, 2004, 2006, 2009, 2011])].copy()
education = education.loc[education["wave"].eq(2006), ["idind", "a11"]].drop_duplicates("idind")
master = master.drop_duplicates("idind")

if pexam.duplicated(["idind", "wave"]).any():
    raise RuntimeError("Duplicate person-wave rows found in CHNS physical examination file")
if surveys.duplicated(["idind", "wave"]).any():
    raise RuntimeError("Duplicate person-wave rows found in CHNS survey participation file")

long = pexam.merge(surveys, on=["idind", "wave"], how="left", suffixes=("", "_survey"))
long = long.merge(master[["idind", "gender", "dod_rpt"]], on="idind", how="left")
long = long.merge(education, on="idind", how="left")
long["person_id"] = pd.to_numeric(long["idind"], errors="coerce").astype("Int64").astype("string")

# Objective BMI and blood pressure cleaning.
long["height"] = pd.to_numeric(long["height"], errors="coerce").where(pd.to_numeric(long["height"], errors="coerce").between(100, 220))
long["weight"] = pd.to_numeric(long["weight"], errors="coerce").where(pd.to_numeric(long["weight"], errors="coerce").between(20, 250))
long["bmi"] = long["weight"] / (long["height"] / 100) ** 2
long["bmi"] = long["bmi"].where(long["bmi"].between(10, 70))

for kind, cols in {
    "sbp": ["systol1", "systol2", "systol3"],
    "dbp": ["diastol1", "diastol2", "diastol3"],
}.items():
    lower, upper = (50, 300) if kind == "sbp" else (30, 200)
    values = long[cols].apply(pd.to_numeric, errors="coerce").where(lambda x: x.ge(lower) & x.le(upper))
    long[f"{kind}_n"] = values.notna().sum(axis=1)
    long[kind] = values.mean(axis=1).where(long[f"{kind}_n"].ge(2))
long.loc[~(long["sbp_n"].ge(2) & long["dbp_n"].ge(2)), ["sbp", "dbp"]] = np.nan

long["med_current"] = np.select(
    [long["u24"].eq(1), long["u24"].eq(0) | long["u22"].eq(0)], [1.0, 0.0], default=np.nan
)
long["state5"] = five_state(long["sbp"], long["dbp"], long["med_current"])
long["current_smoker"] = np.select(
    [long["u27"].eq(1), long["u27"].eq(0) | long["u25"].eq(0)], [1.0, 0.0], default=np.nan
)
long["current_drinker"] = as_binary(long["u40"])
long["diabetes"] = as_binary(long["u24a"])
long["mi_history"] = as_binary(long["u24j"])
long["stroke_history"] = as_binary(long["u24l"])
long["cvd_history"] = long[["mi_history", "stroke_history"]].max(axis=1)

# Wide person-level file.
value_cols = [
    "bmi", "state5", "sbp", "dbp", "med_current", "age", "urban", "stratum", "commid", "t1", "t2",
    "current_smoker", "current_drinker", "diabetes", "cvd_history", "u48a", "a11",
]
wide = long.pivot(index="person_id", columns="wave", values=value_cols)
wide.columns = [f"{name}_{int(year)}" for name, year in wide.columns]
wide = wide.reset_index()

for year in [2000, 2004, 2006]:
    if f"bmi_{year}" not in wide:
        wide[f"bmi_{year}"] = np.nan

bmi_values = wide[["bmi_2000", "bmi_2004", "bmi_2006"]].apply(pd.to_numeric, errors="coerce")
wide[["bmi_2000", "bmi_2004", "bmi_2006"]] = bmi_values
wide["exposure_complete"] = bmi_values.notna().all(axis=1).astype(int)
wide["bmi_mean"] = bmi_values.mean(axis=1)
wide["bmi_sd"] = bmi_values.std(axis=1, ddof=1)
wide["bmi_cv"] = 100 * wide["bmi_sd"] / wide["bmi_mean"]
wide["bmi_arv"] = (abs(wide["bmi_2004"] - wide["bmi_2000"]) + abs(wide["bmi_2006"] - wide["bmi_2004"])) / 2
wide["bmi_auc"] = (wide["bmi_2000"] + wide["bmi_2004"]) / 2 * 4 + (wide["bmi_2004"] + wide["bmi_2006"]) / 2 * 2
wide["bmi_twmean"] = wide["bmi_auc"] / 6
excess = (bmi_values - 25).clip(lower=0)
wide["excess_bmi25"] = (excess["bmi_2000"] + excess["bmi_2004"]) / 2 * 4 + (excess["bmi_2004"] + excess["bmi_2006"]) / 2 * 2
wide["bmi_change_pct"] = (wide["bmi_2006"] - wide["bmi_2000"]) / wide["bmi_2000"] * 100
wide["bmi_change_cat_code"] = np.select([wide["bmi_change_pct"].lt(-5), wide["bmi_change_pct"].ge(5)], [-1, 1], default=0)

# VIM: SD / mean^beta, rescaled to the cohort mean level.
vim_source = wide.loc[wide["exposure_complete"].eq(1) & wide["bmi_sd"].gt(0), ["bmi_mean", "bmi_sd"]]
vim_beta = float(np.polyfit(np.log(vim_source["bmi_mean"]), np.log(vim_source["bmi_sd"]), 1)[0])
vim_scale = float(wide.loc[wide["exposure_complete"].eq(1), "bmi_mean"].mean() ** vim_beta)
wide["bmi_vim"] = wide["bmi_sd"] / (wide["bmi_mean"] ** vim_beta) * vim_scale

# Baseline covariates.
wide["female"] = wide["person_id"].map(master.assign(person_id=pd.to_numeric(master["idind"], errors="coerce").astype("Int64").astype("string")).set_index("person_id")["gender"]).eq(2).astype(int)
wide["age2006"] = wide["age_2006"]
wide["educ_years"] = wide["a11_2006"]
wide["province"] = wide["t1_2006"]
wide["urban"] = wide["urban_2006"]
wide["commid"] = wide["commid_2006"]
wide["self_health"] = wide["u48a_2006"]
for v in ["current_smoker", "current_drinker", "diabetes", "cvd_history"]:
    wide[v] = wide[f"{v}_2006"]
wide["age10"] = wide["age2006"] / 10
wide["age2006_sq"] = wide["age2006"] ** 2

# Exposure target requires actual participation in 2000 to preserve positivity.
participated_2000 = set(surveys.loc[surveys["wave"].eq(2000), "idind"].dropna().astype(int).astype(str))
wide["participated_2000"] = wide["person_id"].isin(participated_2000).astype(int)
target_mask = wide["age2006"].ge(50) & wide["state5_2006"].notna() & wide["bmi_2006"].notna() & wide["participated_2000"].eq(1)
target = wide.loc[target_mask].copy()

continuous = ["age2006", "age2006_sq", "educ_years", "self_health", "bmi_2006"]
categorical = ["female", "current_smoker", "current_drinker", "diabetes", "cvd_history", "province", "urban", "state5_2006"]
p_exp, sw_exp, diag_exp = fit_sw("complete_bmi_history", target, target["exposure_complete"], continuous, categorical)
target["sw_exposure"] = sw_exp
sw_map = target.set_index("person_id")["sw_exposure"]

cohort_mask = wide["age2006"].ge(50) & wide["state5_2006"].notna() & wide["exposure_complete"].eq(1)
cohort = wide.loc[cohort_mask].copy()
cohort["sw_exposure"] = cohort["person_id"].map(sw_map)
if cohort["sw_exposure"].isna().any():
    raise RuntimeError("CHNS analysis cohort contains participants outside the exposure-weight target")

cohort["state_observed_2009"] = cohort["state5_2009"].notna().astype(int)
cohort["state_observed_2011"] = cohort["state5_2011"].notna().astype(int)
follow_cont = ["age2006", "educ_years", "self_health", "bmi_mean", "bmi_vim", "excess_bmi25"]
follow_cat = ["female", "current_smoker", "current_drinker", "diabetes", "cvd_history", "province", "urban", "state5_2006"]
p09, sw09, diag09 = fit_sw("state_observed_2009", cohort, cohort["state_observed_2009"], follow_cont, follow_cat)
cohort["sw_state2009"] = sw09
cond11 = cohort["state_observed_2009"].eq(1)
cohort["state5_2009_model"] = cohort["state5_2009"]
p11, sw11, diag11 = fit_sw(
    "state_observed_2011_conditional_on_2009", cohort.loc[cond11], cohort.loc[cond11, "state_observed_2011"],
    follow_cont, follow_cat + ["state5_2009_model"],
)
cohort["sw_state2011_cond"] = np.nan
cohort.loc[cond11, "sw_state2011_cond"] = sw11
cohort["analysis_weight_2009"] = cohort["sw_exposure"] * cohort["sw_state2009"]
cohort["analysis_weight_2011"] = cohort["sw_exposure"] * cohort["sw_state2009"] * cohort["sw_state2011_cond"]

# Standardize exposures in the CHNS replication cohort.
for var in ["excess_bmi25", "bmi_auc", "bmi_mean", "bmi_vim", "bmi_arv", "bmi_cv", "bmi_sd", "bmi_change_pct"]:
    cohort[f"{var}_z"] = (cohort[var] - cohort[var].mean()) / cohort[var].std(ddof=0)

# Build interval data with harmonized endpoints.
interval_rows = []
for start, end in [(2006, 2009), (2009, 2011)]:
    for _, base in cohort.iterrows():
        o, dest = base[f"state5_{start}"], base[f"state5_{end}"]
        row = {
            "person_id": base["person_id"], "interval": f"{start}-{end}", "start_year": start, "end_year": end,
            "origin": o, "destination": dest, "age_start10": (base["age2006"] + start - 2006) / 10,
            "model_weight": base[f"analysis_weight_{end}"], "base_weight": base["sw_exposure"],
            "commid": base["commid"], "female": base["female"], "educ_years": base["educ_years"],
            "current_smoker": base["current_smoker"], "current_drinker": base["current_drinker"],
            "diabetes": base["diabetes"], "cvd_history": base["cvd_history"], "self_health": base["self_health"],
            "province": base["province"], "urban": base["urban"], "bmi_change_cat_code": base["bmi_change_cat_code"],
        }
        for var in ["excess_bmi25", "bmi_auc", "bmi_mean", "bmi_vim", "bmi_arv", "bmi_cv", "bmi_sd", "bmi_change_pct"]:
            row[var] = base[var]
            row[f"{var}_z"] = base[f"{var}_z"]
        row["hypertension_onset"] = int(dest in [3, 4, 5]) if pd.notna(o) and pd.notna(dest) and o in [1, 2] else np.nan
        row["treatment_initiation"] = int(dest in [4, 5]) if pd.notna(o) and pd.notna(dest) and o == 3 else np.nan
        row["untreated_bp_improvement"] = int(dest in [1, 2]) if pd.notna(o) and pd.notna(dest) and o == 3 and dest in [1, 2, 3] else np.nan
        row["control_loss"] = int(dest == 5) if pd.notna(o) and pd.notna(dest) and o == 4 and dest in [4, 5] else np.nan
        row["control_achievement"] = int(dest == 4) if pd.notna(o) and pd.notna(dest) and o == 5 and dest in [4, 5] else np.nan
        interval_rows.append(row)
intervals = pd.DataFrame(interval_rows)

endpoint_names = ["hypertension_onset", "treatment_initiation", "untreated_bp_improvement", "control_loss", "control_achievement"]
event_gate = []
for endpoint in endpoint_names:
    eligible = intervals[endpoint].notna()
    event_gate.append({
        "endpoint": endpoint, "eligible_intervals": int(eligible.sum()), "events": int(intervals.loc[eligible, endpoint].sum()),
        "nonevents": int(eligible.sum() - intervals.loc[eligible, endpoint].sum()),
        "participants": int(intervals.loc[eligible, "person_id"].nunique()),
        "weighted_eligible": int((eligible & intervals["model_weight"].gt(0)).sum()),
    })
event_gate = pd.DataFrame(event_gate)

transition_rows = []
for interval, subset in intervals.loc[intervals["origin"].notna() & intervals["destination"].notna()].groupby("interval"):
    matrix = pd.crosstab(subset["origin"].astype(int), subset["destination"].astype(int))
    for o in range(1, 6):
        for dest in range(1, 6):
            transition_rows.append({"interval": interval, "origin": o, "destination": dest,
                                    "n": int(matrix.loc[o, dest]) if o in matrix.index and dest in matrix.columns else 0})
transition_counts = pd.DataFrame(transition_rows)

attrition = pd.DataFrame([
    {"year": 2006, "baseline_cohort": len(cohort), "state_observed": int(cohort["state5_2006"].notna().sum())},
    {"year": 2009, "baseline_cohort": len(cohort), "state_observed": int(cohort["state5_2009"].notna().sum())},
    {"year": 2011, "baseline_cohort": len(cohort), "state_observed": int(cohort["state5_2011"].notna().sum())},
])

cohort.to_csv(OUT / "chns_baseline_2006.csv", index=False, encoding="utf-8")
intervals.to_csv(OUT / "chns_intervals.csv", index=False, encoding="utf-8")
pd.DataFrame([diag_exp, diag09, diag11]).to_csv(OUT / "chns_ipw_diagnostics.csv", index=False, encoding="utf-8-sig")
event_gate.to_csv(OUT / "chns_endpoint_event_gate.csv", index=False, encoding="utf-8-sig")
transition_counts.to_csv(OUT / "chns_transition_counts.csv", index=False, encoding="utf-8-sig")
attrition.to_csv(OUT / "chns_attrition.csv", index=False, encoding="utf-8-sig")

summary = {
    "exposure_target_n": int(len(target)), "baseline_n": int(len(cohort)),
    "state_observed_2009": int(cohort["state5_2009"].notna().sum()),
    "state_observed_2011": int(cohort["state5_2011"].notna().sum()),
    "both_followups": int(cohort[["state5_2009", "state5_2011"]].notna().all(axis=1).sum()),
    "vim_beta": vim_beta, "event_gate": event_gate.to_dict("records"),
    "weight_note": "CHNS does not provide nationally representative sampling weights; community-clustered inference is required.",
}
(OUT / "chns_replication_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps(summary, indent=2, ensure_ascii=False))
