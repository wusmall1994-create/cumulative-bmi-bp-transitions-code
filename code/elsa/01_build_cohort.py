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
DATA = Path(os.getenv("ELSA_DATA_DIR", ROOT / "data" / "elsa"))
OUT = OUTPUT_ROOT / "elsa_analysis"
OUT.mkdir(parents=True, exist_ok=True)


def read_dta(name: str, columns: list[str]) -> pd.DataFrame:
    frame, _ = pyreadstat.read_dta(DATA / name, usecols=columns, apply_value_formats=False)
    frame.columns = frame.columns.str.lower()
    return frame


def plausible(series: pd.Series, low: float, high: float) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    return values.where(values.between(low, high))


def allowed(series: pd.Series, values: list[int]) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    return numeric.where(numeric.isin(values))


def indicator_from_value(series: pd.Series, positive_values: list[int]) -> pd.Series:
    result = pd.Series(np.nan, index=series.index, dtype=float)
    result.loc[series.notna()] = series.loc[series.notna()].isin(positive_values).astype(float)
    return result


def treatment_status(frame: pd.DataFrame) -> pd.Series:
    hemda = pd.to_numeric(frame["hemda"], errors="coerce")
    hemdab = pd.to_numeric(frame["hemdab"], errors="coerce")
    treated = (hemda.eq(1) | hemdab.eq(1)).astype(float)
    ambiguous = (~hemda.isin([-1, 1, 2])) | (~hemdab.isin([-1, 1, 2]))
    treated.loc[ambiguous & ~treated.eq(1)] = np.nan
    return treated


def five_state(sbp: pd.Series, dbp: pd.Series, treated: pd.Series) -> pd.Series:
    state = pd.Series(pd.NA, index=sbp.index, dtype="Int64")
    valid = sbp.notna() & dbp.notna() & treated.isin([0, 1])
    normal = (sbp < 120) & (dbp < 80)
    below_htn = (sbp < 140) & (dbp < 90)
    htn = (sbp >= 140) | (dbp >= 90)
    state.loc[valid & treated.eq(0) & normal] = 1
    state.loc[valid & treated.eq(0) & below_htn & ~normal] = 2
    state.loc[valid & treated.eq(0) & htn] = 3
    state.loc[valid & treated.eq(1) & below_htn] = 4
    state.loc[valid & treated.eq(1) & htn] = 5
    return state


def accaha_state(sbp: pd.Series, dbp: pd.Series, treated: pd.Series) -> pd.Series:
    state = pd.Series(pd.NA, index=sbp.index, dtype="Int64")
    valid = sbp.notna() & dbp.notna() & treated.isin([0, 1])
    normal = (sbp < 120) & (dbp < 80)
    elevated = (sbp >= 120) & (sbp < 130) & (dbp < 80)
    htn = (sbp >= 130) | (dbp >= 80)
    state.loc[valid & treated.eq(0) & normal] = 1
    state.loc[valid & treated.eq(0) & elevated] = 2
    state.loc[valid & treated.eq(0) & htn] = 3
    state.loc[valid & treated.eq(1) & ~htn] = 4
    state.loc[valid & treated.eq(1) & htn] = 5
    return state


def build_design(frame: pd.DataFrame, continuous: list[str], categorical: list[str]) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
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
    values = pd.to_numeric(weight, errors="coerce").dropna()
    denominator = float((values**2).sum())
    return float(values.sum() ** 2 / denominator) if len(values) and denominator else np.nan


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
    probability = pd.Series(np.nan, index=frame.index, dtype=float)
    probability.loc[eligible] = np.asarray(fit.predict(x), dtype=float)
    probability = probability.clip(0.02, 0.98)
    numerator = float(yy.mean())
    sw = pd.Series(np.nan, index=frame.index, dtype=float)
    selected = eligible & y.eq(1)
    sw.loc[selected] = numerator / probability.loc[selected]
    low, high = sw.dropna().quantile([0.01, 0.99]).tolist()
    sw_truncated = sw.clip(low, high)
    diagnostics = {
        "model": name,
        "model_type": model_type,
        "eligible_n": int(eligible.sum()),
        "selected_n": int(selected.sum()),
        "selected_rate": numerator,
        "probability_min": float(probability.loc[eligible].min()),
        "probability_p01": float(probability.loc[eligible].quantile(0.01)),
        "probability_median": float(probability.loc[eligible].median()),
        "probability_p99": float(probability.loc[eligible].quantile(0.99)),
        "probability_max": float(probability.loc[eligible].max()),
        "sw_p99": float(sw.dropna().quantile(0.99)),
        "sw_max": float(sw.dropna().max()),
        "truncate_low": float(low),
        "truncate_high": float(high),
        "ess_truncated": effective_sample_size(sw_truncated.dropna()),
        "design_columns": int(x.shape[1]),
    }
    return probability, sw_truncated, diagnostics


# Repeated measured BMI and baseline/outcome health visits.
w2 = read_dta("wave_2_nurse_data_v2.dta", ["idauniq", "bmival"])
w2["w2_nurse_present"] = 1
w2["bmi_w2"] = plausible(w2.pop("bmival"), 10, 70)

w4 = read_dta("wave_4_nurse_data.dta", ["idauniq", "bmival"])
w4["bmi_w4"] = plausible(w4.pop("bmival"), 10, 70)

w6n = read_dta(
    "wave_6_elsa_nurse_data_v2.dta",
    ["idauniq", "BMIVAL", "SYSVAL", "DIAVAL", "w6nurwt"],
)
w6n["bmi_w6"] = plausible(w6n.pop("bmival"), 10, 70)
w6n["sbp_w6"] = plausible(w6n.pop("sysval"), 60, 260)
w6n["dbp_w6"] = plausible(w6n.pop("diaval"), 30, 150)
w6n["w6nurwt"] = pd.to_numeric(w6n["w6nurwt"], errors="coerce").where(lambda x: x > 0)

w6c = read_dta(
    "wave_6_elsa_data_eul.dta",
    ["idauniq", "hemda", "hemdab", "gor", "idahhw6"],
)
w6c["treated_w6"] = treatment_status(w6c)
w6c["region"] = w6c["gor"].astype("string").str.strip()
w6c.loc[~w6c["region"].str.match(r"^E1200000[1-9]$", na=False), "region"] = pd.NA
w6c = w6c.drop(columns=["hemda", "hemdab", "gor"])

w89n = read_dta(
    "elsa_nurse_w8w9_data_eul.dta",
    ["idauniq", "wave", "sysval", "diaval", "w89nurwt_N", "w89nurwt_G", "idahhw8w9"],
)
w89n["sbp_out"] = plausible(w89n.pop("sysval"), 60, 260)
w89n["dbp_out"] = plausible(w89n.pop("diaval"), 30, 150)
w89n["outcome_wave"] = pd.to_numeric(w89n.pop("wave"), errors="coerce")
w89n["w89nurwt_n"] = pd.to_numeric(w89n["w89nurwt_n"], errors="coerce").where(lambda x: x > 0)
w89n["w89nurwt_g"] = pd.to_numeric(w89n["w89nurwt_g"], errors="coerce").where(lambda x: x > 0)
duplicate_ids = int(w89n.loc[w89n.duplicated("idauniq", keep=False), "idauniq"].nunique())
w89n = w89n.sort_values(["idauniq", "outcome_wave"]).drop_duplicates("idauniq", keep="first")

outcome_core = []
for wave, filename in [(8, "wave_8_elsa_data_eul_v2.dta"), (9, "wave_9_elsa_data_eul_v2.dta")]:
    core = read_dta(filename, ["idauniq", "hemda", "hemdab"])
    core["treated_out"] = treatment_status(core)
    core["outcome_wave"] = wave
    outcome_core.append(core[["idauniq", "outcome_wave", "treated_out"]])
outcome_core = pd.concat(outcome_core, ignore_index=True)
w89n = w89n.merge(outcome_core, on=["idauniq", "outcome_wave"], how="left", validate="one_to_one")

# Harmonized baseline covariates and public-use urban/rural indicator.
harmonized = read_dta(
    "gh_elsa_h.dta",
    [
        "idauniq", "ragender", "raeducl", "r6agey", "r6mstat", "r6shlt", "r6smoken", "r6drink",
        "r6diabe", "r6stroke", "r6hrtatte", "r6mdactx_e", "r6wtresp", "r6lwtresp",
    ],
)
harmonized["age_w6"] = plausible(harmonized.pop("r6agey"), 18, 110)
gender = allowed(harmonized.pop("ragender"), [1, 2])
harmonized["female"] = indicator_from_value(gender, [2])
harmonized["education"] = allowed(harmonized.pop("raeducl"), [1, 2, 3])
marital = allowed(harmonized.pop("r6mstat"), [1, 2, 3, 4, 5, 7, 8])
harmonized["married_partnered"] = indicator_from_value(marital, [1, 2, 3])
harmonized["self_health"] = allowed(harmonized.pop("r6shlt"), [1, 2, 3, 4, 5])
harmonized["current_smoker"] = allowed(harmonized.pop("r6smoken"), [0, 1])
harmonized["current_drinker"] = allowed(harmonized.pop("r6drink"), [0, 1])
harmonized["diabetes"] = allowed(harmonized.pop("r6diabe"), [0, 1])
stroke = allowed(harmonized.pop("r6stroke"), [0, 1])
heart_attack = allowed(harmonized.pop("r6hrtatte"), [0, 1])
harmonized["cvd_history"] = pd.concat([stroke, heart_attack], axis=1).max(axis=1, skipna=False)
harmonized["physical_activity"] = allowed(harmonized.pop("r6mdactx_e"), [2, 3, 4, 5])
harmonized["core_weight_w6"] = pd.to_numeric(harmonized.pop("r6wtresp"), errors="coerce").where(lambda x: x > 0)
harmonized["long_weight_w6"] = pd.to_numeric(harmonized.pop("r6lwtresp"), errors="coerce").where(lambda x: x > 0)

geography = read_dta("elsa_geog_urindewr_2011_eul.dta", ["idauniq", "w6_urindewr_2011"])
urban_rural = allowed(geography.pop("w6_urindewr_2011"), [1, 2])
geography["rural"] = indicator_from_value(urban_rural, [2])

# W6 target population with cross-wave exposure and follow-up information.
wide = (
    w6n.merge(w6c, on="idauniq", how="left", validate="one_to_one")
    .merge(harmonized, on="idauniq", how="left", validate="one_to_one")
    .merge(geography, on="idauniq", how="left", validate="one_to_one")
    .merge(w2[["idauniq", "w2_nurse_present", "bmi_w2"]], on="idauniq", how="left", validate="one_to_one")
    .merge(w4[["idauniq", "bmi_w4"]], on="idauniq", how="left", validate="one_to_one")
    .merge(w89n, on="idauniq", how="left", validate="one_to_one")
)
wide["person_id"] = pd.to_numeric(wide["idauniq"], errors="coerce").astype("Int64").astype("string")
wide["state5_w6"] = five_state(wide["sbp_w6"], wide["dbp_w6"], wide["treated_w6"])
wide["state5_out"] = five_state(wide["sbp_out"], wide["dbp_out"], wide["treated_out"])
wide["state5_accaha_w6"] = accaha_state(wide["sbp_w6"], wide["dbp_w6"], wide["treated_w6"])
wide["state5_accaha_out"] = accaha_state(wide["sbp_out"], wide["dbp_out"], wide["treated_out"])

bmi_values = wide[["bmi_w2", "bmi_w4", "bmi_w6"]]
wide["exposure_complete"] = bmi_values.notna().all(axis=1).astype(int)
wide["bmi_mean"] = bmi_values.mean(axis=1)
wide["bmi_sd"] = bmi_values.std(axis=1, ddof=1)
wide["bmi_cv"] = 100 * wide["bmi_sd"] / wide["bmi_mean"]
wide["bmi_arv"] = (abs(wide["bmi_w4"] - wide["bmi_w2"]) + abs(wide["bmi_w6"] - wide["bmi_w4"])) / 2
wide["bmi_auc"] = (wide["bmi_w2"] + wide["bmi_w4"]) / 2 * 4 + (wide["bmi_w4"] + wide["bmi_w6"]) / 2 * 4
wide["bmi_twmean"] = wide["bmi_auc"] / 8
excess = (bmi_values - 25).clip(lower=0)
wide["excess_bmi25"] = (excess["bmi_w2"] + excess["bmi_w4"]) / 2 * 4 + (excess["bmi_w4"] + excess["bmi_w6"]) / 2 * 4
wide["bmi_change_pct"] = (wide["bmi_w6"] - wide["bmi_w2"]) / wide["bmi_w2"] * 100
wide["bmi_change_cat_code"] = np.select([wide["bmi_change_pct"].lt(-5), wide["bmi_change_pct"].ge(5)], [-1, 1], default=0)
wide["age10"] = wide["age_w6"] / 10
wide["age65"] = wide["age_w6"].ge(65).astype(int)
wide["obese_w6"] = wide["bmi_w6"].ge(30).astype(int)

# Target requires actual W2 nurse participation, objective W6 BMI/state and a valid W6 nurse weight.
target_mask = (
    wide["w2_nurse_present"].eq(1)
    & wide["age_w6"].ge(50)
    & wide["bmi_w6"].notna()
    & wide["state5_w6"].notna()
    & wide["w6nurwt"].gt(0)
)
target = wide.loc[target_mask].copy()

ipw_continuous = ["age_w6", "self_health", "bmi_w6"]
ipw_categorical = [
    "female", "education", "married_partnered", "current_smoker", "current_drinker", "diabetes",
    "cvd_history", "physical_activity", "region", "rural", "state5_w6",
]
_, sw_exposure, diag_exposure = fit_sw(
    "complete_W2_W4_W6_BMI_history", target, target["exposure_complete"], ipw_continuous, ipw_categorical
)
target["sw_exposure"] = sw_exposure
sw_exposure_map = target.set_index("person_id")["sw_exposure"]

cohort = target.loc[target["exposure_complete"].eq(1)].copy()
cohort["sw_exposure"] = cohort["person_id"].map(sw_exposure_map)
vim_source = cohort.loc[cohort["bmi_sd"].gt(0), ["bmi_mean", "bmi_sd"]]
vim_beta = float(np.polyfit(np.log(vim_source["bmi_mean"]), np.log(vim_source["bmi_sd"]), 1)[0])
vim_scale = float(cohort["bmi_mean"].mean() ** vim_beta)
cohort["bmi_vim"] = cohort["bmi_sd"] / (cohort["bmi_mean"] ** vim_beta) * vim_scale

cohort["state_observed_out"] = cohort["state5_out"].notna().astype(int)
follow_continuous = ["age_w6", "self_health", "bmi_mean", "bmi_vim", "excess_bmi25"]
follow_categorical = ipw_categorical
_, sw_outcome, diag_outcome = fit_sw(
    "W8_or_W9_five_state_observed", cohort, cohort["state_observed_out"], follow_continuous, follow_categorical
)
cohort["sw_outcome"] = sw_outcome
cohort["base_weight"] = cohort["w6nurwt"] * cohort["sw_exposure"]
cohort["model_weight"] = cohort["base_weight"] * cohort["sw_outcome"]
cohort["outcome_nurse_weight"] = cohort["w89nurwt_g"] * cohort["sw_exposure"]

exposure_vars = [
    "excess_bmi25", "bmi_auc", "bmi_twmean", "bmi_mean", "bmi_vim", "bmi_arv", "bmi_cv", "bmi_sd", "bmi_change_pct"
]
scaling = []
for variable in exposure_vars:
    mean = float(cohort[variable].mean())
    sd = float(cohort[variable].std(ddof=0))
    cohort[f"{variable}_z"] = (cohort[variable] - mean) / sd
    scaling.append({"variable": variable, "mean": mean, "sd": sd})
cohort["burden_q"] = pd.qcut(cohort["excess_bmi25"], 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop").astype("string")
cohort["vim_q"] = pd.qcut(cohort["bmi_vim"], 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop").astype("string")

# One panel interval per person; W8 and W9 health visits have different durations.
intervals = cohort.copy()
intervals["origin"] = intervals["state5_w6"]
intervals["destination"] = intervals["state5_out"]
intervals["origin_accaha"] = intervals["state5_accaha_w6"]
intervals["destination_accaha"] = intervals["state5_accaha_out"]
intervals["start_year"] = 2012
intervals["end_year"] = intervals["outcome_wave"].map({8: 2016, 9: 2018})
intervals["interval_years"] = intervals["end_year"] - intervals["start_year"]
intervals["interval"] = intervals["outcome_wave"].map({8: "W6-W8", 9: "W6-W9"})
intervals["age_start10"] = intervals["age_w6"] / 10

def add_endpoints(frame: pd.DataFrame, origin: str, destination: str, suffix: str = "") -> None:
    o, d = frame[origin], frame[destination]
    frame[f"hypertension_onset{suffix}"] = np.where(o.isin([1, 2]) & d.notna(), d.isin([3, 4, 5]).astype(float), np.nan)
    frame[f"treatment_initiation{suffix}"] = np.where(o.eq(3) & d.notna(), d.isin([4, 5]).astype(float), np.nan)
    frame[f"untreated_bp_improvement{suffix}"] = np.where(o.eq(3) & d.isin([1, 2, 3]), d.isin([1, 2]).astype(float), np.nan)
    frame[f"control_loss{suffix}"] = np.where(o.eq(4) & d.isin([4, 5]), d.eq(5).astype(float), np.nan)
    frame[f"control_achievement{suffix}"] = np.where(o.eq(5) & d.isin([4, 5]), d.eq(4).astype(float), np.nan)


add_endpoints(intervals, "origin", "destination")
add_endpoints(intervals, "origin_accaha", "destination_accaha", "_accaha")

endpoint_names = [
    "hypertension_onset", "treatment_initiation", "untreated_bp_improvement", "control_loss", "control_achievement"
]
event_rows = []
for endpoint in endpoint_names:
    eligible = intervals[endpoint].notna()
    event_rows.append(
        {
            "endpoint": endpoint,
            "eligible_intervals": int(eligible.sum()),
            "events": int(intervals.loc[eligible, endpoint].sum()),
            "nonevents": int(eligible.sum() - intervals.loc[eligible, endpoint].sum()),
            "participants": int(intervals.loc[eligible, "person_id"].nunique()),
            "weighted_eligible": int((eligible & intervals["model_weight"].gt(0)).sum()),
        }
    )
event_gate = pd.DataFrame(event_rows)

transition_rows = []
observed = intervals.loc[intervals["origin"].notna() & intervals["destination"].notna()]
for interval, subset in observed.groupby("interval", observed=True):
    matrix = pd.crosstab(subset["origin"].astype(int), subset["destination"].astype(int))
    for origin in range(1, 6):
        for destination in range(1, 6):
            transition_rows.append(
                {
                    "interval": interval,
                    "origin": origin,
                    "destination": destination,
                    "n": int(matrix.loc[origin, destination]) if origin in matrix.index and destination in matrix.columns else 0,
                }
            )
transition_counts = pd.DataFrame(transition_rows)

flow = pd.DataFrame(
    [
        {"step": "W6 nurse participants", "n": int(len(w6n))},
        {"step": "W2 nurse participant + age >=50 + valid W6 BMI/state/weight", "n": int(len(target))},
        {"step": "Complete W2/W4/W6 measured BMI exposure", "n": int(len(cohort))},
        {"step": "Observed W8/9 five-state outcome", "n": int(cohort["state_observed_out"].sum())},
        {"step": "Positive primary analysis weight", "n": int(cohort["model_weight"].gt(0).sum())},
    ]
)

covariates = [
    "age_w6", "female", "education", "married_partnered", "current_smoker", "current_drinker", "diabetes",
    "cvd_history", "self_health", "physical_activity", "region", "rural", "idahhw6", "w6nurwt",
]
covariate_missingness = pd.DataFrame(
    [
        {"variable": variable, "missing_n": int(cohort[variable].isna().sum()), "missing_percent": float(cohort[variable].isna().mean() * 100)}
        for variable in covariates
    ]
)
exposure_summary = cohort[exposure_vars].describe(percentiles=[0.25, 0.5, 0.75]).T.reset_index(names="variable")
exposure_correlations = cohort[exposure_vars].corr().reset_index(names="variable")

baseline_columns = [
    "person_id", "idauniq", "idahhw6", "age_w6", "female", "education", "married_partnered", "current_smoker",
    "current_drinker", "diabetes", "cvd_history", "self_health", "physical_activity", "region", "rural",
    "bmi_w2", "bmi_w4", "bmi_w6", "bmi_mean", "bmi_twmean", "bmi_auc", "excess_bmi25", "bmi_sd", "bmi_cv",
    "bmi_arv", "bmi_vim", "bmi_change_pct", "bmi_change_cat_code", "state5_w6", "state5_out", "outcome_wave",
    "base_weight", "model_weight", "outcome_nurse_weight", "burden_q", "vim_q",
] + [f"{variable}_z" for variable in exposure_vars]

cohort[baseline_columns].to_csv(OUT / "elsa_baseline_w6.csv", index=False, encoding="utf-8")
intervals.to_csv(OUT / "elsa_intervals.csv", index=False, encoding="utf-8")
flow.to_csv(OUT / "elsa_cohort_flow.csv", index=False, encoding="utf-8-sig")
pd.DataFrame([diag_exposure, diag_outcome]).to_csv(OUT / "elsa_ipw_diagnostics.csv", index=False, encoding="utf-8-sig")
event_gate.to_csv(OUT / "elsa_endpoint_event_gate.csv", index=False, encoding="utf-8-sig")
transition_counts.to_csv(OUT / "elsa_transition_counts.csv", index=False, encoding="utf-8-sig")
covariate_missingness.to_csv(OUT / "elsa_covariate_missingness.csv", index=False, encoding="utf-8-sig")
exposure_summary.to_csv(OUT / "elsa_exposure_summary.csv", index=False, encoding="utf-8-sig")
exposure_correlations.to_csv(OUT / "elsa_exposure_correlations.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(scaling).to_csv(OUT / "elsa_exposure_scaling.csv", index=False, encoding="utf-8-sig")

summary = {
    "exposure_target_n": int(len(target)),
    "baseline_exposure_complete_n": int(len(cohort)),
    "outcome_state_observed_n": int(cohort["state_observed_out"].sum()),
    "primary_weight_positive_n": int(cohort["model_weight"].gt(0).sum()),
    "vim_beta": vim_beta,
    "duplicate_w8_w9_ids_retained_at_first_visit": duplicate_ids,
    "event_gate": event_gate.to_dict("records"),
    "design_note": "W6 nurse weight multiplied by stabilized exposure-history and outcome-observation weights; household-clustered inference.",
}
(OUT / "elsa_replication_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps(summary, indent=2, ensure_ascii=False))
