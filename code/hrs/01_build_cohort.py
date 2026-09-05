from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadstat
import scipy.linalg
import statsmodels.api as sm


REPO_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = Path(os.getenv("BMI_BP_OUTPUT_DIR", REPO_ROOT / "outputs"))
DATA_ROOT = Path(os.getenv("HRS_DATA_DIR", REPO_ROOT / "data" / "hrs"))
RAND_LONG = DATA_ROOT / "01_rand_longitudinal" / "extracted" / "randhrs1992_2022v1_STATA" / "randhrs1992_2022v1.dta"
TRACKER = DATA_ROOT / "04_tracker" / "extracted" / "trk2022tr_r.dta"
REGION = DATA_ROOT / "05_region_mobility" / "extracted" / "HRSXREGION22.dta"
FAT = {
    2014: DATA_ROOT / "02_rand_fat" / "2014" / "extracted" / "h14f2b.dta",
    2018: DATA_ROOT / "02_rand_fat" / "2018" / "extracted" / "h18f2c.dta",
    2022: DATA_ROOT / "02_rand_fat" / "2022" / "extracted" / "h22e3a.dta",
}
OUT = OUTPUT_ROOT / "hrs_formal_cohort"
OUT.mkdir(parents=True, exist_ok=True)

WAVE = {2014: 12, 2018: 14, 2022: 16}
PREFIX = {2014: "O", 2018: "Q", 2022: "S"}
MED_VARS = {
    2014: ("oc005", "oc006"),
    2018: ("qc005", "qc006"),
    2022: ("sc005", "sc006"),
}
STATE_NAMES = {
    1: "Normal untreated",
    2: "Elevated untreated",
    3: "Untreated hypertension",
    4: "Treated controlled",
    5: "Treated uncontrolled",
    6: "Death",
}


def normalize_id(frame: pd.DataFrame) -> pd.Series:
    hhid = frame["hhid"].astype("string").str.strip().str.replace(r"\.0$", "", regex=True).str.zfill(6)
    pn = frame["pn"].astype("string").str.strip().str.replace(r"\.0$", "", regex=True).str.zfill(3)
    return hhid + pn


def numeric(frame: pd.DataFrame, columns: list[str]) -> None:
    for column in columns:
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")


def current_medication(raw_med: pd.Series, diagnosis: pd.Series) -> tuple[pd.Series, pd.Series]:
    med = pd.Series(pd.NA, index=raw_med.index, dtype="Int64")
    source = pd.Series("unresolved", index=raw_med.index, dtype="string")
    med.loc[raw_med.eq(1)] = 1
    source.loc[raw_med.eq(1)] = "observed_yes"
    med.loc[raw_med.isin([4, 5, 6])] = 0
    source.loc[raw_med.isin([4, 5, 6])] = "observed_no"
    infer_no = raw_med.isna() & diagnosis.isin([4, 5, 6])
    med.loc[infer_no] = 0
    source.loc[infer_no] = "inferred_no_from_diagnosis_skip"
    source.loc[raw_med.isin([-8, 8, 9])] = "item_nonresponse"
    source.loc[raw_med.isna() & diagnosis.isin([1, 3])] = "diagnosed_but_med_blank"
    source.loc[raw_med.isna() & diagnosis.isna()] = "both_blank"
    return med, source


def five_state(sbp: pd.Series, dbp: pd.Series, med: pd.Series) -> pd.Series:
    state = pd.Series(pd.NA, index=sbp.index, dtype="Int64")
    valid = sbp.notna() & dbp.notna() & med.notna()
    normal = (sbp < 120) & (dbp < 80)
    below_htn = (sbp < 140) & (dbp < 90)
    htn = (sbp >= 140) | (dbp >= 90)
    state.loc[valid & med.eq(0) & normal] = 1
    state.loc[valid & med.eq(0) & below_htn & ~normal] = 2
    state.loc[valid & med.eq(0) & htn] = 3
    state.loc[valid & med.eq(1) & below_htn] = 4
    state.loc[valid & med.eq(1) & htn] = 5
    return state


def region4(division: pd.Series) -> pd.Series:
    out = pd.Series(pd.NA, index=division.index, dtype="Int64")
    out.loc[division.isin([1, 2])] = 1
    out.loc[division.isin([3, 4])] = 2
    out.loc[division.isin([5, 6, 7])] = 3
    out.loc[division.isin([8, 9])] = 4
    return out


def effective_sample_size(weight: pd.Series) -> float:
    weight = weight.dropna().astype(float)
    if weight.empty or float((weight**2).sum()) == 0:
        return np.nan
    return float(weight.sum() ** 2 / (weight**2).sum())


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
    # Some wave-specific survivor subsets make otherwise distinct missingness or
    # category indicators exactly collinear. Keep an independent column basis so
    # the GLM coefficients and predicted probabilities are uniquely determined.
    if np.linalg.matrix_rank(design.to_numpy()) < design.shape[1]:
        _, _, pivot = scipy.linalg.qr(design.to_numpy(), mode="economic", pivoting=True)
        rank = np.linalg.matrix_rank(design.to_numpy())
        keep = sorted(pivot[:rank])
        design = design.iloc[:, keep]
    return design


def fit_stabilized_weight(
    name: str,
    frame: pd.DataFrame,
    outcome: pd.Series,
    continuous: list[str],
    categorical: list[str],
) -> tuple[pd.Series, pd.Series, dict, pd.Series]:
    y = pd.to_numeric(outcome, errors="coerce")
    eligible = y.isin([0, 1])
    data = frame.loc[eligible].copy()
    y_fit = y.loc[eligible].astype(float)
    design = build_design(data, continuous, categorical)
    model_type = "GLM"
    try:
        fitted = sm.GLM(y_fit, design, family=sm.families.Binomial()).fit(maxiter=200, disp=0)
    except Exception:
        fitted = sm.GLM(y_fit, design, family=sm.families.Binomial()).fit_regularized(alpha=1e-6, L1_wt=0.0, maxiter=500)
        model_type = "ridge_regularized_GLM"
    probability = pd.Series(np.nan, index=frame.index, dtype=float)
    probability.loc[eligible] = np.asarray(fitted.predict(design), dtype=float)
    probability = probability.clip(0.02, 0.98)
    numerator = float(y_fit.mean())
    sw = pd.Series(np.nan, index=frame.index, dtype=float)
    sw.loc[eligible & y.eq(1)] = numerator / probability.loc[eligible & y.eq(1)]
    selected = sw.dropna()
    if selected.empty:
        truncated = sw.copy()
        lower = upper = np.nan
    else:
        lower, upper = selected.quantile([0.01, 0.99]).tolist()
        truncated = sw.clip(lower, upper)
    selected_t = truncated.dropna()
    diag = {
        "model": name,
        "model_type": model_type,
        "eligible_n": int(eligible.sum()),
        "selected_n": int((eligible & y.eq(1)).sum()),
        "selected_rate": numerator,
        "probability_min": float(probability.loc[eligible].min()),
        "probability_p01": float(probability.loc[eligible].quantile(0.01)),
        "probability_median": float(probability.loc[eligible].median()),
        "probability_p99": float(probability.loc[eligible].quantile(0.99)),
        "probability_max": float(probability.loc[eligible].max()),
        "sw_min": float(selected.min()) if not selected.empty else np.nan,
        "sw_p01": float(selected.quantile(0.01)) if not selected.empty else np.nan,
        "sw_median": float(selected.median()) if not selected.empty else np.nan,
        "sw_p99": float(selected.quantile(0.99)) if not selected.empty else np.nan,
        "sw_max": float(selected.max()) if not selected.empty else np.nan,
        "truncate_low": lower,
        "truncate_high": upper,
        "ess_untruncated": effective_sample_size(selected),
        "ess_truncated": effective_sample_size(selected_t),
        "design_columns": int(design.shape[1]),
    }
    coefficients = pd.Series(np.asarray(fitted.params), index=design.columns, name=name)
    return probability, truncated, diag, coefficients


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce")
    weights = pd.to_numeric(weights, errors="coerce")
    ok = values.notna() & weights.notna() & weights.gt(0)
    if not ok.any():
        return np.nan
    return float(np.average(values.loc[ok], weights=weights.loc[ok]))


def balance_rows(
    model: str,
    target: pd.DataFrame,
    selected: pd.Series,
    weight: pd.Series,
    variables: list[str],
) -> list[dict]:
    rows = []
    target_weight = pd.Series(1.0, index=target.index)
    for variable in variables:
        values = pd.to_numeric(target[variable], errors="coerce")
        target_mean = weighted_mean(values, target_weight)
        target_sd = float(values.std(ddof=0))
        if not np.isfinite(target_sd) or target_sd == 0:
            target_sd = 1.0
        before_mean = weighted_mean(values.loc[selected], pd.Series(1.0, index=target.index).loc[selected])
        after_mean = weighted_mean(values.loc[selected], weight.loc[selected])
        rows.append(
            {
                "model": model,
                "variable": variable,
                "target_mean": target_mean,
                "selected_mean_before": before_mean,
                "selected_mean_after": after_mean,
                "abs_smd_before": abs(before_mean - target_mean) / target_sd,
                "abs_smd_after": abs(after_mean - target_mean) / target_sd,
            }
        )
    return rows


# RAND longitudinal variables used for exposure construction and baseline IPW predictors.
long_cols = [
    "hhid", "pn", "hacohort", "ragender", "raracem", "raeduc", "raedyrs",
    "r12agey_e", "r12mstat", "r12shlt", "r12smoken", "r12drink",
    "r12diabe", "r12hearte", "r12stroke", "r12cesd", "r12lbrf", "r12work",
    "h12atotb", "h12itot",
    "r8bmi", "r10bmi", "r12bmi",
    "r8pmbmi", "r10pmbmi", "r12pmbmi",
    "r12bpsys", "r12bpdia", "r14bpsys", "r14bpdia", "r16bpsys", "r16bpdia",
]
wide, _ = pyreadstat.read_dta(RAND_LONG, usecols=long_cols, apply_value_formats=False)
wide.columns = [c.lower() for c in wide.columns]
numeric(wide, [c for c in wide.columns if c not in {"hhid", "pn"}])
wide["person_id"] = normalize_id(wide)
if wide["person_id"].duplicated().any():
    raise RuntimeError("Duplicate IDs in RAND longitudinal file")

# Current antihypertensive medication from version-matched Fat Files.
for year, path in FAT.items():
    diagnosis, medication = MED_VARS[year]
    fat, _ = pyreadstat.read_dta(path, usecols=["hhid", "pn", diagnosis, medication], apply_value_formats=False)
    fat.columns = [c.lower() for c in fat.columns]
    numeric(fat, [diagnosis, medication])
    fat["person_id"] = normalize_id(fat)
    med_current, med_source = current_medication(fat[medication], fat[diagnosis])
    fat[f"med_current_{year}"] = med_current
    fat[f"med_source_{year}"] = med_source
    fat[f"in_fat_{year}"] = 1
    wide = wide.merge(
        fat[["person_id", diagnosis, medication, f"med_current_{year}", f"med_source_{year}", f"in_fat_{year}"]],
        on="person_id", how="left", validate="one_to_one"
    )

# Tracker merge.
tracker_cols = ["HHID", "PN", "EFTFASSIGN", "SECU", "STRATUM", "WTCOHORT"]
for prefix in ["K", "M", "O", "Q", "S"]:
    tracker_cols += [f"{prefix}IWWAVE", f"{prefix}IWTYPE", f"{prefix}ALIVE", f"{prefix}WGTR", f"{prefix}PMWGTR"]
tracker, _ = pyreadstat.read_dta(TRACKER, usecols=tracker_cols, apply_value_formats=False)
tracker.columns = [c.lower() for c in tracker.columns]
tracker["person_id"] = normalize_id(tracker)
numeric(tracker, [c for c in tracker.columns if c not in {"hhid", "pn", "person_id"}])
if tracker["person_id"].duplicated().any():
    raise RuntimeError("Duplicate IDs in Tracker")
tracker_keep = [c for c in tracker.columns if c not in {"hhid", "pn"}]
wide = wide.merge(tracker[tracker_keep], on="person_id", how="left", validate="one_to_one")

# Public region/division and urbanicity data.
region_cols = ["HHID", "PN", "REGION14", "REGION18", "REGION22", "BEALE2023_14", "BEALE2023_18", "BEALE2023_22"]
region, _ = pyreadstat.read_dta(REGION, usecols=region_cols, apply_value_formats=False)
region.columns = [c.lower() for c in region.columns]
region["person_id"] = normalize_id(region)
numeric(region, [c for c in region.columns if c not in {"hhid", "pn", "person_id"}])
if region["person_id"].duplicated().any():
    raise RuntimeError("Duplicate IDs in Region file")
wide = wide.merge(region.drop(columns=["hhid", "pn"]), on="person_id", how="left", validate="one_to_one")

# Valid measurements and five-state outcomes.
for column in ["r8pmbmi", "r10pmbmi", "r12pmbmi", "r8bmi", "r10bmi", "r12bmi"]:
    wide[column] = wide[column].where(wide[column].between(12, 60))
for year, wave in WAVE.items():
    wide[f"sbp_{year}"] = wide[f"r{wave}bpsys"].where(wide[f"r{wave}bpsys"].between(70, 250))
    wide[f"dbp_{year}"] = wide[f"r{wave}bpdia"].where(wide[f"r{wave}bpdia"].between(40, 150))
    wide[f"state5_{year}"] = five_state(wide[f"sbp_{year}"], wide[f"dbp_{year}"], wide[f"med_current_{year}"])

# Exposure metrics use three objective BMI measurements over a common eight-year window.
objective = wide[["r8pmbmi", "r10pmbmi", "r12pmbmi"]].to_numpy(dtype=float)
years = np.array([2006.0, 2010.0, 2014.0])
exposure_complete = np.isfinite(objective).all(axis=1)
wide["exposure_complete"] = exposure_complete.astype(int)
wide["bmi_mean"] = np.where(exposure_complete, objective.mean(axis=1), np.nan)
wide["bmi_auc"] = np.where(exposure_complete, np.trapezoid(objective, x=years, axis=1), np.nan)
wide["bmi_twmean"] = wide["bmi_auc"] / 8.0
wide["excess_bmi25"] = np.where(exposure_complete, np.trapezoid(np.maximum(objective - 25, 0), x=years, axis=1), np.nan)
wide["bmi_sd"] = np.where(exposure_complete, objective.std(axis=1, ddof=1), np.nan)
wide["bmi_cv"] = wide["bmi_sd"] / wide["bmi_mean"] * 100
wide["bmi_arv"] = np.where(exposure_complete, np.abs(np.diff(objective, axis=1)).mean(axis=1), np.nan)
wide["bmi_self_mean"] = wide[["r8bmi", "r10bmi", "r12bmi"]].mean(axis=1, skipna=True)

initial_candidate = wide["r12agey_e"].ge(50) & wide["exposure_complete"].eq(1) & wide["state5_2014"].notna()
vim_fit = initial_candidate & wide["bmi_mean"].gt(0) & wide["bmi_sd"].gt(0)
vim_beta = float(np.polyfit(np.log(wide.loc[vim_fit, "bmi_mean"]), np.log(wide.loc[vim_fit, "bmi_sd"]), 1)[0])
wide["bmi_vim"] = wide["bmi_sd"] / np.power(wide["bmi_mean"], vim_beta) * np.power(wide.loc[initial_candidate, "bmi_mean"].mean(), vim_beta)

# Baseline region and socioeconomic transformations for weighting models.
wide["region4_2014"] = region4(wide["region14"])
wide["urban3_2014"] = wide["beale2023_14"].where(wide["beale2023_14"].isin([1, 2, 3]))
wide["assets_asinh"] = np.arcsinh(wide["h12atotb"] / 10000.0)
wide["income_asinh"] = np.arcsinh(wide["h12itot"] / 10000.0)
wide["female"] = wide["ragender"].eq(2).astype(int)
wide["age2014_sq"] = wide["r12agey_e"] ** 2

# Tracker-derived status and death state.
for year in [2014, 2018, 2022]:
    prefix = PREFIX[year].lower()
    wide[f"core_obtained_{year}"] = wide[f"{prefix}iwtype"].eq(1).astype(int)
    wide[f"alive_presumed_{year}"] = wide[f"{prefix}alive"].isin([1, 2]).astype(int)
    wide[f"died_by_{year}"] = wide[f"{prefix}alive"].isin([5, 6]).astype(int)
    wide[f"state_observed_{year}"] = wide[f"state5_{year}"].notna().astype(int)
    death_state = wide[f"state5_{year}"].astype("Float64")
    death_state.loc[wide[f"died_by_{year}"].eq(1)] = 6
    wide[f"state6_{year}"] = death_state

baseline = wide["r12agey_e"].ge(50) & wide["exposure_complete"].eq(1) & wide["state5_2014"].notna()
cohort = wide.loc[baseline].copy()
if len(cohort) != 3740:
    raise RuntimeError(f"Expected 3,740 baseline participants, found {len(cohort)}")

# Candidate stabilized weight for having the complete three-measure objective BMI history.
exposure_target_mask = (
    wide["r12agey_e"].ge(50)
    & wide["state5_2014"].notna()
    & wide["eftfassign"].eq(1)
    & wide["kiwtype"].ne(99)
)
exposure_target = wide.loc[exposure_target_mask].copy()
exposure_outcome = exposure_target["exposure_complete"]
continuous_exposure = ["r12agey_e", "age2014_sq", "raedyrs", "r12shlt", "r12cesd", "bmi_self_mean", "assets_asinh", "income_asinh"]
categorical_common = ["ragender", "raracem", "r12mstat", "r12smoken", "r12drink", "r12diabe", "r12hearte", "r12stroke", "region4_2014", "urban3_2014"]
p_exp, sw_exp, diag_exp, coef_exp = fit_stabilized_weight(
    "complete_objective_bmi_history", exposure_target, exposure_outcome, continuous_exposure, categorical_common
)
exposure_target["p_exposure_complete"] = p_exp
exposure_target["sw_exposure"] = sw_exp
exp_weight_map = exposure_target.set_index("person_id")["sw_exposure"]
cohort["sw_exposure"] = cohort["person_id"].map(exp_weight_map)

# Baseline physical-measurement weight, normalized within the analysis cohort.
cohort["opm_weight"] = cohort["opmwgtr"].where(cohort["opmwgtr"].gt(0))
cohort["opm_weight_norm"] = cohort["opm_weight"] / cohort["opm_weight"].mean()
cohort["base_weight_candidate"] = cohort["opm_weight_norm"] * cohort["sw_exposure"]
cohort["base_weight_candidate"] = cohort["base_weight_candidate"] / cohort["base_weight_candidate"].mean()

# Follow-up observation/censoring models.
continuous_followup = ["r12agey_e", "raedyrs", "r12shlt", "r12cesd", "bmi_mean", "bmi_vim", "assets_asinh", "income_asinh"]
categorical_followup = categorical_common + ["state5_2014"]

p18, sw18, diag18, coef18 = fit_stabilized_weight(
    "state_observed_2018_including_death_as_censoring",
    cohort,
    cohort["state_observed_2018"],
    continuous_followup,
    categorical_followup,
)
cohort["p_state2018"] = p18
cohort["sw_state2018"] = sw18

# Among survivors, separate response/measurement missingness from death.
survivor18_mask = cohort["alive_presumed_2018"].eq(1)
p18_surv, sw18_surv, diag18_surv, coef18_surv = fit_stabilized_weight(
    "state_observed_2018_among_survivors",
    cohort.loc[survivor18_mask],
    cohort.loc[survivor18_mask, "state_observed_2018"],
    continuous_followup,
    categorical_followup,
)
cohort["p_state2018_survivor"] = np.nan
cohort["sw_state2018_survivor"] = np.nan
cohort.loc[survivor18_mask, "p_state2018_survivor"] = p18_surv
cohort.loc[survivor18_mask, "sw_state2018_survivor"] = sw18_surv

# Conditional 2022 observation among those with an observed 2018 state.
cond22_mask = cohort["state_observed_2018"].eq(1)
cohort["state5_2018_for_model"] = cohort["state5_2018"]
continuous_2022 = continuous_followup
categorical_2022 = categorical_followup + ["state5_2018_for_model"]
p22, sw22, diag22, coef22 = fit_stabilized_weight(
    "state_observed_2022_conditional_on_2018_observed",
    cohort.loc[cond22_mask],
    cohort.loc[cond22_mask, "state_observed_2022"],
    continuous_2022,
    categorical_2022,
)
cohort["p_state2022_cond"] = np.nan
cohort["sw_state2022_cond"] = np.nan
cohort.loc[cond22_mask, "p_state2022_cond"] = p22
cohort.loc[cond22_mask, "sw_state2022_cond"] = sw22

survivor22_cond_mask = cond22_mask & cohort["alive_presumed_2022"].eq(1)
p22_surv, sw22_surv, diag22_surv, coef22_surv = fit_stabilized_weight(
    "state_observed_2022_among_survivors_conditional_on_2018_observed",
    cohort.loc[survivor22_cond_mask],
    cohort.loc[survivor22_cond_mask, "state_observed_2022"],
    continuous_2022,
    categorical_2022,
)
cohort["p_state2022_survivor"] = np.nan
cohort["sw_state2022_survivor"] = np.nan
cohort.loc[survivor22_cond_mask, "p_state2022_survivor"] = p22_surv
cohort.loc[survivor22_cond_mask, "sw_state2022_survivor"] = sw22_surv

cohort["ipcw_2018"] = cohort["sw_state2018"]
cohort["ipcw_2022_cumulative"] = cohort["sw_state2018"] * cohort["sw_state2022_cond"]
cohort["analysis_weight_2018"] = cohort["base_weight_candidate"] * cohort["ipcw_2018"]
cohort["analysis_weight_2022"] = cohort["base_weight_candidate"] * cohort["ipcw_2022_cumulative"]

# Missingness reason hierarchy.
missingness_rows = []
for year in [2018, 2022]:
    observed = cohort[f"state_observed_{year}"].eq(1)
    died = cohort[f"died_by_{year}"].eq(1)
    alive = cohort[f"alive_presumed_{year}"].eq(1)
    core = cohort[f"core_obtained_{year}"].eq(1)
    reasons = pd.Series("other_or_version_mismatch", index=cohort.index, dtype="string")
    reasons.loc[observed] = "state_observed"
    reasons.loc[~observed & died] = "dead_by_wave"
    reasons.loc[~observed & alive & ~core] = "alive_no_core_interview"
    reasons.loc[~observed & alive & core] = "core_interview_but_state_missing"
    for reason, count in reasons.value_counts().items():
        missingness_rows.append({"year": year, "reason": reason, "n": int(count), "percent": float(count / len(cohort))})
    cohort[f"missing_reason_{year}"] = reasons

# Wave-specific long-form analysis dataset.
long_rows = []
for year, wave in WAVE.items():
    prefix = PREFIX[year].lower()
    piece = pd.DataFrame(index=cohort.index)
    piece["person_id"] = cohort["person_id"]
    piece["hhid"] = cohort["hhid"].astype("string")
    piece["pn"] = cohort["pn"].astype("string")
    piece["year"] = year
    piece["wave"] = wave
    piece["age_2014"] = cohort["r12agey_e"]
    piece["gender"] = cohort["ragender"]
    piece["race"] = cohort["raracem"]
    piece["educ_years"] = cohort["raedyrs"]
    piece["bmi_2006"] = cohort["r8pmbmi"]
    piece["bmi_2010"] = cohort["r10pmbmi"]
    piece["bmi_2014"] = cohort["r12pmbmi"]
    for name in ["bmi_mean", "bmi_twmean", "bmi_auc", "excess_bmi25", "bmi_sd", "bmi_cv", "bmi_arv", "bmi_vim"]:
        piece[name] = cohort[name]
    piece["sbp"] = cohort[f"sbp_{year}"]
    piece["dbp"] = cohort[f"dbp_{year}"]
    piece["med_current"] = cohort[f"med_current_{year}"]
    piece["state5"] = cohort[f"state5_{year}"]
    piece["state6_death"] = cohort[f"state6_{year}"]
    piece["state_observed"] = cohort[f"state_observed_{year}"]
    piece["alive_presumed"] = cohort[f"alive_presumed_{year}"]
    piece["died_by_wave"] = cohort[f"died_by_{year}"]
    piece["core_obtained"] = cohort[f"core_obtained_{year}"]
    piece["iwtype"] = cohort[f"{prefix}iwtype"]
    piece["iwwave"] = cohort[f"{prefix}iwwave"]
    piece["alive_code"] = cohort[f"{prefix}alive"]
    piece["core_weight"] = cohort[f"{prefix}wgtr"]
    piece["pm_weight"] = cohort[f"{prefix}pmwgtr"]
    piece["region_division"] = cohort[f"region{str(year)[-2:]}"]
    piece["region4"] = region4(piece["region_division"])
    piece["urban3"] = cohort[f"beale2023_{str(year)[-2:]}"]
    piece["eftf_assign"] = cohort["eftfassign"]
    piece["secu"] = cohort["secu"]
    piece["stratum"] = cohort["stratum"]
    piece["opm_weight_norm"] = cohort["opm_weight_norm"]
    piece["sw_exposure"] = cohort["sw_exposure"]
    piece["base_weight_candidate"] = cohort["base_weight_candidate"]
    if year == 2014:
        piece["ipcw_wave"] = 1.0
        piece["analysis_weight_ipcw"] = cohort["base_weight_candidate"]
        piece["missing_reason"] = "baseline_observed"
    elif year == 2018:
        piece["ipcw_wave"] = cohort["ipcw_2018"]
        piece["analysis_weight_ipcw"] = cohort["analysis_weight_2018"]
        piece["missing_reason"] = cohort["missing_reason_2018"]
    else:
        piece["ipcw_wave"] = cohort["ipcw_2022_cumulative"]
        piece["analysis_weight_ipcw"] = cohort["analysis_weight_2022"]
        piece["missing_reason"] = cohort["missing_reason_2022"]
    long_rows.append(piece)
long = pd.concat(long_rows, ignore_index=True)

# Attrition/death flow and version-alignment checks.
flow_rows = []
for year in [2014, 2018, 2022]:
    flow_rows.append(
        {
            "year": year,
            "baseline_cohort": len(cohort),
            "core_obtained": int(cohort[f"core_obtained_{year}"].sum()),
            "alive_or_presumed": int(cohort[f"alive_presumed_{year}"].sum()),
            "dead_by_wave": int(cohort[f"died_by_{year}"].sum()),
            "five_state_observed": int(cohort[f"state_observed_{year}"].sum()),
            "positive_pm_weight": int(cohort[f"{PREFIX[year].lower()}pmwgtr"].gt(0).sum()),
        }
    )
attrition_flow = pd.DataFrame(flow_rows)

version_alignment = pd.DataFrame(
    [
        {
            "check": "2022 early Fat state observed but final Tracker core status not 1",
            "n": int((cohort["state_observed_2022"].eq(1) & ~cohort["siwtype"].eq(1)).sum()),
        },
        {
            "check": "2022 final Tracker core status 1 but five-state missing (primarily no physical measurement)",
            "n": int((cohort["siwtype"].eq(1) & cohort["state_observed_2022"].eq(0)).sum()),
        },
        {
            "check": "baseline cohort not matched to final Tracker",
            "n": int(cohort["eftfassign"].isna().sum()),
        },
        {
            "check": "baseline cohort not matched to Region file",
            "n": int(cohort["region14"].isna().sum()),
        },
    ]
)

weight_diagnostics = pd.DataFrame([diag_exp, diag18, diag18_surv, diag22, diag22_surv])
coefficients = pd.concat([coef_exp, coef18, coef18_surv, coef22, coef22_surv], axis=1).reset_index(names="term")

balance = []
balance_vars_exp = ["r12agey_e", "female", "raedyrs", "r12shlt", "bmi_self_mean", "r12diabe", "r12hearte", "r12stroke"]
balance += balance_rows(
    "complete_objective_bmi_history",
    exposure_target,
    exposure_target["exposure_complete"].eq(1),
    exposure_target["sw_exposure"],
    balance_vars_exp,
)
balance_vars_follow = ["r12agey_e", "female", "raedyrs", "r12shlt", "bmi_mean", "bmi_vim", "r12diabe", "r12hearte", "r12stroke"]
balance += balance_rows(
    "state_observed_2018_including_death_as_censoring",
    cohort,
    cohort["state_observed_2018"].eq(1),
    cohort["sw_state2018"],
    balance_vars_follow,
)
balance += balance_rows(
    "state_observed_2022_conditional_on_2018_observed",
    cohort.loc[cond22_mask],
    cohort.loc[cond22_mask, "state_observed_2022"].eq(1),
    cohort.loc[cond22_mask, "sw_state2022_cond"],
    balance_vars_follow,
)
balance = pd.DataFrame(balance)

# Compact codebook for the individual-level deliverables.
codebook_rows = [
    ["person_id", "HHID+PN person identifier"], ["year", "Observation year"],
    ["bmi_mean", "Arithmetic mean of measured BMI in 2006/2010/2014"],
    ["bmi_twmean", "Time-weighted mean measured BMI, 2006-2014"],
    ["bmi_auc", "Trapezoidal BMI area under curve, BMI-years"],
    ["excess_bmi25", "BMI-years above BMI 25 kg/m2"],
    ["bmi_vim", f"Variation independent of mean; beta={vim_beta:.6f}"],
    ["state5", "Five-state blood pressure/treatment status; codes 1-5"],
    ["state6_death", "Five-state status plus death as absorbing code 6"],
    ["region_division", "Public Census division of residence, codes 1-9"],
    ["urban3", "Beale 2023 collapsed urbanicity: 1 urban, 2 suburban, 3 ex-urban"],
    ["pm_weight", "Wave-specific HRS physical-measurement respondent weight"],
    ["opm_weight_norm", "Normalized 2014 physical-measurement weight"],
    ["sw_exposure", "Candidate stabilized weight for complete 3-wave measured-BMI history"],
    ["ipcw_wave", "Candidate stabilized inverse probability of wave observation/censoring weight"],
    ["analysis_weight_ipcw", "Candidate combined baseline and IPCW weight; not yet final for modeling"],
]
dataset_codebook = pd.DataFrame(codebook_rows, columns=["variable", "description"])

# Save aggregate audits.
attrition_flow.to_csv(OUT / "attrition_death_flow.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(missingness_rows).to_csv(OUT / "missingness_reasons.csv", index=False, encoding="utf-8-sig")
weight_diagnostics.to_csv(OUT / "ipw_diagnostics.csv", index=False, encoding="utf-8-sig")
coefficients.to_csv(OUT / "ipw_model_coefficients.csv", index=False, encoding="utf-8-sig")
balance.to_csv(OUT / "ipw_balance.csv", index=False, encoding="utf-8-sig")
version_alignment.to_csv(OUT / "version_alignment.csv", index=False, encoding="utf-8-sig")
dataset_codebook.to_csv(OUT / "analysis_dataset_codebook.csv", index=False, encoding="utf-8-sig")

# Stata-ready individual-level datasets. Convert nullable integers to floats for portability.
baseline_columns = [
    "person_id", "hhid", "pn", "r12agey_e", "ragender", "raracem", "raedyrs",
    "r12mstat", "r12shlt", "r12smoken", "r12drink", "r12diabe", "r12hearte", "r12stroke", "r12cesd",
    "r8pmbmi", "r10pmbmi", "r12pmbmi", "bmi_mean", "bmi_twmean", "bmi_auc", "excess_bmi25",
    "bmi_sd", "bmi_cv", "bmi_arv", "bmi_vim", "state5_2014", "state5_2018", "state5_2022",
    "state6_2014", "state6_2018", "state6_2022", "region14", "region18", "region22",
    "beale2023_14", "beale2023_18", "beale2023_22", "eftfassign", "secu", "stratum",
    "opmwgtr", "qpmwgtr", "spmwgtr", "sw_exposure", "sw_state2018", "sw_state2022_cond",
    "base_weight_candidate", "analysis_weight_2018", "analysis_weight_2022",
]
baseline_out = cohort[baseline_columns].copy()
for column in baseline_out.columns:
    if str(baseline_out[column].dtype) in {"Int64", "boolean"}:
        baseline_out[column] = baseline_out[column].astype(float)

for frame in [baseline_out, long]:
    for column in frame.columns:
        if str(frame[column].dtype) in {"Int64", "boolean"}:
            frame[column] = frame[column].astype(float)

baseline_out.to_stata(OUT / "hrs_o1_baseline_2014.dta", write_index=False, version=118)
long.to_stata(OUT / "hrs_o1_person_wave_2014_2022.dta", write_index=False, version=118)
baseline_out.to_csv(OUT / "hrs_o1_baseline_2014.csv", index=False, encoding="utf-8-sig")
long.to_csv(OUT / "hrs_o1_person_wave_2014_2022.csv", index=False, encoding="utf-8-sig")

# Summary and report.
diagnostic_summary = {
    "baseline_n": int(len(cohort)),
    "person_wave_rows": int(len(long)),
    "observed_state_2018": int(cohort["state_observed_2018"].sum()),
    "observed_state_2022": int(cohort["state_observed_2022"].sum()),
    "both_followups_observed": int((cohort["state_observed_2018"].eq(1) & cohort["state_observed_2022"].eq(1)).sum()),
    "dead_by_2018": int(cohort["died_by_2018"].sum()),
    "dead_by_2022": int(cohort["died_by_2022"].sum()),
    "positive_opm_weight": int(cohort["opmwgtr"].gt(0).sum()),
    "tracker_unmatched": int(cohort["eftfassign"].isna().sum()),
    "region_unmatched": int(cohort["region14"].isna().sum()),
    "vim_beta": vim_beta,
    "max_abs_smd_before": float(balance["abs_smd_before"].max()),
    "max_abs_smd_after": float(balance["abs_smd_after"].max()),
    "ipw_decision": "Candidate weights acceptable for sensitivity analysis; death-state model remains required.",
}
(OUT / "formal_cohort_summary.json").write_text(json.dumps(diagnostic_summary, ensure_ascii=False, indent=2), encoding="utf-8")

diag_table = weight_diagnostics[["model", "eligible_n", "selected_n", "selected_rate", "probability_min", "probability_max", "sw_p99", "sw_max", "ess_truncated"]]
report = f"""# HRS正式分析队列、失访死亡与权重预检

## 完成情况

已按HHID+PN合并RAND纵向文件、2014/2018/2022 Fat Files、2022 Final Tracker及公开Region/Mobility文件，并生成2014基线数据和2014/2018/2022人—波次长表。

## 正式O1队列

- 基线人数：{len(cohort):,}。
- 人—波次记录：{len(long):,}。
- 2018五级状态可观察：{diagnostic_summary['observed_state_2018']:,}。
- 2022五级状态可观察：{diagnostic_summary['observed_state_2022']:,}。
- 两次随访均观察：{diagnostic_summary['both_followups_observed']:,}。
- 2018前/截至2018已死亡：{diagnostic_summary['dead_by_2018']:,}。
- 截至2022已死亡：{diagnostic_summary['dead_by_2022']:,}。

## 权重

- Tracker提供EFTFASSIGN、SECU、STRATUM、核心访谈权重和体格测量权重。
- 基线对象中2014体格测量权重大于0者：{diagnostic_summary['positive_opm_weight']:,}/{len(cohort):,}。
- 已生成三类候选权重：重复客观BMI完整性权重、随访观察权重、死亡与失访合并的IPCW。
- 权重在1%和99%分位截尾，详细概率范围、极端值和有效样本量见ipw_diagnostics.csv。
- 最大绝对标准化差异由{diagnostic_summary['max_abs_smd_before']:.3f}降至{diagnostic_summary['max_abs_smd_after']:.3f}。

## 重要边界

IPW不能完全解决死亡带来的信息性选择，因此正式主模型应保留五级血压状态；死亡作为第六吸收状态至少作为关键敏感性分析。2022结局来自与RAND纵向文件版本一致的Early Fat；Final Tracker用于死亡和权重，版本不一致个案已在version_alignment.csv单列。

## 下一步

进入协变量逐项缺失率与编码核查，然后建立多状态Markov主模型和死亡六状态敏感性模型。
"""
(OUT / "formal_cohort_report.md").write_text(report, encoding="utf-8")

print(json.dumps(diagnostic_summary, ensure_ascii=False, indent=2))
print("\nATTRITION")
print(attrition_flow.to_string(index=False))
print("\nMISSINGNESS")
print(pd.DataFrame(missingness_rows).to_string(index=False))
print("\nIPW DIAGNOSTICS")
print(diag_table.to_string(index=False))
print("\nVERSION ALIGNMENT")
print(version_alignment.to_string(index=False))
print(f"output={OUT}")
