from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = Path(os.getenv("BMI_BP_OUTPUT_DIR", ROOT / "outputs"))
OUT = OUTPUT_ROOT / "dynamic_analysis"
OUT.mkdir(parents=True, exist_ok=True)


def numeric(frame: pd.DataFrame, columns: list[str]) -> None:
    for column in columns:
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")


def add_core_and_buffered_endpoints(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    numeric(frame, ["origin", "destination", "start_sbp", "start_dbp", "end_sbp", "end_dbp"])
    frame["origin_core3"] = np.select(
        [frame["origin"].eq(1), frame["origin"].eq(2), frame["origin"].isin([3, 4, 5])],
        [1.0, 2.0, 3.0],
        default=np.nan,
    )
    frame["destination_core3"] = np.select(
        [frame["destination"].eq(1), frame["destination"].eq(2), frame["destination"].isin([3, 4, 5])],
        [1.0, 2.0, 3.0],
        default=np.nan,
    )
    known_core = frame["origin_core3"].notna() & frame["destination_core3"].notna()
    frame["core3_hypertension_onset"] = np.where(
        known_core & frame["origin_core3"].isin([1, 2]), frame["destination_core3"].eq(3).astype(float), np.nan
    )
    frame["core3_hypertension_improvement"] = np.where(
        known_core & frame["origin_core3"].eq(3), frame["destination_core3"].isin([1, 2]).astype(float), np.nan
    )

    # A positive hypertension margin means that at least one BP component is
    # above the 140/90 threshold.  A positive lower-state margin means that
    # both components are below the threshold by the stated buffer.
    frame["start_htn_margin"] = np.maximum(frame["start_sbp"] - 140, frame["start_dbp"] - 90)
    frame["end_htn_margin"] = np.maximum(frame["end_sbp"] - 140, frame["end_dbp"] - 90)
    frame["start_lower_margin"] = np.minimum(140 - frame["start_sbp"], 90 - frame["start_dbp"])
    frame["end_lower_margin"] = np.minimum(140 - frame["end_sbp"], 90 - frame["end_dbp"])
    for buffer in [5, 10]:
        # Retain treatment-defined incident hypertension (states 4/5), because
        # it cannot be attributed to a borderline BP reading.  For untreated
        # incident hypertension, require the destination BP to clear the
        # diagnostic boundary by the prespecified buffer.  Ambiguous, near-
        # threshold destination readings are excluded rather than recoded as
        # non-events.
        onset_risk = frame["origin"].isin([1, 2]) & frame["start_lower_margin"].ge(buffer)
        onset_event = frame["destination"].isin([4, 5]) | (
            frame["destination"].eq(3) & frame["end_htn_margin"].ge(buffer)
        )
        onset_nonevent = frame["destination"].isin([1, 2])
        onset_observed = onset_risk & (onset_event | onset_nonevent)
        frame[f"buffer{buffer}_hypertension_onset"] = np.where(onset_observed, onset_event.astype(float), np.nan)

        improvement_risk = frame["origin"].eq(3) & frame["start_htn_margin"].ge(buffer)
        # States 1 and 2 are untreated by construction; treated destination
        # states and near-threshold apparent improvements are excluded rather
        # than counted as non-improvement.
        improvement_event = frame["destination"].isin([1, 2]) & frame["end_lower_margin"].ge(buffer)
        improvement_nonevent = frame["destination"].eq(3)
        improvement_observed = improvement_risk & (improvement_event | improvement_nonevent)
        frame[f"buffer{buffer}_untreated_improvement"] = np.where(
            improvement_observed, improvement_event.astype(float), np.nan
        )
    return frame


def build_hrs() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    base_path = OUTPUT_ROOT / "hrs_analysis"
    source_path = OUTPUT_ROOT / "hrs_formal_cohort"
    intervals = pd.read_csv(base_path / "analysis_intervals.csv", dtype={"person_id": "string"}, low_memory=False)
    baseline = pd.read_csv(base_path / "analysis_baseline.csv", dtype={"person_id": "string"}, low_memory=False)
    long = pd.read_csv(source_path / "hrs_o1_person_wave_2014_2022.csv", dtype={"person_id": "string"}, low_memory=False)
    numeric(long, ["year", "sbp", "dbp", "med_current", "state5", "died_by_wave", "cvd_history"])
    bp = long.pivot(index="person_id", columns="year", values=["sbp", "dbp", "med_current", "state5", "died_by_wave"])
    bp.columns = [f"{variable}_{int(year)}" for variable, year in bp.columns]
    bp = bp.reset_index()
    intervals = intervals.merge(bp, on="person_id", how="left")
    intervals["start_sbp"] = np.where(intervals["start_year"].eq(2014), intervals["sbp_2014"], intervals["sbp_2018"])
    intervals["start_dbp"] = np.where(intervals["start_year"].eq(2014), intervals["dbp_2014"], intervals["dbp_2018"])
    intervals["end_sbp"] = np.where(intervals["end_year"].eq(2018), intervals["sbp_2018"], intervals["sbp_2022"])
    intervals["end_dbp"] = np.where(intervals["end_year"].eq(2018), intervals["dbp_2018"], intervals["dbp_2022"])
    intervals = add_core_and_buffered_endpoints(intervals)
    intervals["cohort"] = "HRS"
    intervals["cluster"] = intervals["secu"]
    intervals["strata_value"] = intervals["stratum"]
    intervals["interval_years"] = intervals["end_year"] - intervals["start_year"]

    baseline = baseline.merge(bp, on="person_id", how="left", suffixes=("", "_long"))
    numeric(baseline, ["state5_2014", "state5_2018", "state5_2022", "analysis_weight_2022"])
    eligible = baseline["state5_2014"].eq(3) & baseline[["state5_2018", "state5_2022"]].notna().all(axis=1)
    baseline["initial_improvement"] = np.where(
        eligible, baseline["state5_2018"].isin([1, 2]).astype(float), np.nan
    )
    baseline["sustained_improvement"] = np.where(
        eligible,
        (baseline["state5_2018"].isin([1, 2]) & baseline["state5_2022"].isin([1, 2])).astype(float),
        np.nan,
    )
    baseline["maintenance_among_initial_improvers"] = np.where(
        eligible & baseline["state5_2018"].isin([1, 2]), baseline["state5_2022"].isin([1, 2]).astype(float), np.nan
    )
    baseline["cohort"] = "HRS"
    baseline["model_weight"] = baseline["analysis_weight_2022"]
    baseline["cluster"] = baseline["secu"]
    baseline["strata_value"] = baseline["stratum"]
    summary = {
        "eligible": int(eligible.sum()),
        "initial_improvement": int(np.nansum(baseline["initial_improvement"])),
        "sustained_improvement": int(np.nansum(baseline["sustained_improvement"])),
        "maintenance_eligible": int(baseline["maintenance_among_initial_improvers"].notna().sum()),
    }
    return intervals, baseline, summary


def build_chns() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    path = OUTPUT_ROOT / "chns_analysis"
    intervals = pd.read_csv(path / "chns_intervals.csv", dtype={"person_id": "string"}, low_memory=False)
    baseline = pd.read_csv(path / "chns_baseline_2006.csv", dtype={"person_id": "string"}, low_memory=False)
    bp_columns = [
        "person_id", "sbp_2006", "dbp_2006", "med_current_2006", "state5_2006",
        "sbp_2009", "dbp_2009", "med_current_2009", "state5_2009",
        "sbp_2011", "dbp_2011", "med_current_2011", "state5_2011",
    ]
    intervals = intervals.merge(baseline[bp_columns], on="person_id", how="left")
    intervals["start_sbp"] = np.where(intervals["start_year"].eq(2006), intervals["sbp_2006"], intervals["sbp_2009"])
    intervals["start_dbp"] = np.where(intervals["start_year"].eq(2006), intervals["dbp_2006"], intervals["dbp_2009"])
    intervals["end_sbp"] = np.where(intervals["end_year"].eq(2009), intervals["sbp_2009"], intervals["sbp_2011"])
    intervals["end_dbp"] = np.where(intervals["end_year"].eq(2009), intervals["dbp_2009"], intervals["dbp_2011"])
    intervals = add_core_and_buffered_endpoints(intervals)
    intervals["cohort"] = "CHNS"
    intervals["cluster"] = intervals["commid"]
    intervals["strata_value"] = np.nan
    intervals["interval_years"] = intervals["end_year"] - intervals["start_year"]

    numeric(baseline, ["state5_2006", "state5_2009", "state5_2011", "analysis_weight_2011"])
    eligible = baseline["state5_2006"].eq(3) & baseline[["state5_2009", "state5_2011"]].notna().all(axis=1)
    baseline["initial_improvement"] = np.where(
        eligible, baseline["state5_2009"].isin([1, 2]).astype(float), np.nan
    )
    baseline["sustained_improvement"] = np.where(
        eligible,
        (baseline["state5_2009"].isin([1, 2]) & baseline["state5_2011"].isin([1, 2])).astype(float),
        np.nan,
    )
    baseline["maintenance_among_initial_improvers"] = np.where(
        eligible & baseline["state5_2009"].isin([1, 2]), baseline["state5_2011"].isin([1, 2]).astype(float), np.nan
    )
    baseline["cohort"] = "CHNS"
    baseline["model_weight"] = baseline["analysis_weight_2011"]
    baseline["cluster"] = baseline["commid"]
    baseline["strata_value"] = np.nan
    summary = {
        "eligible": int(eligible.sum()),
        "initial_improvement": int(np.nansum(baseline["initial_improvement"])),
        "sustained_improvement": int(np.nansum(baseline["sustained_improvement"])),
        "maintenance_eligible": int(baseline["maintenance_among_initial_improvers"].notna().sum()),
    }
    return intervals, baseline, summary


def build_elsa() -> tuple[pd.DataFrame, dict[str, int]]:
    path = OUTPUT_ROOT / "elsa_analysis"
    intervals = pd.read_csv(path / "elsa_intervals.csv", dtype={"person_id": "string"}, low_memory=False)
    intervals["start_sbp"] = intervals["sbp_w6"]
    intervals["start_dbp"] = intervals["dbp_w6"]
    intervals["end_sbp"] = intervals["sbp_out"]
    intervals["end_dbp"] = intervals["dbp_out"]
    intervals = add_core_and_buffered_endpoints(intervals)
    intervals["cohort"] = "ELSA"
    intervals["cluster"] = intervals["idahhw6"]
    intervals["strata_value"] = np.nan
    return intervals, {"eligible": 0, "initial_improvement": 0, "sustained_improvement": 0, "maintenance_eligible": 0}


hrs_intervals, hrs_sustained, hrs_summary = build_hrs()
chns_intervals, chns_sustained, chns_summary = build_chns()
elsa_intervals, elsa_summary = build_elsa()

hrs_intervals.to_csv(OUT / "hrs_extended_intervals.csv", index=False, encoding="utf-8")
chns_intervals.to_csv(OUT / "chns_extended_intervals.csv", index=False, encoding="utf-8")
elsa_intervals.to_csv(OUT / "elsa_extended_intervals.csv", index=False, encoding="utf-8")
hrs_sustained.to_csv(OUT / "hrs_sustained_improvement.csv", index=False, encoding="utf-8")
chns_sustained.to_csv(OUT / "chns_sustained_improvement.csv", index=False, encoding="utf-8")

event_rows: list[dict[str, object]] = []
for cohort, frame in [("HRS", hrs_intervals), ("CHNS", chns_intervals), ("ELSA", elsa_intervals)]:
    for endpoint in [
        "hypertension_onset", "untreated_bp_improvement", "core3_hypertension_onset",
        "core3_hypertension_improvement", "buffer5_hypertension_onset", "buffer5_untreated_improvement",
        "buffer10_hypertension_onset", "buffer10_untreated_improvement",
    ]:
        eligible = frame[endpoint].notna()
        event_rows.append(
            {
                "cohort": cohort,
                "endpoint": endpoint,
                "eligible_intervals": int(eligible.sum()),
                "events": int(frame.loc[eligible, endpoint].sum()),
                "participants": int(frame.loc[eligible, "person_id"].nunique()),
            }
        )
pd.DataFrame(event_rows).to_csv(OUT / "dynamic_endpoint_event_gate.csv", index=False, encoding="utf-8-sig")

sustained_rows = []
for cohort, summary in [("HRS", hrs_summary), ("CHNS", chns_summary), ("ELSA", elsa_summary)]:
    sustained_rows.append({"cohort": cohort, **summary, "confirmatory_gate": "PASS" if summary["sustained_improvement"] >= 20 else "FAIL"})
pd.DataFrame(sustained_rows).to_csv(OUT / "sustained_improvement_event_gate.csv", index=False, encoding="utf-8-sig")

# Strict temporal-sequence feasibility gate for a hard-outcome extension.
hrs_long = pd.read_csv(
    OUTPUT_ROOT / "hrs_formal_cohort" / "hrs_o1_person_wave_2014_2022.csv",
    dtype={"person_id": "string"}, low_memory=False,
)
numeric(hrs_long, ["year", "state5", "died_by_wave"])
hrs18 = hrs_long.loc[hrs_long["year"].eq(2018), ["person_id", "state5"]].rename(columns={"state5": "state_2018"})
hrs22 = hrs_long.loc[hrs_long["year"].eq(2022), ["person_id", "died_by_wave"]].rename(columns={"died_by_wave": "death_2022"})
hrs_hard = hrs18.merge(hrs22, on="person_id", how="inner")
hrs_hard_eligible = hrs_hard["state_2018"].notna() & hrs_hard["death_2022"].notna()

numeric(chns_sustained, ["state5_2009", "state5_2011", "cvd_history_2009", "cvd_history_2011"])
chns_hard_eligible = (
    chns_sustained["state5_2009"].notna()
    & chns_sustained["cvd_history_2009"].eq(0)
    & chns_sustained["cvd_history_2011"].notna()
)
hard_gate = pd.DataFrame(
    [
        {
            "cohort": "HRS", "candidate_outcome": "Death during 2018-2022", "strict_three_windows": "YES",
            "eligible": int(hrs_hard_eligible.sum()), "events": int(hrs_hard.loc[hrs_hard_eligible, "death_2022"].eq(1).sum()),
            "harmonized_ready": "YES", "gate": "COHORT-SPECIFIC ONLY",
        },
        {
            "cohort": "CHNS", "candidate_outcome": "Incident MI/stroke history during 2009-2011", "strict_three_windows": "YES",
            "eligible": int(chns_hard_eligible.sum()),
            "events": int(chns_sustained.loc[chns_hard_eligible, "cvd_history_2011"].eq(1).sum()),
            "harmonized_ready": "LIMITED", "gate": "EVENT COUNT REVIEW",
        },
        {
            "cohort": "ELSA", "candidate_outcome": "Post-W8/W9 cardiovascular event or death", "strict_three_windows": "POTENTIALLY",
            "eligible": np.nan, "events": np.nan, "harmonized_ready": "NO", "gate": "REQUIRES NEW LATER-WAVE BUILD",
        },
    ]
)
hard_gate.to_csv(OUT / "hard_outcome_feasibility_gate.csv", index=False, encoding="utf-8-sig")

summary = {
    "interval_rows": {"HRS": len(hrs_intervals), "CHNS": len(chns_intervals), "ELSA": len(elsa_intervals)},
    "sustained_improvement": {"HRS": hrs_summary, "CHNS": chns_summary, "ELSA": elsa_summary},
    "hard_outcome_gate": hard_gate.replace({np.nan: None}).to_dict("records"),
}
(OUT / "dynamic_data_build_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False, indent=2))
