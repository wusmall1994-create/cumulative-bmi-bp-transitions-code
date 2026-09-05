"""Run the code-only HRS, CHNS, ELSA, and cross-cohort analysis pipeline."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent

STAGES = [
    ("hrs-build", [sys.executable, "code/hrs/01_build_cohort.py"]),
    ("hrs-prepare", [sys.executable, "code/hrs/02_prepare_analysis.py"]),
    ("hrs-model", ["Rscript", "code/hrs/03_fit_transition_models.R"]),
    ("chns-build", [sys.executable, "code/chns/01_build_cohort.py"]),
    ("chns-model", ["Rscript", "code/chns/02_fit_transition_models.R"]),
    ("elsa-build", [sys.executable, "code/elsa/01_build_cohort.py"]),
    ("elsa-model", ["Rscript", "code/elsa/02_fit_transition_models.R"]),
    ("meta-two", ["Rscript", "code/meta/01_harmonized_meta_hrs_chns.R"]),
    ("meta-three", ["Rscript", "code/meta/02_harmonized_meta_three_cohort.R"]),
    ("dynamic-build", [sys.executable, "code/meta/03_build_dynamic_data.py"]),
    ("dynamic-model", ["Rscript", "code/meta/04_dynamic_transition_models.R"]),
    ("robust-meta", ["Rscript", "code/meta/05_robust_meta_analysis.R"]),
    ("incremental-value", ["Rscript", "code/meta/06_incremental_value_analysis.R"]),
    ("method-sensitivity", ["Rscript", "code/meta/07_methodological_sensitivity.R"]),
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from-stage", choices=[name for name, _ in STAGES])
    parser.add_argument("--only", choices=[name for name, _ in STAGES])
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    stages = STAGES
    if args.from_stage:
        start = [name for name, _ in stages].index(args.from_stage)
        stages = stages[start:]
    if args.only:
        stages = [item for item in stages if item[0] == args.only]

    env = os.environ.copy()
    env["BMI_BP_PROJECT_ROOT"] = str(ROOT)
    for name, command in stages:
        print(f"[{name}] {' '.join(command)}", flush=True)
        if not args.dry_run:
            subprocess.run(command, cwd=ROOT, env=env, check=True)


if __name__ == "__main__":
    main()

