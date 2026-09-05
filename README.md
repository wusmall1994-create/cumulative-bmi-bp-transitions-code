# Cumulative excess BMI burden and blood pressure transitions

This repository contains the statistical analysis code for a harmonized longitudinal study of cumulative excess BMI burden and blood pressure progression and improvement in the Health and Retirement Study (HRS), China Health and Nutrition Survey (CHNS), and English Longitudinal Study of Ageing (ELSA).

Only code and documentation are included. No manuscript, raw cohort data, participant-level derived data, model outputs, credentials, local absolute paths, email addresses, telephone numbers, or other personal information are distributed here.

## Analysis structure

The numbered scripts reproduce the analysis in dependency order:

1. `code/hrs/`: construct the HRS cohort, prepare interval data, and fit transition-specific models.
2. `code/chns/`: construct the CHNS cohort and fit the harmonized models.
3. `code/elsa/`: construct the ELSA cohort and fit the harmonized models.
4. `code/meta/01`–`02`: combine cohort-specific estimates.
5. `code/meta/03`–`05`: construct the dynamic transition datasets, fit progression/improvement models, and apply small-sample meta-analytic inference.
6. `code/meta/06`: compare cumulative burden with baseline, last, mean, and time-weighted mean BMI.
7. `code/meta/07`: run the prespecified methodological sensitivity analyses, including baseline blood pressure adjustment, CHNS-specific BMI thresholds, percentile-standardized probabilities, continuous blood pressure change, weighting sensitivity, and robustness exclusions.

These are transition-specific discrete-time models; the code does not claim a continuous-time Markov model.

## Data are not included

HRS, CHNS, and ELSA data must be obtained directly from their custodians under the applicable terms of use. Do not place restricted data inside the Git repository. The `.gitignore` file blocks common participant-level and document formats, but users remain responsible for complying with each data-use agreement.

Set the following environment variables to the directories containing the authorized, extracted source files:

```text
HRS_DATA_DIR=/path/to/hrs
CHNS_DATA_DIR=/path/to/chns
ELSA_DATA_DIR=/path/to/elsa
BMI_BP_OUTPUT_DIR=/path/to/private/output   # optional
```

The HRS builder expects the RAND longitudinal file, the 2022 cross-wave tracker, the public region file, and the 2014/2018/2022 RAND Fat Files in the relative layout documented at the top of `code/hrs/01_build_cohort.py`. The CHNS builder expects `pexam_00.sas7bdat`, `surveys_pub_12.sas7bdat`, `educ_12.sas7bdat`, and `mast_pub_12.sas7bdat`. The ELSA builder expects the Wave 2, 4, 6, and 8/9 nurse/core files, the Gateway harmonized file, and the public urban/rural file named in `code/elsa/01_build_cohort.py`.

## Software

- Python 3.11 or later
- R 4.3 or later

Install Python and R dependencies:

```bash
python -m pip install -r requirements.txt
Rscript install_r_packages.R
```

## Reproduction

From the repository root, preview the complete execution plan:

```bash
python run_pipeline.py --dry-run
```

Run all stages after configuring authorized source data:

```bash
python run_pipeline.py
```

Resume at a named stage or run one stage only:

```bash
python run_pipeline.py --from-stage dynamic-build
python run_pipeline.py --only method-sensitivity
```

All participant-level intermediate files and results are written under `BMI_BP_OUTPUT_DIR`, or under the ignored `outputs/` directory by default.

## Privacy and reproducibility safeguards

- Input and output data are never tracked.
- Repository history is scanned for restricted file types, local paths, credentials, and direct contact details before release.
- Random seeds are fixed where resampling is used.
- Cohort-specific models are fitted independently before meta-analysis.
- The code preserves treatment-aware five-state outcomes and harmonized three-state sensitivity definitions.

## CHNS acknowledgment

This research uses data from China Health and Nutrition Survey (CHNS). We thank the National Institute of Nutrition and Food Safety, China Center for Disease Control and Prevention, Carolina Population Center (5 R24 HD050924), the University of North Carolina at Chapel Hill, the NIH (R01-HD30880, DK056350, R24 HD050924, and R01-HD38700) and the Fogarty International Center, NIH for financial support for the CHNS data collection and analysis files from 1989 to 2011 and future surveys, and the China-Japan Friendship Hospital, Ministry of Health for support for CHNS 2009.

## Citation and license

Please cite the archived software release using `CITATION.cff`. The analysis code is released under the MIT License. The license applies only to this repository's code and does not apply to HRS, CHNS, ELSA, or any derived participant-level data.

