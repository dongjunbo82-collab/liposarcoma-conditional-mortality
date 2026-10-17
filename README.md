# Statistical analysis code for conditional mortality in liposarcoma

Version 1.0.0. This repository contains statistical analysis source code and its
technical documentation. It contains no SEER records, derived patient-level
files, mortality-rate files, manuscript documents, figures, credentials, author
contact information or institutional addresses. The creator label is collective
to respect the requested privacy boundary.

## What the code does

- Defines an adult liposarcoma cohort from an authorized 23-column SEER export.
- Estimates conditional cumulative incidence with competing events at survival
  landmarks, using shared patient bootstrap resamples.
- Describes histology, anatomical site, attained age, grade and recorded treatment.
- Performs diagnosis-period and other supporting sensitivity analyses.
- Fits supporting standardized and continuous-age models.
- Compares observed cardiovascular mortality with public population rates,
  including exact Poisson intervals and date/age assumption checks.
- Reconstructs event/censoring support and generates statistical figures locally.

The analyses are exploratory. Population-level estimates do not represent
recurrence risk, independent treatment effects or individualized predictions.
Intervals are pointwise; small groups can have limited support.

## Requirements

Python 3.12 or later with packages in `requirements.txt`; R with data.table,
survival, splines, ggplot2 and svglite. The development environment used R 4.6.1
and SEER*Stat 9.0.43.0 for export. Python scripts invoke `Rscript` on PATH;
set `RSCRIPT` to another executable path if necessary. Dependency checks do not
install software. No R or Python installation is bundled here.

## Inputs and local-only outputs

Obtain SEER access through the National Cancer Institute and export the fields
in `docs/input_schema.md` in the specified order, as tab-delimited labels with
a header. Place the export at `inputs/seer_export.tsv` or set `SEER_EXPORT` to
your authorized local file. This code does not confer data access or authorize
redistribution. Keep SEER data under the applicable data-use agreement.

The output directory defaults to `outputs/population_study`; set
`LPS_OUTPUT_DIR` to an absolute local directory to keep generated data elsewhere.
Some generated files contain patient-level information and must remain local.
All input and output directories are git-ignored. Never force-add those files.

The death-certificate decomposition needs the official SEER COD dictionary.
The SMR analysis also needs CDC WONDER reference rates. These inputs are not
bundled. Download public references by running:

```sh
python analysis/download_reference_inputs.py
python analysis/validate_wonder_reference.py
```

For only the COD dictionary, use `--seer-only`. The downloader sends fixed public
aggregate queries only and does not read or send SEER records. CDC endpoints may
change or return rate limits; the downloader spaces requests by 16 seconds and
stops on errors rather than inventing replacement rates. See
`docs/reference_inputs.md` for the expected queries and manual alternative.

## Run order

From the repository root, after installing dependencies and preparing inputs:

```sh
Rscript tests/test_estimators.R
python run_analysis.py --stage primary
python run_analysis.py --stage secondary
python run_analysis.py --stage models
python run_analysis.py --stage smr
python run_analysis.py --stage figures
```

Alternatively, `python run_analysis.py --stage all` runs these stages in order.
The primary and supporting nonparametric analyses use 1,000 bootstrap samples;
standardization uses 500. Full models may take appreciable time. Random seeds
are fixed in each script. Scripts intentionally retain assertions for the
original cohort and event coding so a changed export fails visibly instead of
silently producing apparently comparable results.

`population_standardized.R` now reads the canonical primary cohort directly,
rather than requiring an older redundant preprocessing export. This is a
packaging change; variable mappings and statistical model specifications are
retained. File locations and executable paths have also been parameterized.

## Code map and verification

See `docs/analysis_map.md`. The small artificial test checks analytic competing
risk probabilities, a separate survival::survfit calculation, positive-landmark
eligibility and event/censor ties. It contains no patient observations.
`analysis/verify_primary_estimates.py` independently checks all primary point
estimates after a local run. Further verification is embedded in the R scripts.
This computational verification is not independent clinical validation.

## License and citation

Source code is released under the MIT License. Data providers retain all rights
and restrictions over their own data; the software license does not cover them.
Use `CITATION.cff` for a code citation. A DOI will be added only after the release
has actually been archived and its DOI resolves.
