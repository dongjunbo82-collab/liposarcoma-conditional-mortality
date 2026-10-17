# Analysis map

| Script | Purpose |
|---|---|
| population_primary.R | Cohort, competing-risk conditional estimates, shared bootstrap and support flags |
| verify_primary_estimates.py | Separate point-estimate reconstruction |
| followup_support.py | Events, censoring and follow-up-window support |
| survivorship_cod_map.py | COD label mapping using official dictionary |
| survivorship_extension.R | Cause decomposition, age bands and age standardization |
| revision_analyses.R | Treatment subsets and sensitivity intervals |
| population_standardized.R | Supporting covariate-standardized cumulative incidence |
| age_definition_sensitivity.R | Continuous attained-age competing-hazard models |
| age_model_checks.R | Linear-age and top-coded-age checks; apparent aggregate calibration |
| clinical_context_extension.R | Grade, treatment descriptions and myxoid coding checks |
| cardiovascular_smr.py | Expected person-time, SMRs and age/calendar sensitivity analyses |
| validate_cardiovascular_smr.R | Independent expected-count and Poisson-interval checks |
| cancer_period_sensitivity.R | All joint histology-site diagnosis-period checks with paired resampling |
| focus_figures.R | Final conditional mortality and age/cause figures |
| population_figures.R; survivorship_figures.R | Additional descriptive plots; run optionally before final figures |
| download_reference_inputs.py; validate_wonder_reference.py | Public reference retrieval and validation |

Inputs and outputs remain local. The code does not create or upload manuscripts.
The date/age allocation and general-population comparisons are approximations
specified in code; there is no external clinical validation or causal treatment
comparison. Bootstrap draws are shared within each analysis family. Independent
families use their stated seeds and can yield slightly different interval
endpoints for the same point estimate.
