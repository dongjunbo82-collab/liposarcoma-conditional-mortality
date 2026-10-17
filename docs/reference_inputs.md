# Public reference inputs

The repository supplies query code, not reference datasets. Local files expected:

- `survivorship/seer_cod_recode.txt` under `LPS_OUTPUT_DIR`: official SEER COD
  recode definition dated 03/01/2018, downloaded from
  https://seer.cancer.gov/codrecode/1969_d03012018/index.txt .
- `revision/wonder_20261016/`: CDC WONDER XML responses and the CSV created by
  `validate_wonder_reference.py`. That validator checks completeness, uniqueness,
  positive denominators and that deaths/population reproduce reported rates.

Queries use US residents, all races, separately by calendar year and sex,
underlying cardiovascular causes ICD-10 I00–I13, I20–I51 and I60–I78.
Use D76 for 2009–2017 and D158 for 2018–2023. Query five-year ages 20–24 through
80–84, and separately obtain 85+ from the ten-year-age option. Rates are crude
stratum-specific rates, not age-adjusted rates. Expected reference cells:
15 calendar years × 2 sexes × 14 age categories.

If the API is unavailable, run the equivalent queries in the official CDC
WONDER interface and retain XML responses matching the expected format. Do not
replace suppressed or unavailable rates with zero. Check the service’s current
terms and query interface. No SEER person-level information is needed by CDC.
Population rate versions were not harmonized to a common census vintage.

Official documentation: https://wonder.cdc.gov/wonder/help/ucd.html and
https://wonder.cdc.gov/wonder/help/ucd-expanded.html .
