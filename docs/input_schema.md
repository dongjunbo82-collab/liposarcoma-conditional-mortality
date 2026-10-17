# Required SEER export schema

Input: tab-separated labels, UTF-8 (UTF-8 BOM accepted), one header row, exactly
23 fields in the following order. The code uses positional fields v1 through
v23; header wording is not used to assign fields. No example data are included.

| Position | Field |
|---|---|
| 1 | Patient ID |
| 2 | Year of diagnosis |
| 3 | Age recode with single ages and 90+ |
| 4 | Sex |
| 5 | Race recode (White, Black, Other) |
| 6 | Primary Site - labeled |
| 7 | Histologic Type ICD-O-3 |
| 8 | Behavior code ICD-O-3 |
| 9 | Diagnostic Confirmation |
| 10 | Sequence number |
| 11 | Grade Recode (thru 2017) |
| 12 | Derived Summary Grade 2018 (2018+) |
| 13 | Summary stage 2000 (1998-2017) |
| 14 | Combined Summary Stage with Expanded Regional Codes (2004+) |
| 15 | RX Summ--Surg Prim Site (1998-2022) |
| 16 | RX Summ--Surg Prim Site 2023 (2023+) |
| 17 | Radiation recode |
| 18 | Chemotherapy recode (yes, no/unk) |
| 19 | Survival months |
| 20 | Vital status recode (study cutoff used) |
| 21 | SEER cause-specific death classification |
| 22 | SEER other cause of death classification |
| 23 | COD to site recode |

Use the SEER Research Data 17 Registries November 2025 submission with diagnoses
through 2023. The original export selects malignant behavior, known age and
histology codes 8850–8858 across available diagnosis years 2000–2023. Scripts
apply the adult, histologic-confirmation, 2004+ and first/only-primary rules.
Unknown stage and recorded nonsurgical-treatment categories are retained as
specified. SEER attribution and certificate-coded causes are separate endpoints.
Zero-month survival is assigned 0.5 month in the primary convention; at positive
landmarks recorded survival must exceed the landmark. A midpoint sensitivity
adds 0.5 month to all recorded survival durations. All times start at diagnosis.

Before running, verify the field order against your own dictionary. Numeric
exports or altered labels need explicit recoding and should not be substituted
without checking the parsing and event definitions.
