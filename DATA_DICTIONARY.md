# Data dictionary

## derived_data/cohort_reports_v2_flags.csv — primary-suspect cohort (4499 reports)
| Variable | Description |
|---|---|
| primaryid, caseid | FAERS/AEMS report and case identifiers (latest version after deduplication) |
| target_drug | sotorasib or adagrasib (primary suspect) |
| drug_seq | Drug sequence number of the KRAS G12C inhibitor in the report |
| ici | Co-reported ICI generic names (semicolon-separated); empty if none |
| ici_exposed | TRUE if any ICI was listed in any role ("ICI co-reported") |
| hepato | TRUE if any of the 36 hepatic PTs (hepato_pts.txt) was reported |
| severe_hep | TRUE if a liver failure PT was reported |
| outc_DE/LT/HO/DS, serious | Death, life-threatening, hospitalization, disability; any of these |
| sex | F, M or empty (not reported) |
| country | Country of occurrence (ISO 2-letter); empty = not reported |
| event_dt | Event date (YYYYMMDD, may be partial) |
| rept_cod | Report type: EXP expedited, PER periodic, DIR direct |
| age_yr | Age in years (converted from age_cod); empty = not reported |
| indication | Indication PTs linked to the KRAS G12C inhibitor (semicolon-separated) |
| tumor | NSCLC, CRC, Pancreatic, Other, Unknown (text-matched from indication). In all analyses, "Other" reports whose indication contains "product used for unknown indication", "K-ras gene mutation", "neoplasm malignant" or "illness" are recoded to "not specified" (Unknown) |
| tto_days, tto | Days from earliest KRAS G12C inhibitor start date to event date (0-3650 d) |
| ici_to_kras_days | Days from latest recorded ICI date to KRAS G12C inhibitor start (superseded by ici_interval_vs_event.csv) |
| occp | Reporter type: MD, PH, OT, HP (health professional), LW, CN (consumer); empty = not reported |
| quarter, year | Data release quarter and year (year of receipt) |
| ici_roles | ICI role codes (PS, SS, I, C) |
| ici_suspect, ici_concomitant_only | ICI recorded as suspect/interacting; as concomitant only |
| trial_coded | KRAS G12C inhibitor recorded under the development code (AMG 510) |
| rpsr, study_report | Report source code; flagged as study report |
| switch | Adagrasib report that also lists sotorasib |
| ivg | Interval group (original grouping; superseded) |

## derived_data/cohort_any_role.csv — any-role cohort (Table 9)
| Variable | Description |
|---|---|
| target_drug, drug_seq | KRAS G12C inhibitor in its highest-ranked role, and its sequence number |
| rk, role | Role rank (0 PS, 1 SS, 2 I, 3 C) and role code |
| other variables | As above |

## derived_data/ici_interval_vs_event.csv — ICI co-reported reports
Classification of ICI start/stop dates and KRAS G12C inhibitor start date relative to the event date (documented overlap, sequential use, ICI start only, ICI after event, stop without start, insufficient dates).

## Model covariates (all adjusted models)
Country (US vs other/not reported), report type (expedited ref.; periodic; direct), reporter type (health professional ref.; consumer; not reported), year of receipt (2021-2022 ref.; 2023; 2024; 2025; 2026), indication (NSCLC ref.; colorectal; other; not specified; omitted in NSCLC-only models), age (<65 y ref.; >=65 y; not reported), sex (male ref.; female; not reported).
