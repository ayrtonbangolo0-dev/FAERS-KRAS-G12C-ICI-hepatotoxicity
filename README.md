# Co-reported immune checkpoint inhibitors and hepatotoxicity with KRAS G12C inhibitors (FAERS/AEMS)

Analysis code for: Bangolo A, *et al.* Co-reported immune checkpoint inhibitors and hepatotoxicity with KRAS G12C inhibitors: A pharmacovigilance study of sotorasib and adagrasib. *World J Gastrointest Pharmacol Ther* (under review).

## Data
FDA Adverse Event Reporting System (FAERS) / Adverse Event Monitoring System (AEMS) quarterly ASCII files, 2021Q2-2026Q3 (22 files, `faers_ascii_YYYYqN.zip`), freely available from
https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html. Place the zip files in one folder; the scripts expect `~/mnt/Downloads` (edit the `DL`/`B` paths at the top of each script) and write intermediate files to `~/fw` and `work/`.

## Source-file checksums
`source_file_provenance.txt` lists each quarterly file with its size, download time and SHA-256 checksum (Supplementary Table 4); verify with `sha256sum faers_ascii_*.zip`.

## Requirements
Python 3.10, pandas 2.3.3, NumPy 2.2.6, matplotlib (figure only), GNU coreutils (sort, join, awk). No other packages: logistic and Firth models are implemented in `11_round2.py` (`fit`), and validated in `14_reviewer2.py` (section D).

## Order of execution
| Step | Script | Purpose |
|---|---|---|
| 1 | `01_build.py <quarter>` (`run_v2.sh` loops over quarters) | Extract reports with sotorasib/adagrasib (any role) and co-reported ICIs; drug, demo, reac, outc, indi, ther tables |
| 2 | `02_dedup.sh <workdir>` | FDA deduplication (latest FDA_DT per CASEID, highest PRIMARYID on ties), deleted cases removed; database-wide PT counts |
| 3 | `03_analysis.py` | Primary-suspect cohort, ROR/PRR/IC against the database, within-cohort comparisons, time to onset (Tables 2, 3, 5) |
| 4 | `04_tables_figure.py` | Table 1 and the original forest plot (now Figure 2; the revised plot was drawn separately) |
| 5 | `05_sensitivity.py`, `09_table6_grid.py` | Sensitivity analyses (Table 6) |
| 6 | `06_hepatic_breakdown.py` | Hepatotoxicity categories and PTs (Table 4, Supplementary Table 1) |
| 7 | `07_rescan.py <quarter>` (`run_rescan.sh`) | Database-wide ICI, docetaxel and development-code flags |
| 8 | `08_reviewer_analyses.py` | Interaction models, Omega (any role) |
| 9 | `10_doce_scan.py <quarter>` | Docetaxel primary-suspect covariates |
| 10 | `11_round2.py` | Docetaxel reference models, pooled NSCLC models, primary-suspect Omega, Firth models (Tables 7, 8) |
| 11 | `12_round3.py` | Immune-mediated hepatitis, adagrasib case review |
| 12 | `13_reviewer1.py` | Treatment dates against the event date (Supplementary Tables 2, 3) |
| 13 | `10b_doce_anyrole.py <quarter>` | Docetaxel in any role (for Table 9) |
| 14 | `14_reviewer2.py` | Broader cohorts, complete cases, K-ras indication coding, Firth validation (Table 9) |
| 15 | `15_editorial.py`, `15b_separation.py` | Study flow counts (Figure 1), model diagnostics, events per parameter and separation checks (Supplementary Table 5) |
| 16 | `fig1_forest_plot.py`, `fig2_study_flow.py` | Revised forest plot (Figure 2) and study flow diagram (Figure 1), drawn with matplotlib from the estimates in the tables |

`hepato_pts.txt` lists the 36 hepatic MedDRA preferred terms (Supplementary Table 1).

## Derived data
The CSV files (`cohort_reports_v2_flags.csv`, `cohort_any_role.csv`, `ici_interval_vs_event.csv`) contain the analysis cohorts (de-identified public FAERS/AEMS report identifiers and derived variables; see `DATA_DICTIONARY.md`) and the JSON files are the outputs of the analysis scripts. All files are in the top-level folder.

## Licence
Code: MIT. FAERS/AEMS data are public-domain US government data.
