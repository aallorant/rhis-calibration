# Statistical calibration of routine health data with survey audits

Code and results for the paper *A statistical calibration framework for subnational
health surveillance from routine data and survey audits*. The framework takes a routine
health register as the estimator of a subnational quantity (here undiagnosed HIV
positivity by district, age and sex) and uses an intermittent probability survey to
estimate and remove the register's selection bias. The paper applies it to Mozambique,
Zimbabwe, Malawi and Lesotho, and studies by simulation the conditions under which the
calibrated estimate improves on current practice.

## Layout

| Folder | Content |
|---|---|
| `R/` | The analysis pipeline, run from the repository root with `Rscript R/<script>.R` |
| `tmb/` | TMB templates for the province-anchored and covariate models |
| `simulation/` | The simulation study (`simulation/R/`) and its outputs (`simulation/output/`) |
| `data/` | Routine HIV testing counts aggregated by country, year, age group and sex |
| `output/` | Model outputs behind every table and figure of the paper, where the scripts read them |
| `figures/` | The four figures of the paper |

## Pipeline

1. `R/01_standardise_rhis.R`, `R/02_standardise_phia.R`, `R/02_standardise_phia_moz.R`:
   standardise the routine registers and the PHIA survey microdata into district and
   province cells (register logit positivity `z` with variance `v`; survey logit
   undiagnosed positivity with design variance `psi`).
2. `R/02d_standardise_phia_round0.R`, `R/02e_standardise_imasida2015.R`: the earlier
   survey round used in the two-round comparison.
3. `R/03a_moz_covariates.R`, `R/03b_country_selftest_covariates.R`,
   `R/07_load_pepfar_covariates.R`: district covariates for the covariate arm of the fit.
4. `R/03_calibration_model_moz.R`, `R/build_moz_district_survey_design.R`: district-level
   survey estimates for Mozambique, used by the transport diagnostic.
5. `R/04_two_stage_calibration.R`, `R/build_conditional_anchors.R`,
   `R/build_conditional_stratum_anchors.R`: per-country inputs of the main fit.
6. `R/RUN_canonical.R`: the main fit (differenced Fay-Herriot model with the block
   covariance, by province and age-sex stratum), with `R/mainfit_agesex.R`,
   `R/refit_block_eblup_canonical.R`, `R/mainfit_national_popweighted.R`,
   `R/mainfit_inla_crosscheck.R`, `R/collapsed_admin2_inla_spatial.R` and
   `R/eblup_vs_latentfield_agreement.R` (with `R/fig_eblup_vs_tmb_scatter.R`) for the appendix comparisons.
7. `R/transport_diagnostic.R`, `R/build_table_transport.R`,
   `R/build_transport_satterthwaite.R`, `R/sigma_w_survey_only.R`: the diagnostic of the
   transport assumption.
8. `R/04d_temporal_lambda.R`, `R/04f_volatility_ratio.R`, `R/mechanical_coupling_null.R`:
   the two-round comparison of the offset.
9. `R/regen_table_districts.R`, `R/build_table_predvar_h.R`, `R/build_table_weights.R`,
   `R/build_retest_bound.R`, `R/national_undiagnosed_positivity.R`: the remaining tables.
10. `R/fig_mainfit_coefficients.R`, `R/fig_transport_provinces.R`: Figures 1 and 2.

Simulation: `simulation/R/00_config.R` to `03_designs.R` define the generative model,
the estimators and the audit designs; `04_run.R`, `07_resolution.R`, `08_coverage.R`,
`08b_coverage_oracle.R`, `16_sigma_u_sweep.R`, `17_sensitivity.R` and
`18_tau_share_sweep.R` produce the results; `13_scenario_paths.R` and `14_epi_realism.R`
draw the figures. Every simulation script is run from the repository root.

## Data

Steps 1 to 5 read the survey microdata and the district-level routine registers, which
are not included:

- PHIA survey microdata are available to registered users from the PHIA Project data
  portal (https://phia-data.icap.columbia.edu/).
- The routine HIV testing data are the property of the national ministries of health.
  `data/routine_hiv_testing_country_year_age_sex.csv` gives the counts of tests and
  positive tests aggregated by country, year, age group and sex.
- The PEPFAR Monitoring, Evaluation and Reporting indicators used in the two-round
  comparison were obtained from the amfAR PEPFAR data portal (https://mer.amfar.org/).

`output/` holds the outputs of steps 6 to 9 at the resolution reported in the paper
(country, province and age-sex stratum), at the paths where the scripts read them, and
`simulation/output/` holds every simulation result. The two figure scripts of step 10 run
from these files alone and reproduce Figures 1 and 2; the tables can be checked against the
files in `output/canonical/estimates/`. The district-level registers and survey microdata
are needed only to rerun steps 1 to 9.

## Software

R (4.4) with TMB, INLA, dplyr, tidyr, data.table, ggplot2, cowplot, survey, sf, spdep,
readxl, readr, stringr and lubridate.

## Licence

MIT.
