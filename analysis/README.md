# Analysis compendium

The scripts and result files behind the tables and figures of
Frutos-Galarza et al. (2026). This directory is not part of the R package:
it is excluded from the build through `.Rbuildignore` and exists so that the
published numbers can be traced back to the code that produced them.

The package itself, in `R/`, implements the method. These scripts use it to
run the Monte Carlo campaigns and the Tennessee Eastman application.

---

## Running the scripts

The scripts were written to run from a working directory containing two
sibling folders, and they have been left exactly as they were when they
produced the published figures. To reproduce a result:

1. Create a working directory anywhere, with two subfolders:

2. Set that directory as the working directory in R:

```r
setwd("path/to/working_directory")
```

3. Run any script by sourcing it from wherever this repository sits:

```r
source("path/to/robustT2AFM/analysis/01_scripts/SCRIPT_FINAL_SIMULATIONS.R")
```

The output lands in `02_resultados/`. Compare it against the corresponding
file in `02_results/` of this repository, which holds the versions the paper
was written from.

Scripts that only run simulations need no data and work with an empty
`04_datos/`. Scripts that read the Tennessee Eastman collection need the
files described next.

## The Tennessee Eastman data

The data are not redistributed here. They belong to the collection published
by Rieth, Amsel, Tran and Cook, available at

  https://doi.org/10.7910/DVN/6C3JR1

Download the four CSV files and place them in `04_datos/` under these names:

The fault-free training set characterises normal operation and is used only
for variable selection. The testing sets supply the runs that make up
Phase 1 and Phase 2. Keeping the two apart is deliberate: it prevents the
selection of variables from being made on the same runs the method is then
evaluated on. Computing the Table 6 summary on the testing set instead
shifts the coefficient of variation of `xmv_9` from 0.056643 to 0.058679.

## Reproducibility

The numeric anchors, fixed parameters and seeds are documented in
`REPRODUCIBILITY.md` at the repository root. Every campaign derives its
random streams from a single base, `SEED_BASE = 2026`, so a rerun on the
same data must return the same numbers.

`21_phase2_stability_resampling.R` takes about 93 minutes on eight parallel
workers. The other Monte Carlo campaigns, which run 2000 replicates each,
take on the order of tens of minutes. The Tennessee Eastman scripts are
dominated by the time spent reading the data files, which are large.

---

## Where each number in the paper comes from

| Script | Writes | Appears in |
|---|---|---|
| `ZTEST_TABLE1_FAITHFUL_v2.R` | `table1_ztest_faithful.csv`, `far_per_replicate.csv` | Table 1, first two rows |
| `KEFF_EFFECT_v2.R` | `table1_row_keff.csv` | Table 1, K_eff row |
| `CENTERED_BOOTSTRAP_v5.R` | `table1_rows_bootstrap_v5.csv` | Table 1, three bootstrap rows |
| `16_ci_paired_far_difference.R` | `ci_paired_far.csv` | Section 3.2, paired confidence interval |
| `15_arl0_convergence_K_v3.R` | `arl0_convergence_K_v3.csv` | Table 2 |
| `20_arl0_clean_vs_contaminated_v2.R` | `arl0_clean_vs_contaminated_v2.csv` | Section 3.2, clean against contaminated Phase 1 |
| `SCRIPT_FINAL_SIMULATIONS.R` | `table_3_3_power.csv`, `table_3_4_sensitivity.csv`, `table_3_4b_contamination.csv` | Tables 3 and 5 |
| `SCRIPT_FINAL_ABLATIONS_2x2.R` | `table_ablation_2x2_final.csv` | Table 4 |
| `22_rho_sensitivity_parallel.R` | `table_3_4c_rho.csv` | Table 5, correlation block |
| `Script_figures_detection_power.R` | reads `table_3_3_power.csv` | Figure 2 |
| `TEP_VARIABLES_SUMMARY.R` | `tep_variables_summary.csv` | Table 6 |
| `19_acf_with_step1.R` | `step_sensitivity_ljungbox_with_step1.csv` | Table 7, Table 8, Section 4.3 |
| `SUBSAMPLING_STEP_SENSITIVITY.R` | `step_sensitivity_detection.csv` | Table 8, false alarm column |
| `SCRIPT_Figure_ACF.R` | figure only | Figure 4 |
| `PHASE1_PIPELINE_STEP30.R` | no file; base pipeline | Section 4.4, and the basis of Figures 5 and 6 |
| `Script_Weights.R` | figure only | Figure 5 |
| `BRIDGE_Script_Control_Chart.R`, `Script_Control_Chart.R` | figure only | Figure 6 |
| `18_phase1_idv1_cross_matrix.R` | `phase1_cross_matrix.csv` | Table 9, Sections 4.4 and 4.6 |
| `TEP_ROBUST_CENTER.R` | `tep_robust_center_by_variable.csv`, `tep_robust_center_summary.csv` | Section 4.4 |
| `21_phase2_stability_resampling.R` | `phase2_stability_summary.csv`, `phase2_stability_by_composition.csv` | Section 4.5, twenty compositions |
| `TEP_EXTENDED_VALIDATION_STEP30.R` | `table_multiple_faults_step30.csv` | Table 10 |
| `23_REFERENCE_CENTER_AUDIT_v4.R` | `weighted_center_audit_v4.csv` | Section 4.6 |
| `24_consolidate_center_audit.R` | `center_audit_summary.csv`, `center_audit_full_detail.csv` | Section 4.6 |
| `TEP_BOOTSTRAP_LIMIT.R` | `tep_reproducibility.csv`, `tep_limits_comparison.csv` | Section 4.4, the 16.14 determinant ratio |

## Two notes on the result files

`tep_limits_comparison.csv` does not appear as a table. It records the
bootstrap limits computed on the Tennessee Eastman data, which rise to 54.98
against the analytic 19.69, and is the empirical basis for the statement in
Section 2.5 that the resampling alternatives were evaluated and discarded.

`center_audit_summary.csv` and `weighted_center_audit_v4.csv` hold the same
results. The first is rounded to four decimals by the consolidation script;
the second keeps full precision and is the one to cite.

`SCRIPT_FINAL_SIMULATIONS.R` also writes `table_3_2_calibration.csv`, which
is an intermediate output that no table of the paper uses. It is not shipped
here.
