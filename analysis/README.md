# Analysis compendium

The scripts and result files behind the tables and figures of
Frutos-Galarza et al. (2026), in its revised version. This directory is not
part of the R package: it is excluded from the build through `.Rbuildignore`
and exists so that the published numbers can be traced back to the code that
produced them.

The package itself, in `R/`, implements the method. These scripts use it to
run the Monte Carlo campaigns and the Tennessee Eastman application. The
compendium that accompanied the original manuscript remains available on
Zenodo at https://doi.org/10.5281/zenodo.22544989.

---

## Running the scripts

The scripts have been left exactly as they were when they produced the
published numbers, and they read and write folders through relative paths.
To reproduce a result:

1. Create a working directory anywhere, with these subfolders:

```
working_directory/
  01_scripts/                      copy PHASE1_PIPELINE_STEP30.R here
  02_resultados/                   copy the five version 0.2.2 files here (see below)
  03_figuras/
  04_datos/                        the Tennessee Eastman data
  05_revision_R1/02_resultados/
  05_revision_R1/03_T2_crudos/
```

2. Set that directory as the working directory in R:

```r
setwd("path/to/working_directory")
```

3. Run any script by sourcing it from wherever this repository sits:

```r
source("path/to/robustT2AFM/analysis/01_scripts/SCRIPT_M2_power.R")
```

The output lands in `05_revision_R1/02_resultados/`, and the figures in
`03_figuras/`. Compare it against the file with the same name in
`02_results/` of this repository, which holds the versions the paper was
written from. When a script needs a file written by another one, either run
that other script first or copy the file from `02_results/`.

The scripts need R 4.5 with robustT2AFM 0.3.0 and the packages MASS,
robustbase, rrcov, sn, ggplot2 and patchwork. Scripts that only run
simulations need no data and work with an empty `04_datos/`.

## The Tennessee Eastman data

The data are not redistributed here. They belong to the collection published
by Rieth, Amsel, Tran and Cook, available at

  https://doi.org/10.7910/DVN/6C3JR1

Download the files and place them in `04_datos/` under these names:

```
TEP_FaultFree_Training.csv
TEP_FaultFree_Testing.csv
TEP_Faulty_Testing.csv
```

The fault-free training set characterises normal operation and is used only
for variable selection (Table 15). The testing sets supply the runs that make
up Phase 1 and Phase 2. Keeping the two apart is deliberate: it prevents the
selection of variables from being made on the same runs the method is then
evaluated on.

## The T² values of each replicate

The long simulation scripts (M1, M2, M2b, M3 and M5) save the T² values of
every replicate in `05_revision_R1/03_T2_crudos/`, and the short scripts
(M1b, M1c, M2c, M3b, M5b, M2M3_diagnostics and V_ucl_simulated_vs_M5) read
them to compute the tables. They take 112 MB, so they are not kept in this
repository but in a separate Zenodo record: [DOI pending]. Unzip them into
that folder. Without them, the long scripts have to be run first, which takes
several hours.

## Reproducibility

Every script fixes its seeds at the top, so a rerun on the same data returns
the same numbers. In addition, 26 of the 31 scripts check numeric anchors,
values already known that must come out unchanged, and stop if any of them is
not reproduced. The exact versions of R and of each package used for every
run are recorded in `02_results/sessionInfo_R3_02_v030.txt`.

---

## Where each number in the paper comes from

| Script | Writes | Appears in |
|---|---|---|
| `SCRIPT_M1_limit.R` | `table_M1_T1_calibration.csv`, `table_M1_T2_K.csv`, `table_M1_T3_seeds.csv`, `table_M1_T3_pooled.csv`, `table_M1_T13_center.csv`, `table_M1_mstar.csv` | Section 2.2 (m* = 14; 16.7 and 15.7 after reweighting); Section 3.2 (1019 and 1326) |
| `SCRIPT_M1b_limit_posthoc.R` | `table_M1b_mstar_arl0.csv`, `table_M1b_empirical_limit.csv` | Sections 2.5 and 3.2 (m* = 13, 16.7 and 20) |
| `SCRIPT_M1c_recount_mstar14.R` | `table_M1c_cells_mstar.csv`, `table_M1c_tests.csv`, `table_M1c_empirical_limit.csv` | Tables 2, 3 and 4; Section 3.2 (tests, confidence interval, correction factor) |
| `SCRIPT_M5_simulated_limit.R` | `table_M5_simulated_limit.csv`, `table_M5_limit_spread.csv`, `table_M5_TEP.csv` | Table 5 (simulated and classical limits); Section 3.2 |
| `SCRIPT_M5b_recount_mstar14.R` | `table_M5b_eq8_mstar14.csv` | Table 5, Equation (8) column |
| `SCRIPT_R1_10_shift_directions.R` | `table_R1_10_shift_directions.csv` | Shift directions used by M2 (Section 3.1) |
| `SCRIPT_M2_power.R` | `table_M2_all.csv`, `table_M2_T4_power.csv`, `table_M2_T5_factorial.csv`, `table_M2_ncal.csv`, `table_M2_competitors.csv`, `table_M2_directions.csv`, `table_M2_paired.csv` | Tables 6, 7 and 13; Section 3.3 (50,000 batches) |
| `SCRIPT_M2M3_diagnostics.R` | `table_M3_W_diagnostic.csv` | Sections 3.3 and 3.7 (shift of the center, first eigenvalue) |
| `SCRIPT_M2b_sensitivity_sameseeds.R` | `table_M2b_all.csv`, `table_M2b_paired.csv`, `table_M2b_T6_sensitivity.csv` | Table 8; Section 3.8 |
| `SCRIPT_M2c_directions_factorial.R` | `table_M2c_directions_factorial.csv`, `table_M2c_directions_paired.csv` | Table 10 |
| `SCRIPT_M3_adverse_scenarios.R` | `table_M3_all.csv`, `table_M3_paired.csv`, `table_M3_generator_check.csv` | Tables 11 and 12; Section 3.1 (skewness) |
| `SCRIPT_M3b_exact_means.R` | `table_M3_all_6dec.csv` | Tables 11 and 12 (third decimal) |
| `SCRIPT_R1_04_scale_dependence.R` | `table_R1_04_scale_dependence.csv` | Table 9 (unscaled and correlation columns) |
| `SCRIPT_R1_04e_sim_MAD.R` | `table_R1_04e_sim_MAD.csv` | Table 9 (standard deviation and MAD columns) |
| `SCRIPT_R1_05_Sw_properties_v030.R` | `table_R1_05_Sw_properties_v030.csv`, `table_R1_05_Sw_elementwise_bias_v030.csv` | Tables A1 and A2; Sections 2.4, 3.3 and 3.6 |
| `TEP_VARIABLES_SUMMARY.R` | `tep_variables_summary.csv` | Table 15 |
| `19_acf_with_step1.R` | `step_sensitivity_ljungbox_with_step1.csv` | Tables 16 and 17; Section 4.3 |
| `SCRIPT_Figure_ACF_v2.R` | figure only | Figure 4 |
| `PHASE1_PIPELINE_STEP30.R` | no file; builds the batches | Called by the Tennessee Eastman scripts |
| `SCRIPT_M4_TEP_final_mstar14.R` | `table_M4_final_A_limits.csv`, `table_M4_final_B_idv7.csv`, `table_M4_final_C_faults.csv`, `table_M4_final_D_step.csv`, `table_M4_final_E_center.csv` | Tables 18 and 19; Table 17, false alarm column; Section 4.5 |
| `SCRIPT_M4b_TEP_D2_D4.R` | `table_M4b_D2_D4.csv`, `table_M4b_D3_summary.csv`, `table_M4b_D3_compositions.csv` | Sections 4.4 and 4.5 (determinant ratio, twenty compositions) |
| `SCRIPT_M4c_TEP_center_weights.R` | `table_M4c_center_by_variable.csv`, `table_M4c_weights.csv`, `table_M4c_center_membership.csv` | Sections 4.4 and 4.7 |
| `SCRIPT_M4d_TEP_competitors_verify.R` | `table_M4d_competitors_TEP.csv`, `table_M4d_mechanism.csv` | Section 4.5 (RMCD and MRCD on the Tennessee Eastman process) |
| `SCRIPT_M4e_TEP_table12_mstar14.R` | `table_M4e_table12.csv` | Table 20 |
| `SCRIPT_R1_04c_TEP_weighting.R` | `table_R1_04c_TEP_weighting.csv` | Section 4.4 (dependence on the units in the basic combination) |
| `SCRIPT_V030_verify_TEP.R` | `table_R1_04f_TEP_final_units.csv` | Section 4.4 (the final method does not depend on the units) |
| `SCRIPT_R3_02_computation_time_v030.R` | `table_R3_02_computation_time_v030.csv`, `sessionInfo_R3_02_v030.txt` | Section 2.5 and Discussion (computation time) |
| `SCRIPT_Fig2_power_final.R` | figure only | Figure 2 |
| `SCRIPT_Fig5_weights_final.R` | figure only | Figure 5 |
| `SCRIPT_Fig6_control_charts_final.R` | `table_Fig6_T2_batches.csv` | Figure 6; Section 4.5 |
| `SCRIPT_V_ucl_simulated_vs_M5.R` | no file | Checks that `ucl_simulated()` in the package reproduces the simulated limit of Table 5 |

Figures 1 and 3 are diagrams and do not come from any script. Table 14
counts the comparisons in `table_M2_paired.csv`, `table_M2b_paired.csv`,
`table_M2c_directions_paired.csv` and `table_M3_paired.csv`.

## Notes on the result files

Five files in `02_results/` are results of version 0.2.2 and do not support
any table of the revised paper: `table_3_4_sensitivity.csv`,
`table_3_4b_contamination.csv`, `table_3_4c_rho.csv`,
`table_ablation_2x2_final.csv` and `phase2_stability_by_composition.csv`.
They are kept because M2 reads them to check that the new code reproduces the
numbers already published, and M4b takes the twenty Phase 2 compositions from
the last one.

`table_M5_TEP.csv` is a by-product of M5, computed with m* = 13, and is not
cited in the paper.

A few numbers in the text are simple operations on these files, for example
the ratio between the first eigenvalue and the mean of the other three in
Table A1.

## Notes on running the scripts

The run order is M1, M1b, M1c; M5, M5b; R1_10, M2, M2M3_diagnostics, M2b,
M2c; M3, M3b; R1_04, R1_04e; R1_05. For the Tennessee Eastman process:
TEP_VARIABLES_SUMMARY, 19_acf_with_step1, M4, M4b, M4c, M4d, M4e, R1_04c,
V030_verify_TEP and R3_02. The figures come last. V_ucl_simulated_vs_M5
needs the T² values of M5.

`19_acf_with_step1.R` compares its results against
`step_sensitivity_ljungbox.csv`, a version 0.2.2 file that is no longer
shipped; if it is not found, the script issues a warning and carries on.

The comments of some scripts cite the table numbering of an earlier draft.
The table above uses the numbering of the paper.

## Label values in the result files

The analysis was carried out in Spanish and some labels remain in that
language, and the methods appear under their working names. They are labels
only: every numeric column is independent of them. The correspondence with
the paper is as follows.

| In the files | In the paper |
|---|---|
| `a_classical` | (a) classical chart |
| `b_MCD_unif` | (b) MCD only |
| `c_AFM_cls` | (c) MFA only |
| `d_V7`, `V7` | (d) MCD + MFA, basic combination |
| `e_NEW`, `NEW` | (e) MFA-MCD, the proposed method |
| `f_MRCD` | MRCD |
| `g_RMCD_pooled` | RMCD |
| `AFM` | MFA (Multiple Factor Analysis) |
| `D1`, `mD1` | joint increase, joint decrease |
| `D2`, `D3`, `D4`, `mD4` | contrast, block contrast, single increase, single decrease |
| `sano` / `fallo` | fault-free batch / faulty batch |
| `limpia` / `contaminada`, `clean` / `contaminated` | clean / contaminated Phase 1 |
| `publicada`; `sorteo_01` … `sorteo_19` | the Phase 2 composition of Table 18; the nineteen redrawn compositions of Section 4.5 |
| `unidades`; `pesos` | units of measurement; how the weights are computed |
| `ratio_sano_cont` | mean weight of the fault-free batches divided by that of the contaminated ones |
| `fallo_en_6_menores`, `en_6_menores` | contaminated batches among the six lowest-weighted |
| `deteccion_de_20`; `falsas_alarmas_de_10` | faulty batches detected out of 20; false alarms out of 10 |
| `T2_med_sano`; `T2_min_fallo`; `T2_max_sano` | median T² of the fault-free batches; minimum T² of the faulty batches; maximum T² of the fault-free batches |
| `max_dif_pesos` | largest difference in the weights when the units change |
