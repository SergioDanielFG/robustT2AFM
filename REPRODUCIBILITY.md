# Reproducibility

Numeric anchors, fixed parameters and seeds for the results reported in
Frutos-Galarza et al. (2026). Any change to the package must leave these
values untouched.

## Numeric anchors

With J = 4, K = 30, I = 20, h = 0.67, alpha = 0.001:

| Quantity | Must yield |
|---|---|
| Robust UCL | 19.69285 |
| Classical Phase II UCL | 19.46440 |
| m* = round(I·h) | 13 |
| Tennessee Eastman, proposed method | 20 of 20, 0 false alarms |
| Tennessee Eastman, classical method | 2 of 20 |
| Determinant inflation | 16.14 |

**Identical, not similar.** `covMcd` is deterministic in this configuration,
verified over ten seeds, so a difference in the sixth decimal means something
broke rather than that randomness occurred.

## Standard parameters

Both must be passed explicitly, because the package defaults differ:

- `mcd_alpha = 0.67` — MCD retention fraction
- `alpha = 0.001` — nominal false alarm rate

These are two unrelated parameters with similar names. `mcd_alpha` governs how
much of each batch the robust estimator keeps; `alpha` governs the control
limit.

## Frozen signatures

The published function signatures are fixed. Names, argument order and
defaults do not change, because the results in the paper were produced against
them:

- `ucl_F_adjusted(calibration, I, alpha)$UCL`
- `hotelling_classical_ucl(K, I, J, alpha, phase = "II")$UCL`
- the monitoring functions take `(new_data, calibration, VARS)` in that order

## Before considering any change finished

```r
devtools::test()
devtools::check()
```

All tests must pass and `check()` must report no errors, warnings or notes. If
`check()` is not clean, the change is not finished.

`check()` needs pandoc on the PATH to rebuild the vignette. The NOTE
"unable to verify current time" is a network artefact of `--as-cran`, not a
defect; a run whose only NOTE is that one counts as clean.

## Seeds

Every Monte Carlo campaign derives its streams from a single base,
`SEED_BASE = 2026`, offset by a per-study multiplier so that no two campaigns
draw from overlapping ranges.

| Study | Seed of replicate | Script |
|---|---|---|
| Table 1, false alarm calibration | `2026 * 100 + rep` | `ZTEST_TABLE1_FAITHFUL_v2.R` |
| Table 1, K_eff row | `2026 * 100 + rep` | `KEFF_EFFECT_v2.R` |
| Table 1, bootstrap rows | `2026 * 100 + rep` | `CENTERED_BOOTSTRAP_v5.R` |
| Table 2, ARL0 against K | `2026 + r` | `15_arl0_convergence_K_v3.R` |
| Clean vs contaminated Phase 1 | bases `2026, 20000, 40000, 60000, 80000`, plus `+ r` | `20_arl0_clean_vs_contaminated_v2.R` |
| Tables 3 and 5, detection | `2026 * 200 / *250 / *300 + rep` | `SCRIPT_FINAL_SIMULATIONS.R` |
| Table 4, 2x2 ablation | `2026 * 200 / *250 / *300 + rep` | `SCRIPT_FINAL_ABLATIONS_2x2.R` |
| Table 5, rho block | `2026 * 200 / *250 / *300 + idx` | `22_rho_sensitivity_parallel.R` |
| Phase 2 stability | `2026` | `21_phase2_stability_resampling.R` |
| Reference center audit | `2026 * 1200 + idx` | `23_REFERENCE_CENTER_AUDIT_v4.R` |
| Tennessee Eastman reproducibility | `2026`, then `1001` to `1010` | `TEP_BOOTSTRAP_LIMIT.R` |

**Bases must stay farther apart than the replicate count.** An earlier version
of the clean-versus-contaminated campaign used bases separated by 1, so the
seed ranges of the five campaigns overlapped and the runs were not independent.
The current script checks the separation before it starts and stops if it is
too small.

The package documentation scenario uses a different seed, `20260425`, chosen
because it is the replicate closest to the mean detection rate over 200
replicates. It is not one of the paper seeds and the two should not be
confused.
