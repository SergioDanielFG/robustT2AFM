# robustT2AFM 0.3.0

Changes made in response to the first round of review of Frutos-Galarza et al.
(2026). The calibration of version 0.2.0 is still available and gives exactly
the same results as before:
`calibrate_afm_mcd(..., scaling = "none", center = "mean")`, and so is its
control limit: `ucl_F_adjusted(..., m_star = "nominal")`.

## New in `calibrate_afm_mcd()` and `run_afm_mcd()`

* `scaling = "mad"` (new default). Before the first eigenvalue of each batch
  is computed, every variable is divided by its median absolute deviation over
  all Phase 1 observations. The AFM weights no longer depend on the units of
  measurement: multiplying a variable by any constant leaves the weights and
  the T2 statistics unchanged. The divisor is common to all batches, so a batch
  with inflated dispersion still receives a small weight. `Sw` and `mu_r` are
  still returned in the original units. (Reviewer 1, comment 4; Reviewer 3,
  note 1.)

* `center = "mcd"` (new default). The reference center is the reweighted MCD
  location of the K batch centers (`robustbase::covMcd`, `center_alpha = 0.5`)
  instead of their simple average. It resists whole batches shifted in mean,
  which pulled the average and raised the false-alarm rate, and it is affine
  equivariant, unlike the spatial median. (Reviewer 1, comments 6 and 7.)
  It runs under a fixed internal seed, so the same Phase 1 always gives the
  same center; the user's random number stream is restored afterwards. It
  needs at least 2J valid batches; with fewer, the function stops and asks
  for `center = "mean"`.

* The calibration now returns `scale`, `scaling`, `center` and `center_alpha`,
  and `summary()` of a study prints the options used.

## New in `ucl_F_adjusted()` and `run_afm_mcd()`

* `m_star = "subset"` (new default). The effective batch size m* in the
  degrees of freedom of the limit is now the size of the subset on which the
  MCD of each batch is computed, `robustbase::h.alpha.n(mcd_alpha, I, J)`
  (the `quan` of `covMcd`): 14 for I = 20, J = 4 and `mcd_alpha = 0.67`.
  Version 0.2.0 used `round(I * mcd_alpha)` = 13, the nominal fraction, which
  is not the number of observations the algorithm works with; it is kept as
  `m_star = "nominal"`. The limit decreases slightly (more degrees of
  freedom). `parameters$m_star_rule` records the rule used, and `summary()`
  prints it. With unequal Phase 1 batches, `I_phase1` is now the smallest size
  whose m* equals the rounded mean of the per-batch m*. (Reviewer 1,
  comments 3 and 14.)

* The whole study of version 0.2.0 is reproduced with
  `run_afm_mcd(..., scaling = "none", center = "mean", m_star = "nominal")`.

## New function `ucl_simulated()`

* A simulated control limit, offered as an alternative to the analytic one of
  `ucl_F_adjusted()`. The whole method, calibration included, is applied to
  `B` synthetic Phase 1 samples drawn from a normal distribution with the
  robust estimates `mu_r` and `Sw`; the limit is the 1 - alpha quantile of
  the T2 of new in-control batches. In the simulation study it brought ARL0
  close to the nominal value with a clean Phase 1, also with K = 100 batches,
  where the analytic limit is conservative; with a contaminated Phase 1 it
  gave the same ARL0 as the analytic limit. `seed` makes it reproducible
  without touching the user's random number stream. Takes about two minutes
  with the defaults for K = 30. (Reviewer 1, comment 2.)

## Unchanged

* The per-batch MCD (reweighted, `mcd_alpha = 0.67`), the weight formula, `Sw`,
  `monitor_afm_mcd()` and the form of the control limit (Equation 8).

# robustT2AFM 0.2.0

* Version accompanying the manuscript as first submitted (tag `v0.2.0`).
