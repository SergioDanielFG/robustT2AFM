# Design decisions

Decisions that cannot be inferred from reading the code, together with the
measurements that support them. Each section answers a question that has been
asked more than once, almost always from the wrong direction.

---

## 1. Masking makes the classical chart deaf, not jumpy

This is the easiest thing in the whole method to get backwards, and it is worth
writing down because the natural intuition points the wrong way.

The wrong intuition says: if Phase 1 is contaminated, the classical chart will
fire constantly. The actual chain runs the other way:

1. Contamination in Phase 1 inflates the pooled covariance `Sp`. On the
   Tennessee Eastman data its determinant is **16.14 times** that of the robust
   reference matrix `Sw`.
2. The classical control limit **does not move**. `hotelling_classical_ucl()`
   is parametric and depends only on K, I, J and alpha: it returns 19.4644 with
   a clean Phase 1 and with a contaminated one.
3. What changes is the statistic. An inflated `Sp` shrinks every batch T² at
   once, and shrunken statistics never reach a limit that has stayed put.

Contamination does not make the classical chart jumpy. It makes it deaf. What
any figure or passage should show is **lost detection**.

The acceptance ellipsoid does widen in data space, which is where the loose
intuition comes from. But the UCL is a specific number and that number does not
change. See also section 10.

### The measurements

On `afm_phase2`, the T² of `F2_B01` drops from 39.19 under clean calibration to
16.92 under contaminated calibration, and the classical alarm count falls from
9 of 20 to 2 of 20, against the same limit of 19.4644.

Across the Tennessee Eastman application and the package example scenario, the
classical method produced **zero false alarms** in every case measured:

- Table 9 of the paper.
- The five faults of Table 10.
- The 20 Phase 2 compositions of the stability analysis.
- The 12 seeds and 200 replicates of the package example scenario.

Not "few". Zero.

### Two necessary qualifications

**In simulation the classical method does raise false alarms.** Table 1 of the
paper reports 146 over 200,000 in-control batches, FAR = 0.000730,
ARL0 = 1370. The zero holds for the Tennessee Eastman application and for the
package scenario; it is not a general property of the classical chart.

**And the proposed method is not free of them either.** When Phase 1 is
calibrated on IDV 1 batches, a far more severe fault than IDV 7, the reference
center drifts by 0.82 standard deviations and six of the ten fault-free Phase 2
batches cross the limit. This is Section 4.6 of the paper. The weighting
protects `Sw`, not `mu_r`.

### What follows from this

No figure, help page or vignette should suggest that the classical method
raises more false alarms under contamination. And where the count is zero, that
zero is not a shortcoming of the chosen scenario to be fixed by finding a
different one. Going looking for such a scenario is the inverted intuition all
over again.

This is what `plot_method_comparison()` exists for, and why its caption opens
with the mechanism rather than with a count.
