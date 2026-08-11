# robustT2AFM

AFM-Weighted Robust T-Squared Control Chart for Batch Processes

## Overview

This package implements a robust multivariate statistical process control
chart for batch manufacturing. The method combines:

- Minimum Covariance Determinant (MCD) estimation per batch
- Multiple Factor Analysis (AFM) weighting of batch covariances (w = 1/lambda1)
- Parametric F-adjusted upper control limit (`ucl_F_adjusted`) computed with
  the effective batch size m* = round(I * h). The paper discussed a bootstrap
  alternative only to discard it, so it is not implemented here. The classical
  Hotelling limit is also available, through `hotelling_classical_ucl()`, as
  the non-robust benchmark.
- Publication-quality control chart with out-of-control batches highlighted
  and labelled for quality engineers, and one-line export to PNG + PDF

The method preserves sensitivity to mean shifts when Phase I calibration
data contains outlier contamination - a scenario where classical Hotelling
T-squared suffers from masking effects.

## Scope

The chart detects shifts of the mean vector. Changes in the covariance
structure require different statistics and are out of scope.

The AFM weights protect the reference covariance further than they protect
the reference centre: they act on dispersion, not on position, so a
calibration batch that is wholly displaced keeps its weight. On the Tennessee
Eastman data a mild displacement moved the centre by 0.06 standard deviations
and monitoring was unaffected; a severe one moved it by 0.82 and six of ten
fault-free batches then crossed the limit. Calibrate on batches you have
reason to believe were on target.

## Minimum workflow

```r
library(robustT2AFM)

# The base configuration of the paper ships with the package: 6 of the
# 30 Phase 1 batches carry outlying observations, and half the Phase 2
# batches are faulty, shifted by 1 sigma.
data(afm_phase1)
data(afm_phase2)

# The whole study in one call.
study <- run_afm_mcd(afm_phase1, afm_phase2)
study            # how many batches are out of control, and which
summary(study)   # the full report

# Or the four steps separately.
vars <- paste0("Var", 1:4)
cal  <- calibrate_afm_mcd(afm_phase1, vars)
ucl  <- ucl_F_adjusted(cal, I = 20)
mon  <- monitor_afm_mcd(afm_phase2, cal, vars, ucl = ucl$UCL)
plot_control_chart(mon, UCL = ucl$UCL,
                   method_label = "AFM-MCD Control Chart - Phase 2",
                   alpha = ucl$parameters$alpha)

# If your batch column is not called "Batch":
#   calibrate_afm_mcd(my_data, vars, batch_col = "Lote")

# Use simulate_batch_process() to build other scenarios.
```

## Installation

remotes::install_github("SergioDanielFG/robustT2AFM")

## Citation

Frutos-Galarza, S. D., Ruiz-Barzola, O., Ramírez, J., Galindo-Villardón, P.
(2026). A Robust Hotelling-Type T2 Control Chart Combining the Minimum
Covariance Determinant Estimator with Multiple Factor Analysis Weighting.
Manuscript in preparation.

## License

MIT (c) 2026 Sergio Daniel Frutos-Galarza, Omar Ruiz-Barzola,
Jhon Ramírez, Purificación Galindo-Villardón
