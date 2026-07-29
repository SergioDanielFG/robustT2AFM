# robustT2AFM

AFM-Weighted Robust T-Squared Control Chart for Batch Processes

## Overview

This package implements a robust multivariate statistical process control
chart for batch manufacturing. The method combines:

- Minimum Covariance Determinant (MCD) estimation per batch
- Multiple Factor Analysis (AFM) weighting of batch covariances (w = 1/lambda1)
- Parametric F-adjusted upper control limit (`ucl_F_adjusted`) computed with
  the effective batch size m* = round(I * h). This is the sole UCL exposed
  by the package; the paper discussed a bootstrap alternative only to
  discard it, so it is not implemented here
- Publication-quality control chart with out-of-control batches highlighted
  and labelled for quality engineers, and one-line export to PNG + PDF

The method preserves sensitivity to mean shifts when Phase I calibration
data contains outlier contamination - a scenario where classical Hotelling
T-squared suffers from masking effects.

## Minimum workflow

```r
library(robustT2AFM)

# 1. Historical data with a "Batch" column and J numeric variables.
sim    <- simulate_batch_process(K1 = 30, K2 = 20, I = 20, J = 4,
                                 prop_ooc_F2 = 0.3, seed = 20260417)
vars   <- paste0("Var", 1:4)

# 2. Calibrate Phase 1.
cal    <- calibrate_afm_mcd(subset(sim, Phase == "Phase 1"), vars)

# 3. Compute the operational UCL.
ucl    <- ucl_F_adjusted(cal, I = 20)

# 4. Monitor Phase 2 and plot the control chart.
mon    <- monitor_afm_mcd(subset(sim, Phase == "Phase 2"), cal, vars)
plot_control_chart(mon, UCL = ucl$UCL,
                   method_label = "AFM-MCD Control Chart - Phase 2",
                   alpha = ucl$parameters$alpha)
```

## Installation

remotes::install_github("SergioDanielFG/robustT2AFM")

## Citation

Frutos Galarza, S. D., Ruiz Barzola, O., Ramirez, J., Galindo Villardon, P.
(2026). AFM-Weighted Robust T-Squared Control Chart with Bootstrap UCL
for Batch Processes under Phase I Contamination. Manuscript in preparation.

## License

MIT (c) 2026 Sergio Daniel Frutos Galarza, Omar Ruiz Barzola,
Jhon Ramirez, Purificacion Galindo Villardon
