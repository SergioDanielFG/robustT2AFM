# robustT2AFM

AFM-Weighted Robust T-Squared Control Chart for Batch Processes

## Overview

This package implements a robust multivariate statistical process control
chart for batch manufacturing. The method combines:

- Minimum Covariance Determinant (MCD) estimation per batch
- Multiple Factor Analysis (AFM) weighting of batch covariances (w = 1/lambda1)
- Bootstrap-calibrated upper control limits
- Visual diagnostic panels for out-of-control root cause analysis

The method preserves sensitivity to mean shifts when Phase I calibration
data contains outlier contamination - a scenario where classical Hotelling
T-squared suffers from masking effects.

## Installation

remotes::install_github("SergioDanielFG/robustT2AFM")

## Citation

Frutos Galarza, S. D., Ruiz Barzola, O., Ramirez, J., Galindo Villardon, P.
(2026). AFM-Weighted Robust T-Squared Control Chart with Bootstrap UCL
for Batch Processes under Phase I Contamination. Manuscript in preparation.

## License

MIT (c) 2026 Sergio Daniel Frutos Galarza, Omar Ruiz Barzola,
Jhon Ramirez, Purificacion Galindo Villardon
