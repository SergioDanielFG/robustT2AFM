#' Parametric Upper Control Limit (UCL) for AFM-MCD using F Distribution
#'
#' Computes the Upper Control Limit for the Phase 2 Hotelling T-squared statistic
#' using the F-distribution approach with effective sample size correction.
#' This is the "Way 1" UCL proposed by Frutos-Galarza et al. (2026), where the
#' effective batch size m* = round(I * h) replaces I in the degrees of freedom
#' to account for the reduction of sample size due to MCD trimming. The
#' Hardin-Rocke correction is NOT applied.
#'
#' @param calibration A list returned by \code{\link{calibrate_afm_mcd}}.
#' @param I Integer. Number of observations per batch in Phase 1 (e.g., 20).
#' @param alpha Numeric in (0, 1). Nominal false alarm rate. Default 0.001
#'   (European convention, in-control ARL0 = 1000), matching the value used
#'   in Frutos-Galarza et al. (2026). Change to 0.0027 to obtain 3-sigma limits.
#'
#' @return A list containing:
#' \describe{
#'   \item{UCL}{Upper Control Limit value.}
#'   \item{method}{Character string: "F-adjusted (parametric)".}
#'   \item{parameters}{Named list with J, K, I, m_star, alpha used.}
#' }
#'
#' @details
#' The UCL is computed as:
#' \deqn{UCL = \frac{J(K+1)(m^*-1)}{K m^* - K - J + 1} F_{J, K m^* - K - J + 1, 1-\alpha}}
#' where:
#' \itemize{
#'   \item J = number of variables
#'   \item K = number of Phase 1 batches
#'   \item m* = round(I * mcd_alpha), effective batch size after MCD
#'   \item F = F-distribution quantile
#' }
#'
#' @references
#' Montgomery, D. C. (2009). Introduction to Statistical Quality Control,
#' 6th edition. John Wiley & Sons, Hoboken, NJ. ISBN 978-0-470-16992-6.
#'
#' @importFrom stats qf
#' @export
#'
#' @examples
#' # Full pipeline: calibrate Phase 1, compute the operational UCL,
#' # monitor Phase 2, and count out-of-control batches.
#' sim    <- simulate_batch_process(
#'   K1 = 30, K2 = 20, I = 20, J = 4, rho = 0.6,
#'   outlier_batches_F1 = 6, outlier_rate = 0.20, outlier_shift = 4,
#'   prop_ooc_F2 = 0.5, shift_ooc = 1.0, seed = 20260425
#' )
#' vars   <- paste0("Var", 1:4)
#' cal    <- calibrate_afm_mcd(subset(sim, Phase == "Phase 1"), vars)
#'
#' # F-adjusted UCL (operational default in Frutos-Galarza et al. 2026).
#' ucl <- ucl_F_adjusted(cal, I = 20)
#' cat(sprintf("UCL = %.2f  (m* = %d, df2 = %d, alpha = %.3f)\n",
#'             ucl$UCL, ucl$parameters$m_star,
#'             ucl$parameters$df2, ucl$parameters$alpha))
#'
#' # Apply against Phase 2: pass the limit to monitor_afm_mcd() and it
#' # flags the out-of-control batches for you.
#' mon <- monitor_afm_mcd(subset(sim, Phase == "Phase 2"), cal, vars,
#'                        ucl = ucl$UCL)
#' sum(mon$is_ooc)                          # number of alarms
ucl_F_adjusted <- function(calibration, I, alpha = 0.001) {

  # --- Input validation ---
  if (!is.list(calibration) || !all(c("Sw", "mcd_alpha", "weights") %in% names(calibration))) {
    stop("'calibration' must be a list from calibrate_afm_mcd().")
  }
  if (!is.numeric(I) || length(I) != 1 || I < 2 || I != round(I)) {
    stop("'I' must be a single positive integer >= 2.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value in (0, 1).")
  }

  # --- Defensive check on the batch size ---
  # Con lotes desiguales no hay un I correcto que comprobar: el limite es una
  # aproximacion se pase el que se pase, y eso es lo que hay que decir.
  sizes <- calibration$batch_sizes
  if (!is.null(sizes) && length(unique(as.integer(sizes))) > 1L) {
    warning("Phase 1 batch sizes were not equal (", min(sizes), " to ",
            max(sizes), "), so this limit is an approximation whichever I is ",
            "passed. Computed with I = ", I, ".", call. = FALSE)
  } else if (!is.null(calibration$I_phase1) && calibration$I_phase1 != I) {
    warning("'I' (", I, ") does not match the Phase 1 batch size recorded ",
            "during calibration (I_phase1 = ", calibration$I_phase1, "). ",
            "m* will be computed with I = ", I, ".")
  }

  # --- Extract parameters from calibration ---
  J <- nrow(calibration$Sw)            # number of variables
  K <- length(calibration$weights)     # number of Phase 1 batches
  h <- calibration$mcd_alpha            # MCD retention fraction
  m_star <- round(I * h)                # effective batch size after MCD

  # --- Check degrees of freedom feasibility ---
  df2 <- K * m_star - K - J + 1
  if (df2 < 1) {
    stop("Degrees of freedom (df2 = K*m* - K - J + 1) must be >= 1. ",
         "Got df2 = ", df2, " with K=", K, ", m*=", m_star, ", J=", J, ". ",
         "Consider increasing I or K.")
  }

  # --- Compute F-adjusted UCL ---
  scale_factor <- (J * (K + 1) * (m_star - 1)) / df2
  F_quantile <- stats::qf(1 - alpha, df1 = J, df2 = df2)
  UCL <- scale_factor * F_quantile

  return(list(
    UCL = UCL,
    method = "F-adjusted (parametric)",
    parameters = list(
      J = J,
      K = K,
      I = I,
      m_star = m_star,
      h = h,
      alpha = alpha,
      df1 = J,
      df2 = df2,
      scale_factor = scale_factor,
      F_quantile = F_quantile
    )
  ))
}
