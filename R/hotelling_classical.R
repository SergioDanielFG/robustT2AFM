#' Phase 1 Calibration of the Classical Hotelling T-squared Control Chart
#'
#' Estimates the in-control grand mean vector and pooled within-batch covariance
#' matrix from a set of K historical batches in Phase 1. These Phase 1 estimates
#' are used as plug-in parameters for the Phase 2 monitoring statistic computed
#' by \code{\link{hotelling_classical_monitor}}. Unlike
#' \code{\link{calibrate_afm_mcd}}, no robust estimation is applied: the
#' function uses the standard sample mean and pooled within-batch sample
#' covariance, which is the classical estimator described in Montgomery (2009,
#' Chapter 11). This function is included in the package as a benchmark for
#' the Monte Carlo comparisons reported in the paper.
#'
#' @param data A data frame containing Phase 1 process data. Must contain a
#'   column named 'Batch' identifying each batch.
#' @param variables Character vector with the names of the process variables.
#'
#' @return A list containing:
#' \describe{
#'   \item{mu_global}{Grand mean vector (length J), computed as the average of
#'         the K batch means.}
#'   \item{Sp}{Pooled within-batch covariance matrix (J x J). With balanced
#'         batches it is the simple average of the per-batch covariances; with
#'         unbalanced batches it is weighted by the within-batch degrees of
#'         freedom (I_k - 1).}
#'   \item{n_batches}{Number of Phase 1 batches K used for calibration.}
#'   \item{batch_sizes}{Integer vector with the size I_k of each calibration batch.}
#'   \item{variables}{The character vector of variable names, stored for
#'         downstream consistency checks in \code{hotelling_classical_monitor}.}
#' }
#'
#' @details
#' The pooled covariance Sp is the standard estimator of the common within-batch
#' covariance Sigma under the assumption that the process is in statistical
#' control during Phase 1. It is the classical (non-robust) counterpart of the
#' AFM-weighted MCD estimator implemented in \code{\link{calibrate_afm_mcd}}.
#'
#' Notation follows Montgomery (2009), translated to the package convention:
#' J = number of variables (p in Montgomery), K = number of preliminary batches
#' (m in Montgomery), I = batch size (n in Montgomery).
#'
#' @references
#' Hotelling, H. (1947). Multivariate quality control, illustrated by the air
#' testing of sample bombsights. In \emph{Techniques of Statistical Analysis}
#' (pp. 111-184). McGraw-Hill.
#'
#' Tracy, N. D., Young, J. C., & Mason, R. L. (1992). Multivariate control
#' charts for individual observations. \emph{Journal of Quality Technology},
#' 24(2), 88-95.
#'
#' Montgomery, D. C. (2009). \emph{Introduction to Statistical Quality Control},
#' 6th edition. John Wiley & Sons, Hoboken, NJ. ISBN 978-0-470-16992-6. Section 11.3.
#'
#' @importFrom stats cov
#' @export
#'
#' @examples
#' # Classical Hotelling calibration on the same Phase 1 used by the
#' # AFM-MCD method, provided as a benchmark. Phase 1 is intentionally
#' # clean here (the classical estimator is not designed to absorb
#' # outliers or shifted batches; see calibrate_afm_mcd for the robust
#' # counterpart).
#' sim    <- simulate_batch_process(
#'   K1 = 30, K2 = 0, I = 20, J = 4, seed = 20260417
#' )
#' phase1 <- subset(sim, Phase == "Phase 1")
#' vars   <- paste0("Var", 1:4)
#'
#' cal <- hotelling_classical_calibrate(phase1, vars)
#' round(cal$mu_global, 3)
#' round(cal$Sp, 3)
hotelling_classical_calibrate <- function(data, variables) {

  # --- Input validation ---
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  if (!"Batch" %in% colnames(data)) {
    stop("'data' must contain a column named 'Batch'.")
  }
  if (!is.character(variables) || length(variables) < 1) {
    stop("'variables' must be a non-empty character vector.")
  }
  if (!all(variables %in% colnames(data))) {
    missing_vars <- setdiff(variables, colnames(data))
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }
  non_numeric <- variables[!vapply(data[variables], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("The following 'variables' are not numeric: ",
         paste(non_numeric, collapse = ", "))
  }

  # --- Setup ---
  batches <- unique(data$Batch)
  K <- length(batches)
  J <- length(variables)

  if (K < 2) {
    stop("At least 2 batches are required for calibration. ",
         "Only ", K, " batch found.")
  }

  # --- Per-batch summaries ---
  batch_means <- matrix(NA_real_, nrow = K, ncol = J,
                        dimnames = list(NULL, variables))
  batch_covs  <- vector("list", K)
  batch_sizes <- integer(K)

  for (k in seq_len(K)) {
    subset_batch <- data[data$Batch == batches[k], variables, drop = FALSE]
    I_k <- nrow(subset_batch)
    if (I_k < 2) {
      stop("Batch '", batches[k],
           "' has fewer than 2 observations (I_k = ", I_k,
           "); covariance cannot be estimated.")
    }
    if (any(!is.finite(as.matrix(subset_batch)))) {
      stop("Batch '", batches[k], "' contains non-finite values (NA/NaN/Inf).")
    }
    batch_means[k, ] <- colMeans(as.matrix(subset_batch))
    batch_covs[[k]]  <- stats::cov(as.matrix(subset_batch))
    batch_sizes[k]   <- I_k
  }

  # --- Grand mean (average of batch means) ---
  mu_global <- colMeans(batch_means)
  names(mu_global) <- variables

  # --- Pooled within-batch covariance ---
  if (length(unique(batch_sizes)) == 1) {
    # Balanced design: simple average of per-batch covariances
    Sp <- Reduce("+", batch_covs) / K
  } else {
    # Unbalanced design: weight by within-batch degrees of freedom (I_k - 1)
    df_w <- batch_sizes - 1L
    Sp   <- Reduce("+", Map(function(S, w) S * w, batch_covs, df_w)) / sum(df_w)
  }
  dimnames(Sp) <- list(variables, variables)

  # --- Sanity check on Sp ---
  if (any(!is.finite(Sp))) {
    stop("Pooled covariance matrix Sp contains non-finite entries.")
  }
  eig_min <- min(eigen(Sp, symmetric = TRUE, only.values = TRUE)$values)
  if (eig_min <= 0) {
    warning("Pooled covariance matrix Sp is not positive-definite ",
            "(min eigenvalue = ", signif(eig_min, 3),
            "). Phase 2 inversion may be unstable.")
  }

  return(list(
    mu_global   = mu_global,
    Sp          = Sp,
    n_batches   = K,
    batch_sizes = batch_sizes,
    variables   = variables
  ))
}


#' Monitor New Batches Using Classical Hotelling Calibration (Phase 2)
#'
#' Evaluates new batches against a reference calibration obtained from
#' \code{\link{hotelling_classical_calibrate}}. For each new batch, computes
#' its sample mean and the Hotelling T-squared statistic using the frozen
#' grand mean (mu_global) and pooled covariance matrix (Sp) from Phase 1.
#' The T-squared statistics can be compared against an Upper Control Limit (UCL)
#' from \code{\link{hotelling_classical_ucl}} to identify out-of-control batches.
#'
#' @param new_data A data frame containing the new batches to evaluate. Must
#'   contain a column named 'Batch' identifying each batch.
#' @param calibration A list returned by \code{\link{hotelling_classical_calibrate}}
#'   containing mu_global and Sp from Phase 1 calibration.
#' @param variables Character vector with the names of the process variables.
#'   Must match the variables used in the calibration.
#'
#' @return A data frame with one row per batch and columns:
#' \describe{
#'   \item{Batch}{Batch identifier.}
#'   \item{I}{Number of observations in the batch.}
#'   \item{T2}{Hotelling T-squared statistic for the batch.}
#' }
#'
#' @details
#' The Phase 2 T-squared statistic for each new batch is computed as:
#' \deqn{T^2 = I \cdot (\bar{x}_{new} - \mu_{global})^T S_p^{-1} (\bar{x}_{new} - \mu_{global})}
#' where \eqn{\bar{x}_{new}} is the sample mean of the new batch,
#' \eqn{\mu_{global}} is the grand mean from Phase 1, \eqn{S_p} is the
#' pooled covariance matrix from Phase 1, and I is the number of
#' observations in the new batch.
#'
#' This is the classical (non-robust) Hotelling T-squared, included in the
#' package as a benchmark for comparisons against the AFM-MCD method
#' implemented in \code{\link{monitor_afm_mcd}}. When Phase 1 batches contain
#' outliers, Sp is inflated and this procedure loses detection power; the
#' AFM-MCD counterpart addresses this limitation.
#'
#' @references
#' Hotelling, H. (1947). Multivariate quality control, illustrated by the air
#' testing of sample bombsights. In \emph{Techniques of Statistical Analysis}
#' (pp. 111-184). McGraw-Hill.
#'
#' Montgomery, D. C. (2009). \emph{Introduction to Statistical Quality Control},
#' 6th edition. John Wiley & Sons, Hoboken, NJ. ISBN 978-0-470-16992-6. Section 11.3, equation 11.19.
#'
#' @export
#'
#' @examples
#' # Classical Hotelling Phase 2 monitoring, shown side-by-side with the
#' # UCL and the OOC flag a quality engineer would use in production.
#' sim  <- simulate_batch_process(
#'   K1 = 30, K2 = 20, I = 20, J = 4,
#'   prop_ooc_F2 = 0.30, shift_ooc = 2, seed = 20260417
#' )
#' vars <- paste0("Var", 1:4)
#' cal  <- hotelling_classical_calibrate(
#'   subset(sim, Phase == "Phase 1"), vars
#' )
#' mon  <- hotelling_classical_monitor(
#'   subset(sim, Phase == "Phase 2"), cal, vars
#' )
#' ucl        <- hotelling_classical_ucl(K = 30, I = 20, J = 4)$UCL
#' mon$is_ooc <- mon$T2 > ucl
#' mon
hotelling_classical_monitor <- function(new_data, calibration, variables) {

  # --- Input validation ---
  if (!is.data.frame(new_data)) {
    stop("'new_data' must be a data frame.")
  }
  if (!"Batch" %in% colnames(new_data)) {
    stop("'new_data' must contain a column named 'Batch'.")
  }
  if (!is.list(calibration) || !all(c("mu_global", "Sp") %in% names(calibration))) {
    stop("'calibration' must be a list from hotelling_classical_calibrate() ",
         "containing 'mu_global' and 'Sp'.")
  }
  if (!all(variables %in% colnames(new_data))) {
    missing_vars <- setdiff(variables, colnames(new_data))
    stop("Variables not found in new_data: ", paste(missing_vars, collapse = ", "))
  }
  if (length(variables) != length(calibration$mu_global)) {
    stop("Number of variables (", length(variables), ") does not match ",
         "calibration mu_global dimension (", length(calibration$mu_global), ").")
  }

  # --- Setup ---
  mu_global <- calibration$mu_global
  Sp_inv <- tryCatch(
    solve(calibration$Sp),
    error = function(e) {
      stop("Pooled covariance matrix is singular and cannot be inverted: ",
           conditionMessage(e))
    }
  )
  batches <- unique(new_data$Batch)

  # --- Compute T-squared per batch ---
  results <- data.frame(
    Batch = character(),
    I = integer(),
    T2 = numeric(),
    stringsAsFactors = FALSE
  )

  for (batch in batches) {
    subset_batch <- new_data[new_data$Batch == batch, variables, drop = FALSE]
    I <- nrow(subset_batch)

    if (I < 1) {
      stop("Batch '", batch, "' has zero observations.")
    }

    # Sample mean (no MCD: classical Hotelling)
    x_bar <- if (I == 1) {
      as.numeric(subset_batch[1, ])
    } else {
      colMeans(as.matrix(subset_batch))
    }

    # T-squared statistic
    diff <- x_bar - mu_global
    T2 <- as.numeric(I * t(diff) %*% Sp_inv %*% diff)

    results <- rbind(results, data.frame(
      Batch = as.character(batch),
      I = I,
      T2 = T2,
      stringsAsFactors = FALSE
    ))
  }

  return(results)
}


#' Upper Control Limit for the Classical Hotelling T-squared Chart
#'
#' Computes the exact F-distribution-based Upper Control Limit for the
#' Hotelling T-squared chart with subgrouped data, in either the Phase 1
#' (retrospective) or Phase 2 (prospective) form derived by Tracy, Young, and
#' Mason (1992) and tabulated by Montgomery (2009). This is the classical
#' (non-robust) counterpart of \code{\link{ucl_F_adjusted}}, included as a
#' benchmark for the Monte Carlo comparisons reported in the paper.
#'
#' @param K Integer. Number of preliminary batches used in Phase 1 calibration.
#'   Must be at least 2.
#' @param I Integer. Number of observations per batch. Must be at least 2.
#' @param J Integer. Number of monitored variables. Must be at least 1.
#' @param alpha Numeric in (0, 1). Family-wise false-alarm probability. Default
#'   0.001 (European convention, in-control ARL0 = 1000), used throughout the
#'   package for consistency with \code{\link{ucl_F_adjusted}} and matching
#'   Frutos-Galarza et al. (2026). Change to 0.0027 to obtain 3-sigma limits.
#' @param phase Character. Either "I" (retrospective limits, Montgomery
#'   eq. 11.20) or "II" (prospective limits, Montgomery eq. 11.21). Default "II".
#'
#' @return A list containing:
#' \describe{
#'   \item{UCL}{Upper Control Limit value.}
#'   \item{method}{Character string: "Hotelling classical (Phase II)" or
#'         "Hotelling classical (Phase I)".}
#'   \item{parameters}{Named list with J, K, I, alpha, phase, df1, df2,
#'         scale_factor, F_quantile used.}
#' }
#'
#' @details
#' The exact UCLs for the Hotelling T-squared statistic on subgrouped data
#' are computed as (Montgomery, 2009, Section 11.3):
#'
#' Phase II (prospective, equation 11.21, default):
#' \deqn{UCL = \frac{J(K+1)(I-1)}{K I - K - J + 1} F_{J, K I - K - J + 1, 1-\alpha}}
#'
#' Phase I (retrospective, equation 11.20):
#' \deqn{UCL = \frac{J(K-1)(I-1)}{K I - K - J + 1} F_{J, K I - K - J + 1, 1-\alpha}}
#'
#' where:
#' \itemize{
#'   \item J = number of variables
#'   \item K = number of Phase 1 batches
#'   \item I = number of observations per batch
#'   \item F = F-distribution quantile
#' }
#'
#' Note that the Phase II UCL equals the Phase I UCL multiplied by (K+1)/(K-1).
#' The denominator degrees of freedom df2 = K I - K - J + 1 must be at least 1
#' for the F quantile to be defined; the function aborts with an informative
#' error otherwise. Lowry and Montgomery (1995) and Jensen et al. (2006)
#' recommend at least K = 20 to 50 preliminary batches for the Phase II UCL
#' to be well estimated.
#'
#' @references
#' Tracy, N. D., Young, J. C., & Mason, R. L. (1992). Multivariate control
#' charts for individual observations. \emph{Journal of Quality Technology},
#' 24(2), 88-95.
#'
#' Montgomery, D. C. (2009). \emph{Introduction to Statistical Quality Control},
#' 6th edition. John Wiley & Sons, Hoboken, NJ. ISBN 978-0-470-16992-6. Section 11.3, equations 11.20 and 11.21.
#'
#' Lowry, C. A., & Montgomery, D. C. (1995). A review of multivariate control
#' charts. \emph{IIE Transactions}, 27(6), 800-810.
#'
#' @importFrom stats qf
#' @export
#'
#' @examples
#' # Phase II UCL with package defaults (K = 30, I = 20, J = 4).
#' ucl_phII <- hotelling_classical_ucl(K = 30, I = 20, J = 4)
#' cat(sprintf("Phase II UCL = %.2f (alpha = %.3f)\n",
#'             ucl_phII$UCL, ucl_phII$parameters$alpha))
#'
#' # Phase I UCL for the same calibration set (retrospective).
#' ucl_phI <- hotelling_classical_ucl(K = 30, I = 20, J = 4, phase = "I")
#' cat(sprintf("Phase I UCL  = %.2f\n", ucl_phI$UCL))
hotelling_classical_ucl <- function(K, I, J, alpha = 0.001, phase = "II") {

  # --- Input validation ---
  if (!is.numeric(K) || length(K) != 1 || K < 2 || K != round(K)) {
    stop("'K' must be a single positive integer >= 2.")
  }
  if (!is.numeric(I) || length(I) != 1 || I < 2 || I != round(I)) {
    stop("'I' must be a single positive integer >= 2.")
  }
  if (!is.numeric(J) || length(J) != 1 || J < 1 || J != round(J)) {
    stop("'J' must be a single positive integer >= 1.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value in (0, 1).")
  }
  if (!is.character(phase) || length(phase) != 1 || !(phase %in% c("I", "II"))) {
    stop("'phase' must be either \"I\" or \"II\".")
  }

  # --- Degrees of freedom ---
  df1 <- J
  df2 <- K * I - K - J + 1

  # --- Check degrees of freedom feasibility ---
  if (df2 < 1) {
    stop("Degrees of freedom (df2 = K*I - K - J + 1) must be >= 1. ",
         "Got df2 = ", df2, " with K=", K, ", I=", I, ", J=", J, ". ",
         "Consider increasing K or I, or reducing J.")
  }

  # --- Compute Hotelling UCL ---
  scale_factor <- if (phase == "II") {
    (J * (K + 1) * (I - 1)) / df2          # Montgomery (2009), eq. 11.21
  } else {
    (J * (K - 1) * (I - 1)) / df2          # Montgomery (2009), eq. 11.20
  }

  F_quantile <- stats::qf(1 - alpha, df1 = df1, df2 = df2)
  UCL <- scale_factor * F_quantile

  method_label <- if (phase == "II") {
    "Hotelling classical (Phase II)"
  } else {
    "Hotelling classical (Phase I)"
  }

  return(list(
    UCL = UCL,
    method = method_label,
    parameters = list(
      J = J,
      K = K,
      I = I,
      alpha = alpha,
      phase = phase,
      df1 = df1,
      df2 = df2,
      scale_factor = scale_factor,
      F_quantile = F_quantile
    )
  ))
}
