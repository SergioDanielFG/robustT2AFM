#' Monitor New Batches Using AFM-MCD Calibration (Phase 2)
#'
#' Evaluates new batches against a reference calibration obtained from
#' \code{\link{calibrate_afm_mcd}}. For each new batch, computes its simple
#' mean (without MCD) and the Hotelling T-squared statistic using the frozen
#' reference center (mu_r) and AFM-weighted covariance matrix (Sw) from Phase 1.
#' The T-squared statistics can be compared against an Upper Control Limit (UCL)
#' to identify out-of-control batches.
#'
#' @param new_data A data frame containing the new batches to evaluate. Must
#'   contain a column named 'Batch' identifying each batch.
#' @param calibration A list returned by \code{\link{calibrate_afm_mcd}}
#'   containing mu_r and Sw from Phase 1 calibration.
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
#' \deqn{T^2 = I \cdot (\bar{x}_{new} - \mu_r)^T S_w^{-1} (\bar{x}_{new} - \mu_r)}
#' where \eqn{\bar{x}_{new}} is the simple mean of the new batch (not MCD),
#' \eqn{\mu_r} is the reference center from Phase 1, \eqn{S_w} is the
#' AFM-weighted covariance matrix from Phase 1, and I is the number of
#' observations in the new batch.
#'
#' In Phase 2, MCD is NOT applied to new batches. The new batch enters the
#' monitoring as-is, and is compared against the robust structure established
#' in Phase 1.
#'
#' @references
#' Montgomery, D. C. (2019). Introduction to Statistical Quality Control
#' (8th ed.). Wiley. Chapter 11: Multivariate Process Monitoring and Control.
#'
#' @export
#'
#' @examples
#' # Simulate a mixed Phase 1 (contamination) + Phase 2 (30% OOC) scenario.
#' sim <- simulate_batch_process(
#'   K1 = 30, K2 = 20, I = 20, J = 4,
#'   outlier_batches_F1 = 2, prop_contam_F1 = 0.07,
#'   prop_ooc_F2 = 0.30, shift_ooc = 2,
#'   seed = 20260417
#' )
#' vars   <- paste0("Var", 1:4)
#' phase1 <- subset(sim, Phase == "Phase 1")
#' phase2 <- subset(sim, Phase == "Phase 2")
#'
#' # Phase 1 calibration and Phase 2 monitoring.
#' cal <- calibrate_afm_mcd(phase1, vars)
#' mon <- monitor_afm_mcd(phase2, cal, vars)
#'
#' # Flag out-of-control batches against the operational F-adjusted UCL.
#' ucl        <- ucl_F_adjusted(cal, I = 20)$UCL
#' mon$is_ooc <- mon$T2 > ucl
#' mon[mon$is_ooc, ]                        # batches to investigate
monitor_afm_mcd <- function(new_data, calibration, variables) {

  # --- Input validation ---
  if (!is.data.frame(new_data)) {
    stop("'new_data' must be a data frame.")
  }
  if (!"Batch" %in% colnames(new_data)) {
    stop("'new_data' must contain a column named 'Batch'.")
  }
  if (!is.list(calibration) || !all(c("mu_r", "Sw") %in% names(calibration))) {
    stop("'calibration' must be a list from calibrate_afm_mcd() ",
         "containing 'mu_r' and 'Sw'.")
  }
  if (!all(variables %in% colnames(new_data))) {
    missing_vars <- setdiff(variables, colnames(new_data))
    stop("Variables not found in new_data: ", paste(missing_vars, collapse = ", "))
  }
  if (length(variables) != length(calibration$mu_r)) {
    stop("Number of variables (", length(variables), ") does not match ",
         "calibration mu_r dimension (", length(calibration$mu_r), ").")
  }

  # --- Setup ---
  mu_r <- calibration$mu_r
  Sw_inv <- solve(calibration$Sw)
  batches <- unique(new_data$Batch)

  # --- Compute T-squared per batch ---
  results <- data.frame(
    Batch = character(),
    I = integer(),
    T2 = numeric(),
    stringsAsFactors = FALSE
  )

  for (batch in batches) {
    subset_batch <- new_data[new_data$Batch == batch, variables]
    I <- nrow(subset_batch)

    # Simple mean (no MCD in Phase 2)
    x_bar <- colMeans(subset_batch)

    # T-squared statistic
    diff <- x_bar - mu_r
    T2 <- as.numeric(I * t(diff) %*% Sw_inv %*% diff)

    results <- rbind(results, data.frame(
      Batch = as.character(batch),
      I = I,
      T2 = T2,
      stringsAsFactors = FALSE
    ))
  }

  return(results)
}
