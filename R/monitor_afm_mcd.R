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
#'   contain a column identifying each batch; see \code{batch_col}.
#' @param calibration A list returned by \code{\link{calibrate_afm_mcd}}
#'   containing mu_r and Sw from Phase 1 calibration.
#' @param variables Character vector with the names of the process variables.
#'   Must match the variables used in the calibration.
#' @param ucl Optional single positive number: the Upper Control Limit to
#'   compare each T-squared against, normally \code{ucl_F_adjusted(cal, I)$UCL}.
#'   Pass the number itself, not the list returned by
#'   \code{\link{ucl_F_adjusted}}. Default \code{NULL}, in which case the
#'   output has exactly the three columns described below.
#' @param batch_col Character. Name of the column that identifies the batch in
#'   \code{new_data}. Default \code{"Batch"}. This renames only what is read:
#'   the returned data frame always calls its first column \code{Batch},
#'   because the plotting functions and the study object depend on that name.
#'
#' @return A data frame with one row per batch and columns:
#' \describe{
#'   \item{Batch}{Batch identifier.}
#'   \item{I}{Number of observations in the batch.}
#'   \item{T2}{Hotelling T-squared statistic for the batch.}
#'   \item{is_ooc}{Logical, present only when \code{ucl} was supplied.
#'         \code{TRUE} marks an out-of-control batch, i.e. one whose T-squared
#'         is strictly greater than the limit. A batch landing exactly on the
#'         limit is not flagged.}
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
#' Computing the statistic requires inverting \eqn{S_w}. If that matrix is
#' singular the function aborts with an explanatory error instead of the raw
#' \code{solve()} message, mirroring \code{\link{hotelling_classical_monitor}}.
#' A singular \eqn{S_w} usually means two monitored variables are redundant,
#' or that the Phase 1 batches were too small for their number of variables.
#'
#' Supplying \code{ucl} only adds the \code{is_ooc} column; the T-squared
#' values are computed identically with and without it.
#'
#' @references
#' Montgomery, D. C. (2009). Introduction to Statistical Quality Control,
#' 6th edition. John Wiley & Sons, Hoboken, NJ. ISBN 978-0-470-16992-6. Chapter 11: Multivariate Process Monitoring and Control.
#'
#' @export
#'
#' @examples
#' data(afm_phase1)
#' data(afm_phase2)
#' vars <- paste0("Var", 1:4)
#'
#' # Phase 1 calibration, operational UCL and Phase 2 monitoring.
#' cal <- calibrate_afm_mcd(afm_phase1, vars)
#' ucl <- ucl_F_adjusted(cal, I = 20)$UCL
#' mon <- monitor_afm_mcd(afm_phase2, cal, vars, ucl = ucl)
#'
#' mon[mon$is_ooc, ]                        # batches to investigate
#' sum(mon$is_ooc)                          # how many alarms fired
#'
#' # Without 'ucl' you get the T-squared values only, and compare them
#' # yourself against whatever limit you prefer.
#' mon_plain <- monitor_afm_mcd(afm_phase2, cal, vars)
#' names(mon_plain)
monitor_afm_mcd <- function(new_data, calibration, variables, ucl = NULL,
                            batch_col = "Batch") {

  # --- Input validation ---
  if (!is.data.frame(new_data)) {
    stop("'new_data' must be a data frame.")
  }
  check_batch_col(
    new_data, batch_col, "new_data",
    "monitor_afm_mcd(new_data, calibration, variables, batch_col = \"<name>\")"
  )
  if (!is.list(calibration) || !all(c("mu_r", "Sw") %in% names(calibration))) {
    stop("'calibration' must be a list from calibrate_afm_mcd() ",
         "containing 'mu_r' and 'Sw'.")
  }
  check_variables(new_data, variables, "new_data",
                  "monitor_afm_mcd(new_data, calibration, variables)")
  if (length(variables) != length(calibration$mu_r)) {
    stop("Number of variables (", length(variables), ") does not match ",
         "calibration mu_r dimension (", length(calibration$mu_r), ").")
  }
  if (!is.null(ucl)) {
    # El error tipico es pasar el objeto de ucl_F_adjusted() en vez de su $UCL.
    if (is.list(ucl) && !is.null(ucl$UCL)) {
      stop("'ucl' must be a single positive number, not the list returned by ",
           "ucl_F_adjusted(). Pass ucl_F_adjusted(calibration, I)$UCL instead.")
    }
    if (!is.numeric(ucl) || length(ucl) != 1 || !is.finite(ucl) || ucl <= 0) {
      stop("'ucl' must be a single positive, finite number (the control ",
           "limit to compare each T-squared against), or NULL to omit the ",
           "out-of-control flag. Typical use: ",
           "ucl_F_adjusted(calibration, I)$UCL.")
    }
  }

  # --- Setup ---
  mu_r <- calibration$mu_r
  # Misma proteccion que hotelling_classical_monitor(): sin ella un Sw singular
  # aborta con el mensaje criptico de solve().
  Sw_inv <- tryCatch(
    solve(calibration$Sw),
    error = function(e) {
      stop("The AFM-weighted covariance matrix Sw is singular and cannot be ",
           "inverted: ", conditionMessage(e),
           "\nThis usually means two or more monitored variables are ",
           "redundant (one is a copy, a sum or a fixed multiple of another), ",
           "or that the Phase 1 batches are too small. Drop the redundant ",
           "variable from 'variables' and calibrate again, or use larger ",
           "Phase 1 batches.")
    }
  )
  batch_id <- new_data[[batch_col]]
  batches <- unique(batch_id)

  # --- Compute T-squared per batch ---
  results <- data.frame(
    Batch = character(),
    I = integer(),
    T2 = numeric(),
    stringsAsFactors = FALSE
  )

  for (batch in batches) {
    subset_batch <- new_data[batch_id == batch, variables]
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

  # --- Optional out-of-control flag ---
  # Estrictamente mayor: un lote justo sobre el limite no se marca.
  if (!is.null(ucl)) {
    results$is_ooc <- results$T2 > ucl
  }

  return(results)
}
