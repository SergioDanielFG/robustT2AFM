#' Calibrate AFM-Weighted MCD Parameters (Phase 1)
#'
#' Estimates robust reference parameters for Phase 1 calibration of an AFM-weighted
#' Hotelling T-squared control chart. For each batch, MCD estimators of location
#' and scatter are computed. Batch covariance matrices are then combined into a
#' global reference matrix using AFM inverse weighting based on first eigenvalues
#' (Escofier & Pages, 1994), which gives less influence to batches with higher
#' dispersion. The global reference center is the average of the MCD batch centers.
#'
#' @param data A data frame containing Phase 1 process data. Must contain a column
#'   named 'Batch' identifying each batch.
#' @param variables Character vector with the names of the process variables.
#' @param mcd_alpha Numeric in (0.60, 0.90). Proportion of observations retained
#'   by MCD. Default 0.67 (breakdown point = 0.33), the value used in
#'   Ruiz-Barzola et al. (2026).
#'
#' @return A list containing:
#' \describe{
#'   \item{mu_r}{Reference center vector (length J).}
#'   \item{Sw}{AFM-weighted reference covariance matrix (J x J).}
#'   \item{weights}{AFM inverse weights per batch (named vector, sum to 1).
#'         Batches with the smallest weight are the ones AFM identified as
#'         most anomalous during Phase 1.}
#'   \item{mcd_centers}{List of MCD centers per valid batch.}
#'   \item{mcd_covariances}{List of MCD covariance matrices per valid batch.}
#'   \item{lambda1}{Named vector of first eigenvalues per batch.}
#'   \item{mcd_alpha}{MCD alpha parameter used (for reference).}
#'   \item{I_phase1}{Observations per Phase 1 batch (size of the first valid
#'         batch), recorded so the UCL can detect a mismatching I.}
#' }
#'
#' @details
#' The AFM inverse weighting is computed as:
#' \deqn{w_k = (1/\lambda_{1,k}) / \sum_i (1/\lambda_{1,i})}
#' where \eqn{\lambda_{1,k}} is the first eigenvalue of the MCD covariance
#' matrix of batch k. Batches with higher dispersion (larger first eigenvalue)
#' receive less weight, following the AFM philosophy.
#'
#' The reference covariance matrix is:
#' \deqn{S_w = \sum_k w_k S_k}
#' where \eqn{S_k} is the MCD covariance of batch k.
#'
#' @references
#' Escofier, B., & Pages, J. (1994). Multiple factor analysis (AFMULT package).
#' Computational Statistics & Data Analysis, 18(1), 121-140.
#'
#' Rousseeuw, P. J., & Van Driessen, K. (1999). A fast algorithm for the
#' minimum covariance determinant estimator. Technometrics, 41(3), 212-223.
#'
#' @importFrom robustbase covMcd
#' @export
#'
#' @examples
#' # Typical quality-engineer workflow:
#' # 1) Assemble Phase 1 historical batches into a data.frame with a
#' #    "Batch" column and the J process variables. In practice this comes
#' #    from read.csv() or your MES. Here we simulate a realistic Phase 1
#' #    that already contains 2 batches with within-batch outliers and
#' #    2 fully-shifted batches, i.e. contamination the calibration must
#' #    absorb without polluting the reference.
#' sim <- simulate_batch_process(
#'   K1 = 30, K2 = 0, I = 20, J = 4,
#'   outlier_batches_F1 = 2, outlier_rate = 0.20, outlier_shift = 4,
#'   prop_contam_F1 = 0.07, shift_contam = 3,
#'   seed = 20260417
#' )
#' phase1 <- subset(sim, Phase == "Phase 1")
#' vars   <- paste0("Var", 1:4)
#'
#' # 2) Calibrate the AFM-weighted MCD reference.
#' cal <- calibrate_afm_mcd(phase1, vars)
#'
#' # 3) Inspect the outputs a quality engineer cares about:
#' round(cal$mu_r, 3)               # robust reference center
#' round(cal$Sw, 3)                 # AFM-weighted covariance
#' round(sort(cal$weights), 4)      # smallest weights = most anomalous batches
calibrate_afm_mcd <- function(data, variables, mcd_alpha = 0.67) {

  # --- Input validation ---
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  if (!"Batch" %in% colnames(data)) {
    stop("'data' must contain a column named 'Batch'.")
  }
  if (!all(variables %in% colnames(data))) {
    missing_vars <- setdiff(variables, colnames(data))
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }
  if (!is.numeric(mcd_alpha) || length(mcd_alpha) != 1 ||
      mcd_alpha < 0.60 || mcd_alpha > 0.90) {
    stop("mcd_alpha must be a single numeric value in [0.60, 0.90]. ",
         "You provided: ", mcd_alpha)
  }

  # --- Setup ---
  batches <- unique(data$Batch)
  J <- length(variables)

  # --- MCD estimation per batch ---
  mcd_centers <- list()
  mcd_covariances <- list()
  batch_sizes <- integer(0)

  for (batch in batches) {
    subset_batch <- data[data$Batch == batch, variables]
    if (nrow(subset_batch) <= J) {
      warning("Batch '", batch, "' has too few observations (",
              nrow(subset_batch), " <= ", J, " variables). Skipping.")
      next
    }
    mcd_est <- robustbase::covMcd(subset_batch, alpha = mcd_alpha)
    mcd_centers[[as.character(batch)]] <- mcd_est$center
    mcd_covariances[[as.character(batch)]] <- mcd_est$cov
    batch_sizes[[as.character(batch)]] <- nrow(subset_batch)
  }

  valid_batches <- names(mcd_centers)
  if (length(valid_batches) < 2) {
    stop("At least 2 valid batches are required after MCD estimation. ",
         "Only ", length(valid_batches), " valid batches found.")
  }
  message("Valid batches used for calibration (", length(valid_batches),
          "): ", paste(valid_batches, collapse = ", "))

  # --- First eigenvalue per batch ---
  lambda1 <- sapply(mcd_covariances, function(S) {
    max(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
  })

  # --- AFM inverse weights: w_k = (1/lambda1_k) / sum(1/lambda1) ---
  # Batches with higher dispersion receive LESS weight (AFM standard formulation)
  inv_lambda1 <- 1 / lambda1
  weights <- inv_lambda1 / sum(inv_lambda1)
  names(weights) <- valid_batches

  # --- AFM-weighted reference covariance matrix ---
  Sw <- Reduce("+", Map(function(w, S) w * S, weights, mcd_covariances))

  # --- Reference center: mean of MCD centers ---
  centers_matrix <- do.call(rbind, mcd_centers)
  mu_r <- colMeans(centers_matrix)

  # --- Phase 1 batch size (first valid batch) ---
  # Recorded so the UCL can warn if a mismatching I is passed later.
  I_phase1 <- as.integer(batch_sizes[[1]])

  return(list(
    mu_r = mu_r,
    Sw = Sw,
    weights = weights,
    mcd_centers = mcd_centers,
    mcd_covariances = mcd_covariances,
    lambda1 = lambda1,
    mcd_alpha = mcd_alpha,
    I_phase1 = I_phase1
  ))
}
