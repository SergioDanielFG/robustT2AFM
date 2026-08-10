#' Calibrate AFM-Weighted MCD Parameters (Phase 1)
#'
#' Estimates robust reference parameters for Phase 1 calibration of an AFM-weighted
#' Hotelling T-squared control chart. For each batch, MCD estimators of location
#' and scatter are computed. Batch covariance matrices are then combined into a
#' global reference matrix using AFM inverse weighting based on first eigenvalues
#' (Abdi, Williams & Valentin, 2013), which gives less influence to batches with higher
#' dispersion. The global reference center is the average of the MCD batch centers.
#'
#' @param data A data frame containing Phase 1 process data. Must contain a column
#'   identifying each batch; see \code{batch_col}.
#' @param variables Character vector with the names of the process variables.
#' @param mcd_alpha Numeric in (0.60, 0.90). Proportion of observations retained
#'   by MCD. Default 0.67 (breakdown point = 0.33), the value used in
#'   Frutos-Galarza et al. (2026).
#' @param verbose Logical. If \code{TRUE}, a message listing the batches that
#'   were actually used for calibration is emitted. Default \code{FALSE}
#'   (silent). Set it to \code{TRUE} when you want to check which batches
#'   survived the minimum-size filter; the same information is always available
#'   afterwards in \code{names(result$weights)}.
#' @param batch_col Character. Name of the column that identifies the batch.
#'   Default \code{"Batch"}. Your export will not necessarily use that name,
#'   and nothing in the method requires it; only this package's own simulated
#'   data does.
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
#'   \item{I_phase1}{Observations per Phase 1 batch, recorded so the UCL can
#'         detect a mismatching I. When the batches are not all the same size
#'         this is the batch size whose m* equals the mean of the per-batch
#'         m*, which is not the mean batch size rounded; see Details.}
#'   \item{batch_sizes}{Named integer vector with the size of every valid
#'         Phase 1 batch, so that downstream functions can check the
#'         equal-size assumption instead of trusting a single number.}
#' }
#'
#' @details
#' AFM comes from \emph{Analyse Factorielle Multiple}; the English literature
#' calls the technique MFA. This package uses AFM throughout.
#'
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
#' \strong{Batches of unequal size.} The calibration itself is well defined
#' whatever the batch sizes: MCD, the weights, Sw and mu_r are computed per
#' batch. The control limit is not: \code{\link{ucl_F_adjusted}} assumes a
#' single common batch size I, through \eqn{m^* = round(I h)}. When the Phase 1
#' batches differ in size the function warns and sets \code{I_phase1} to the
#' batch size whose \eqn{m^*} equals the mean of the per-batch \eqn{m^*_k},
#' that is \code{round(mean(round(I_k * h)) / h)}. This is \strong{not} the
#' same as rounding the mean batch size: the two disagree for roughly one
#' calibration in five.
#'
#' @references
#' Abdi, H., Williams, L. J., & Valentin, D. (2013). Multiple factor analysis:
#' principal component analysis for multitable and multiblock data sets.
#' \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, 5(2),
#' 149-179. \doi{10.1002/wics.1246}
#'
#' Rousseeuw, P. J., & Van Driessen, K. (1999). A fast algorithm for the
#' minimum covariance determinant estimator. Technometrics, 41(3), 212-223.
#'
#' @importFrom robustbase covMcd
#' @export
#'
#' @examples
#' # 30 historical batches, 6 of them carrying outlying observations.
#' # In production this data frame comes from read.csv() or your MES.
#' data(afm_phase1)
#' cal <- calibrate_afm_mcd(afm_phase1, paste0("Var", 1:4))
#'
#' # The reference centre lands near zero, the true process mean, even though
#' # 6 of the 30 calibration batches carry outliers.
#' round(cal$mu_r, 3)               # robust reference center
#' round(cal$Sw, 3)                 # AFM-weighted covariance
#' round(sort(cal$weights), 4)      # smallest weights = most dispersed batches
#'
#' # If your batch identifier is not called "Batch":
#' lots <- afm_phase1
#' names(lots)[names(lots) == "Batch"] <- "Lote"
#' cal2 <- calibrate_afm_mcd(lots, paste0("Var", 1:4), batch_col = "Lote")
#'
#' # To see which batches were actually used, ask for it:
#' cal_v <- calibrate_afm_mcd(afm_phase1, paste0("Var", 1:4), verbose = TRUE)
calibrate_afm_mcd <- function(data, variables, mcd_alpha = 0.67,
                              verbose = FALSE, batch_col = "Batch") {

  # --- Input validation ---
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  check_batch_col(data, batch_col, "data",
                  "calibrate_afm_mcd(data, variables, batch_col = \"<name>\")")
  check_variables(data, variables, "data",
                  "calibrate_afm_mcd(data, variables = c(\"Var1\", \"Var2\"))")
  if (length(variables) < 2) {
    stop("The T-squared chart needs at least 2 process variables; you passed ",
         length(variables), " ('", variables, "'). With a single variable ",
         "there is no covariance structure to weight, and T-squared reduces ",
         "to the square of a standardised mean. Chart that variable on its ",
         "own with a univariate Shewhart X-bar chart, or pass the rest of ",
         "the process variables: ",
         "calibrate_afm_mcd(data, variables = c(\"Var1\", \"Var2\")).",
         call. = FALSE)
  }
  non_numeric <- variables[!vapply(data[variables], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("The following 'variables' are not numeric: ",
         paste(non_numeric, collapse = ", "))
  }
  if (!is.numeric(mcd_alpha) || length(mcd_alpha) != 1 ||
      mcd_alpha < 0.60 || mcd_alpha > 0.90) {
    stop("mcd_alpha must be a single numeric value in [0.60, 0.90]. ",
         "You provided: ", mcd_alpha)
  }
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("'verbose' must be a single logical value (TRUE or FALSE).")
  }

  # --- Setup ---
  batch_id <- data[[batch_col]]
  batches <- unique(batch_id)
  J <- length(variables)

  # --- MCD estimation per batch ---
  mcd_centers <- list()
  mcd_covariances <- list()
  batch_sizes <- integer(0)

  for (batch in batches) {
    subset_batch <- data[batch_id == batch, variables, drop = FALSE]
    if (nrow(subset_batch) <= J) {
      warning("Batch '", batch, "' has too few observations (",
              nrow(subset_batch), " <= ", J, " variables). Skipping.")
      next
    }
    if (any(!is.finite(as.matrix(subset_batch)))) {
      stop("Batch '", batch, "' contains non-finite values (NA/NaN/Inf).")
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
  # El listado de lotes validos solo se emite bajo peticion explicita:
  # la mayoria de los scripts de analisis lo envolvian en suppressMessages().
  if (verbose) {
    message("Valid batches used for calibration (", length(valid_batches),
            "): ", paste(valid_batches, collapse = ", "))
  }

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

  # --- Phase 1 batch size ---
  # Con lotes iguales es ese tamaño. Con lotes desiguales no existe un I
  # exacto: se avisa y se usa el I cuyo m* iguala la media de los m*_k, que
  # es lo que hace coincidir sum_k (m*_k - 1) con los K(m*-1) de la formula.
  # OJO: no es lo mismo que redondear la media de los tamaños; las dos formas
  # difieren en torno a una calibracion de cada cinco.
  batch_sizes <- vapply(batch_sizes, as.integer, integer(1))
  sizes <- as.integer(batch_sizes)
  if (length(unique(sizes)) == 1L) {
    I_phase1 <- sizes[1]
  } else {
    m_target <- round(mean(round(sizes * mcd_alpha)))
    I_phase1 <- as.integer(round(m_target / mcd_alpha))
    modal <- as.integer(names(sort(table(sizes), decreasing = TRUE))[1])
    warning("Phase 1 batches do not all have the same number of observations: ",
            "they range from ", min(sizes), " to ", max(sizes), " (",
            sum(sizes == modal), " of the ", length(sizes),
            " batches have ", modal, "). The control limit formula assumes ",
            "one common batch size, so no exact limit exists for this ",
            "calibration. I_phase1 is set to ", I_phase1, ", the size whose ",
            "m* = round(I * ", mcd_alpha, ") equals ", m_target,
            ", the average of the per-batch m*, which is what matches the ",
            "degrees of freedom of the covariance estimate. To choose it ",
            "yourself: ucl_F_adjusted(calibration, I = <your value>).",
            call. = FALSE)
  }

  return(list(
    mu_r = mu_r,
    Sw = Sw,
    weights = weights,
    mcd_centers = mcd_centers,
    mcd_covariances = mcd_covariances,
    lambda1 = lambda1,
    mcd_alpha = mcd_alpha,
    I_phase1 = I_phase1,
    batch_sizes = batch_sizes
  ))
}
