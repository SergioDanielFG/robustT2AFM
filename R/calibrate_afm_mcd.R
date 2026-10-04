#' Calibrate AFM-Weighted MCD Parameters (Phase 1)
#'
#' Estimates robust reference parameters for Phase 1 calibration of an AFM-weighted
#' Hotelling T-squared control chart. For each batch, reweighted MCD estimators of
#' location and scatter are computed. Batch covariance matrices are then combined
#' into a global reference matrix using AFM inverse weighting based on first
#' eigenvalues (Abdi, Williams & Valentin, 2013), which gives less influence to
#' batches with higher dispersion. The eigenvalues are taken after putting the
#' variables on a common robust scale (MAD), so the weights do not depend on the
#' units of measurement. The global reference center is the reweighted MCD of the
#' batch centers, which resists whole batches shifted in mean.
#'
#' @param data A data frame containing Phase 1 process data. Must contain a column
#'   identifying each batch; see \code{batch_col}.
#' @param variables Character vector with the names of the process variables.
#' @param mcd_alpha Numeric in (0.60, 0.90). Proportion of observations retained
#'   by MCD. Default 0.67 (breakdown point = 0.33), the value used in
#'   Frutos-Galarza et al. (2026).
#' @param scaling Character, \code{"mad"} (default) or \code{"none"}. With
#'   \code{"mad"}, each variable is divided by its median absolute deviation
#'   over all Phase 1 observations before the first eigenvalues are computed,
#'   so the AFM weights do not change when a variable is expressed in other
#'   units. The divisor is common to all batches, so a batch with inflated
#'   dispersion still stands out after scaling. \code{Sw} and \code{mu_r}
#'   are returned in the original units. \code{"none"} reproduces version
#'   0.2.0.
#' @param center Character, \code{"mcd"} (default) or \code{"mean"}. With
#'   \code{"mcd"}, the reference center is the reweighted MCD location of the
#'   K batch centers, which is robust to whole batches shifted in mean and
#'   affine equivariant. \code{"mean"} takes the simple average of the batch
#'   centers and reproduces version 0.2.0.
#' @param center_alpha Numeric in [0.5, 1). Proportion of batch centers kept by
#'   the MCD of the centers when \code{center = "mcd"}. Default 0.5, the value
#'   with the highest breakdown point.
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
#'         this is the smallest batch size whose m* equals the mean of the
#'         per-batch m*, which is not the mean batch size rounded; see
#'         Details.}
#'   \item{batch_sizes}{Named integer vector with the size of every valid
#'         Phase 1 batch, so that downstream functions can check the
#'         equal-size assumption instead of trusting a single number.}
#'   \item{scale}{Divisor applied to each variable before the eigenvalues
#'         (all ones when \code{scaling = "none"}).}
#'   \item{scaling, center, center_alpha}{The options used, recorded so every
#'         calibration states how it was obtained.}
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
#' \strong{Scaling.} With \code{scaling = "mad"}, \eqn{\lambda_{1,k}} is the
#' first eigenvalue of \eqn{D^{-1} S_k D^{-1}}, where \eqn{D} is the diagonal
#' matrix of the MAD of each variable over all Phase 1 observations. Because
#' the MCD is affine equivariant this equals fitting the MCD to the scaled
#' data, and the weights are invariant to a change of units of any variable.
#' Scaling each batch by its own dispersion (for example, using the batch
#' correlation matrices) would remove the variance inflation that the weights
#' rely on to recognise a contaminated batch; that is why the divisor is
#' common to all batches.
#'
#' \strong{Reference center.} With \code{center = "mcd"}, \code{mu_r} is the
#' reweighted MCD location of the K batch centers, computed with
#' \code{robustbase::covMcd(alpha = center_alpha)}. It needs at least twice as
#' many valid batches as variables; with fewer the function stops and asks
#' for \code{center = "mean"}, instead of changing the method silently.
#'
#' Both MCD steps use FastMCD, which draws random subsets. The MCD of the
#' centers runs under a fixed internal seed, so the same Phase 1 always gives
#' the same \code{mu_r}, and the user's random number stream is restored
#' afterwards. The seed only decides where the search starts, not what the
#' optimum is; with 500 starts (the default of \code{covMcd}) different seeds
#' reach the same subset or differ in the third decimal of the center. The
#' per-batch MCD is not seeded internally, as in version 0.2.0: call
#' \code{set.seed()} before the calibration to fix it.
#'
#' \strong{Batches of unequal size.} The calibration itself is well defined
#' whatever the batch sizes: MCD, the weights, Sw and mu_r are computed per
#' batch. The control limit is not: \code{\link{ucl_F_adjusted}} assumes a
#' single common batch size I, through \eqn{m^*}, the size of the MCD subset
#' of a batch of that size (\code{robustbase::h.alpha.n(mcd_alpha, I, J)}).
#' When the Phase 1 batches differ in size the function warns and sets
#' \code{I_phase1} to the smallest batch size whose \eqn{m^*} equals the
#' rounded mean of the per-batch \eqn{m^*_k}; several sizes can share the
#' same \eqn{m^*}, and all of them give the same limit. This is
#' \strong{not} the same as rounding the mean batch size, and the two can
#' disagree.
#'
#' @references
#' Abdi, H., Williams, L. J., & Valentin, D. (2013). Multiple factor analysis:
#' principal component analysis for multitable and multiblock data sets.
#' \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, 5(2),
#' 149-179. \doi{10.1002/wics.1246}
#'
#' Hubert, M., Debruyne, M., & Rousseeuw, P. J. (2018). Minimum covariance
#' determinant and extensions. \emph{Wiley Interdisciplinary Reviews:
#' Computational Statistics}, 10(3), e1421. \doi{10.1002/wics.1421}
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
#'
#' # The calibration of version 0.2.0 (raw weights, average of the centers):
#' cal_020 <- calibrate_afm_mcd(afm_phase1, paste0("Var", 1:4),
#'                              scaling = "none", center = "mean")
calibrate_afm_mcd <- function(data, variables, mcd_alpha = 0.67,
                              scaling = c("mad", "none"),
                              center = c("mcd", "mean"),
                              center_alpha = 0.5,
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
  scaling <- match.arg(scaling)
  center  <- match.arg(center)
  if (!is.numeric(center_alpha) || length(center_alpha) != 1 ||
      is.na(center_alpha) || center_alpha < 0.5 || center_alpha >= 1) {
    stop("center_alpha must be a single numeric value in [0.5, 1). ",
         "You provided: ", center_alpha)
  }

  # --- Setup ---
  batch_id <- data[[batch_col]]
  batches <- unique(batch_id)
  J <- length(variables)

  # --- Common robust scale of each variable (used only for the weights) ---
  if (scaling == "mad") {
    scale_j <- apply(as.matrix(data[, variables, drop = FALSE]), 2,
                     stats::mad, na.rm = TRUE)
    zero <- variables[!is.finite(scale_j) | scale_j <= 0]
    if (length(zero) > 0) {
      stop("The following variables have zero median absolute deviation in ",
           "Phase 1, so they cannot be put on a common scale: ",
           paste(zero, collapse = ", "), ". A variable that barely moves ",
           "carries no information for the chart; remove it, or use ",
           "scaling = \"none\".", call. = FALSE)
    }
  } else {
    scale_j <- rep(1, J)
  }
  names(scale_j) <- variables

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

  # --- First eigenvalue per batch, on the common robust scale ---
  # D^-1 S_k D^-1 es la covarianza MCD del lote con las variables divididas
  # por su MAD; con scaling = "none" D es la identidad (version 0.2.0).
  D_inv <- diag(1 / scale_j, nrow = J)
  lambda1 <- sapply(mcd_covariances, function(S) {
    max(eigen(D_inv %*% S %*% D_inv, symmetric = TRUE,
              only.values = TRUE)$values)
  })

  # --- AFM inverse weights: w_k = (1/lambda1_k) / sum(1/lambda1) ---
  # Batches with higher dispersion receive LESS weight (AFM standard formulation)
  inv_lambda1 <- 1 / lambda1
  weights <- inv_lambda1 / sum(inv_lambda1)
  names(weights) <- valid_batches

  # --- AFM-weighted reference covariance matrix ---
  Sw <- Reduce("+", Map(function(w, S) w * S, weights, mcd_covariances))

  # --- Reference center ---
  centers_matrix <- do.call(rbind, mcd_centers)
  if (center == "mcd") {
    if (nrow(centers_matrix) < 2 * J) {
      stop("center = \"mcd\" needs at least twice as many valid batches as ",
           "variables (", 2 * J, " here): there are ", nrow(centers_matrix),
           " batches and ", J, " variables. With so few batches the MCD of ",
           "the batch centers cannot be computed reliably. Add Phase 1 ",
           "batches, or use center = \"mean\".", call. = FALSE)
    }
    # MCD reponderado de los K centros (FastMCD, como dentro de cada lote).
    # Se fija una semilla local para que la misma Fase 1 de siempre el mismo
    # centro, y se restaura el estado del generador del usuario al salir, de
    # modo que sus simulaciones no se ven alteradas.
    mu_r <- with_local_seed(CENTER_SEED, {
      robustbase::covMcd(centers_matrix, alpha = center_alpha)$center
    })
  } else {
    mu_r <- colMeans(centers_matrix)             # version 0.2.0
  }
  names(mu_r) <- variables

  # --- Phase 1 batch size ---
  # Con lotes iguales es ese tamaño. Con lotes desiguales no existe un I
  # exacto: se avisa y se usa el menor I cuyo m* (tamano del subconjunto del
  # MCD) iguala la media redondeada de los m*_k, que es lo que hace coincidir
  # sum_k (m*_k - 1) con los K(m*-1) de la formula. Varios I comparten el
  # mismo m* y dan el mismo limite. OJO: no es lo mismo que redondear la
  # media de los tamaños.
  batch_sizes <- vapply(batch_sizes, as.integer, integer(1))
  sizes <- as.integer(batch_sizes)
  if (length(unique(sizes)) == 1L) {
    I_phase1 <- sizes[1]
  } else {
    m_of <- function(n) robustbase::h.alpha.n(mcd_alpha, n, J)
    m_target <- round(mean(vapply(sizes, m_of, numeric(1))))
    cand <- seq.int(J + 1L, 2L * max(sizes))
    I_phase1 <- as.integer(cand[which(vapply(cand, m_of, numeric(1)) >= m_target)[1]])
    modal <- as.integer(names(sort(table(sizes), decreasing = TRUE))[1])
    warning("Phase 1 batches do not all have the same number of observations: ",
            "they range from ", min(sizes), " to ", max(sizes), " (",
            sum(sizes == modal), " of the ", length(sizes),
            " batches have ", modal, "). The control limit formula assumes ",
            "one common batch size, so no exact limit exists for this ",
            "calibration. I_phase1 is set to ", I_phase1, ", the smallest ",
            "size whose m* (the MCD subset size) equals ", m_target,
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
    batch_sizes = batch_sizes,
    scale = scale_j,
    scaling = scaling,
    center = center,
    center_alpha = center_alpha
  ))
}
