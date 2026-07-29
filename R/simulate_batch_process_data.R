#' Simulate Batch Process Data with Two Types of Contamination
#'
#' Generates synthetic batch process data with configurable contamination:
#' outliers within batches (to test MCD robustness), and entire shifted
#' batches (to test AFM weighting robustness). Both Phase 1 (calibration)
#' and Phase 2 (monitoring) can be configured.
#'
#' @param K1 Integer. Number of Phase 1 batches. Default 30.
#' @param K2 Integer. Number of Phase 2 batches. Default 20.
#' @param I Integer. Observations per batch. Default 20 (value used in
#'   Ruiz-Barzola et al. 2026).
#' @param J Integer. Number of process variables. Default 4.
#' @param mu Numeric vector of length J. Default zero vector.
#' @param Sigma Covariance matrix J x J. If NULL, equicorrelation with correlation rho.
#' @param rho Numeric in (-1, 1). Base correlation. Default 0.6.
#' @param outlier_batches_F1 Integer in from 0 to K1. Number of Phase 1 batches
#'   containing outlier observations. Default 0. Tests MCD robustness.
#' @param outlier_rate Numeric in (0, 0.5). Fraction of observations within
#'   an affected batch that are outliers. Default 0.20.
#' @param outlier_shift Numeric. Shift magnitude (in sigmas) applied to
#'   outlier observations. Default 4.
#' @param prop_contam_F1 Numeric in (0, 1). Proportion of Phase 1 batches
#'   entirely shifted. Default 0. Tests AFM weighting robustness.
#' @param shift_contam Numeric. Shift magnitude (in sigmas) for shifted
#'   batches. Default 3.
#' @param prop_ooc_F2 Numeric in (0, 1). Proportion of Phase 2 batches OOC.
#'   Default 0.
#' @param shift_ooc Numeric. Shift magnitude (in sigmas) for OOC batches.
#'   Default 2.
#' @param seed Optional integer for reproducibility.
#'
#' @return Data frame with columns Batch, Phase, Status, ContaminationType,
#'   and Var1..VarJ. ContaminationType describes each batch: "Clean",
#'   "Outliers", "Shifted", or "OOC".
#'
#' @details
#' Two types of Phase 1 contamination are supported:
#'
#' Type 1 - Outliers within batches (tests MCD): A fraction \code{outlier_rate}
#' of observations within \code{outlier_batches_F1} Phase 1 batches are shifted
#' by \code{outlier_shift} sigmas. MCD should identify and trim these observations.
#'
#' Type 2 - Entire shifted batches (tests AFM): A proportion \code{prop_contam_F1}
#' of Phase 1 batches are entirely shifted by \code{shift_contam} sigmas. AFM
#' weighting should give these batches lower weight in Sw.
#'
#' Both types can coexist. Phase 2 OOC batches are independent of Phase 1
#' contamination settings.
#'
#' @importFrom MASS mvrnorm
#' @export
#'
#' @examples
#' # Realistic Phase 1 with within-batch outliers (2 batches) and fully-
#' # shifted batches (~7%), plus Phase 2 with 30% OOC. Package defaults
#' # (K1=30, K2=20, I=20, J=4) match Ruiz-Barzola et al. (2026).
#' sim <- simulate_batch_process(
#'   outlier_batches_F1 = 2, outlier_rate = 0.20, outlier_shift = 4,
#'   prop_contam_F1 = 0.07, shift_contam = 3,
#'   prop_ooc_F2 = 0.30,   shift_ooc = 2,
#'   seed = 20260417
#' )
#' table(sim$Phase, sim$ContaminationType)   # composition by phase
#' head(sim, 3)
simulate_batch_process <- function(K1 = 30, K2 = 20, I = 20, J = 4,
                                   mu = NULL, Sigma = NULL, rho = 0.6,
                                   outlier_batches_F1 = 0,
                                   outlier_rate = 0.20,
                                   outlier_shift = 4,
                                   prop_contam_F1 = 0, shift_contam = 3,
                                   prop_ooc_F2 = 0, shift_ooc = 2,
                                   seed = NULL) {

  # --- Input validation ---
  if (K1 < 2 || K1 != round(K1)) stop("'K1' must be integer >= 2.")
  if (K2 < 0 || K2 != round(K2)) stop("'K2' must be non-negative integer.")
  if (I < 2 || I != round(I)) stop("'I' must be integer >= 2.")
  if (J < 2 || J != round(J)) stop("'J' must be integer >= 2.")
  if (outlier_batches_F1 < 0 || outlier_batches_F1 > K1) {
    stop("'outlier_batches_F1' must be in [0, K1].")
  }
  if (outlier_rate < 0 || outlier_rate > 0.5) {
    stop("'outlier_rate' must be in [0, 0.5].")
  }
  if (prop_contam_F1 < 0 || prop_contam_F1 >= 1) {
    stop("'prop_contam_F1' must be in [0, 1).")
  }
  if (prop_ooc_F2 < 0 || prop_ooc_F2 >= 1) {
    stop("'prop_ooc_F2' must be in [0, 1).")
  }

  if (!is.null(seed)) set.seed(seed)

  # --- Default mu and Sigma ---
  if (is.null(mu)) mu <- rep(0, J)
  if (length(mu) != J) stop("'mu' must have length J.")

  if (is.null(Sigma)) {
    Sigma <- matrix(rho, J, J)
    diag(Sigma) <- 1
  }

  var_names <- paste0("Var", seq_len(J))
  sigma_vec <- sqrt(diag(Sigma))

  # --- Determine which batches are affected ---
  n_shifted_F1 <- round(K1 * prop_contam_F1)
  shifted_batches_F1 <- if (n_shifted_F1 > 0) {
    sample(seq_len(K1), n_shifted_F1)
  } else integer(0)

  # Outlier batches: from those NOT already shifted
  available_for_outliers <- setdiff(seq_len(K1), shifted_batches_F1)
  outlier_batches_actual <- if (outlier_batches_F1 > 0 && length(available_for_outliers) > 0) {
    sample(available_for_outliers, min(outlier_batches_F1, length(available_for_outliers)))
  } else integer(0)

  n_ooc_F2 <- round(K2 * prop_ooc_F2)
  ooc_batches_F2 <- if (n_ooc_F2 > 0) sample(seq_len(K2), n_ooc_F2) else integer(0)

  # --- Helper: generate a clean batch ---
  gen_batch <- function(I_local, mu_vec, Sigma_mat) {
    X <- MASS::mvrnorm(n = I_local, mu = mu_vec, Sigma = Sigma_mat)
    colnames(X) <- var_names
    as.data.frame(X)
  }

  # --- Helper: generate a batch with outliers within it ---
  gen_batch_with_outliers <- function(I_local, mu_vec, Sigma_mat,
                                      outlier_rate_local, outlier_shift_local) {
    X <- MASS::mvrnorm(n = I_local, mu = mu_vec, Sigma = Sigma_mat)
    n_outliers <- max(1, round(I_local * outlier_rate_local))
    outlier_indices <- sample(seq_len(I_local), n_outliers)

    X[outlier_indices, ] <- X[outlier_indices, ] +
      matrix(outlier_shift_local * sigma_vec,
             nrow = n_outliers, ncol = J, byrow = TRUE)

    colnames(X) <- var_names
    as.data.frame(X)
  }

  # --- Phase 1 ---
  phase1_list <- lapply(seq_len(K1), function(k) {
    is_shifted <- k %in% shifted_batches_F1
    has_outliers <- k %in% outlier_batches_actual

    if (is_shifted) {
      mu_k <- mu + shift_contam * sigma_vec
      X <- gen_batch(I, mu_k, Sigma)
      contam_type <- "Shifted"
      status <- "Out of Control"
    } else if (has_outliers) {
      X <- gen_batch_with_outliers(I, mu, Sigma, outlier_rate, outlier_shift)
      contam_type <- "Outliers"
      status <- "Out of Control"
    } else {
      X <- gen_batch(I, mu, Sigma)
      contam_type <- "Clean"
      status <- "Under Control"
    }

    data.frame(
      Batch = paste0("F1_B", sprintf("%02d", k)),
      Phase = "Phase 1",
      Status = status,
      ContaminationType = contam_type,
      X,
      stringsAsFactors = FALSE
    )
  })

  # --- Phase 2 ---
  phase2_list <- if (K2 > 0) {
    lapply(seq_len(K2), function(k) {
      is_ooc <- k %in% ooc_batches_F2
      mu_k <- if (is_ooc) mu + shift_ooc * sigma_vec else mu
      X <- gen_batch(I, mu_k, Sigma)
      data.frame(
        Batch = paste0("F2_B", sprintf("%02d", k)),
        Phase = "Phase 2",
        Status = if (is_ooc) "Out of Control" else "Under Control",
        ContaminationType = if (is_ooc) "OOC" else "Clean",
        X,
        stringsAsFactors = FALSE
      )
    })
  } else {
    list()
  }

  # --- Combine ---
  sim_data <- do.call(rbind, c(phase1_list, phase2_list))
  sim_data$Batch <- factor(sim_data$Batch)
  sim_data$Phase <- factor(sim_data$Phase, levels = c("Phase 1", "Phase 2"))
  sim_data$Status <- factor(sim_data$Status, levels = c("Under Control", "Out of Control"))
  sim_data$ContaminationType <- factor(sim_data$ContaminationType,
                                       levels = c("Clean", "Outliers", "Shifted", "OOC"))

  return(sim_data)
}
