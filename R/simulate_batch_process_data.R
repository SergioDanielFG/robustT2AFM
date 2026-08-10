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
#'   Frutos-Galarza et al. 2026).
#' @param J Integer. Number of process variables. Default 4.
#' @param mu Numeric vector of length J. Default zero vector.
#' @param Sigma Covariance matrix J x J. If NULL, equicorrelation with correlation rho.
#' @param rho Numeric in (-1, 1). Base correlation. Default 0.6.
#' @param outlier_batches_F1 Integer from 0 to K1. Number of Phase 1 batches
#'   containing outlier observations. Default 0. Tests MCD robustness.
#' @param outlier_rate Numeric in \[0, 0.5\]. Fraction of observations within
#'   an affected batch that are outliers. Default 0.20.
#' @param outlier_shift Numeric. Shift magnitude (in sigmas) applied to
#'   outlier observations. Default 4.
#' @param prop_contam_F1 Numeric in \[0, 1). Proportion of Phase 1 batches
#'   entirely shifted. Default 0. Tests AFM weighting robustness.
#' @param shift_contam Numeric. Shift magnitude (in sigmas) for shifted
#'   batches. Default 3.
#' @param prop_ooc_F2 Numeric in \[0, 1). Proportion of Phase 2 batches OOC.
#'   Default 0.
#' @param shift_ooc Numeric. Shift magnitude (in sigmas) for OOC batches.
#'   Default 2.
#' @param seed Optional integer for reproducibility.
#'
#' @return Data frame with columns Batch, Phase, Status, ContaminationType,
#'   and Var1..VarJ. \code{Status} is "Faulty" or "Fault-free" and
#'   \code{ContaminationType} says how a batch was spoiled: "Clean",
#'   "Outliers" or "Shifted". A wholly shifted batch is labelled "Shifted" in
#'   both phases, because it is the same mechanism: the mean is displaced and
#'   the covariance left alone. What differs is the role, contamination to be
#'   absorbed in Phase 1 and a signal to be detected in Phase 2, and the
#'   \code{Phase} column already says which.
#'
#'   Both columns are \strong{ground truth}, available only because the data
#'   was generated on purpose, and they use the vocabulary of truth rather
#'   than the vocabulary of the chart's verdict: a batch is "Faulty", it is
#'   not "out of control", because being out of control is what a chart
#'   decides and not what a batch is. Real process data carries neither
#'   column. Their intended use is to supply the \code{faulty} argument of
#'   \code{\link{plot_method_comparison}} and to score detections in a
#'   simulation study.
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
#' # (K1=30, K2=20, I=20, J=4) match Frutos-Galarza et al. (2026).
#' sim <- simulate_batch_process(
#'   outlier_batches_F1 = 2, outlier_rate = 0.20, outlier_shift = 4,
#'   prop_contam_F1 = 0.07, shift_contam = 3,
#'   prop_ooc_F2 = 0.30,   shift_ooc = 2,
#'   seed = 20260731
#' )
#' table(sim$Phase, sim$ContaminationType)   # composition by phase
#'
#' # Both contamination types in one data set: apply the method to it.
#' p1 <- sim[sim$Phase == "Phase 1", ]
#' p2 <- sim[sim$Phase == "Phase 2", ]
#' summary(run_afm_mcd(p1, p2, paste0("Var", 1:4), plot = FALSE))
simulate_batch_process <- function(K1 = 30, K2 = 20, I = 20, J = 4,
                                   mu = NULL, Sigma = NULL, rho = 0.6,
                                   outlier_batches_F1 = 0,
                                   outlier_rate = 0.20,
                                   outlier_shift = 4,
                                   prop_contam_F1 = 0, shift_contam = 3,
                                   prop_ooc_F2 = 0, shift_ooc = 2,
                                   seed = NULL) {

  # --- Input validation ---
  if (K1 < 2 || K1 != round(K1)) {
    stop("'K1' must be a whole number >= 2 (Phase 1 batches). You provided: ",
         K1, ".", call. = FALSE)
  }
  if (K2 < 0 || K2 != round(K2)) {
    stop("'K2' must be a whole number >= 0 (Phase 2 batches; 0 generates no ",
         "Phase 2). You provided: ", K2, ".", call. = FALSE)
  }
  if (I < 2 || I != round(I)) {
    stop("'I' must be a whole number >= 2 (observations per batch). You ",
         "provided: ", I, ".", call. = FALSE)
  }
  if (J < 2 || J != round(J)) {
    stop("'J' must be a whole number >= 2 (process variables). You provided: ",
         J, ".", call. = FALSE)
  }
  if (outlier_batches_F1 < 0 || outlier_batches_F1 > K1) {
    stop("'outlier_batches_F1' is the number of Phase 1 batches carrying ",
         "outliers, so it must be between 0 and K1 = ", K1,
         ". You provided: ", outlier_batches_F1, ".", call. = FALSE)
  }
  if (outlier_rate < 0 || outlier_rate > 0.5) {
    stop("'outlier_rate' is the fraction of observations spoiled inside an ",
         "affected batch, so it must be between 0 and 0.5. You provided: ",
         outlier_rate, ".", call. = FALSE)
  }
  if (prop_contam_F1 < 0 || prop_contam_F1 >= 1) {
    stop("'prop_contam_F1' is the proportion of wholly shifted Phase 1 ",
         "batches, so it must be in [0, 1). You provided: ", prop_contam_F1,
         ".", call. = FALSE)
  }
  if (prop_ooc_F2 < 0 || prop_ooc_F2 >= 1) {
    stop("'prop_ooc_F2' is the proportion of off-target Phase 2 batches, so ",
         "it must be in [0, 1); 1 is excluded, use 0.95 for 19 of 20. You ",
         "provided: ", prop_ooc_F2, ".", call. = FALSE)
  }

  if (!is.null(rho) && (!is.numeric(rho) || length(rho) != 1 ||
                        !is.finite(rho) || rho <= -1 || rho >= 1)) {
    stop("'rho' is the base correlation between variables, so it must be a ",
         "single number in (-1, 1). You provided: ", rho, ".", call. = FALSE)
  }
  if (!is.null(Sigma)) {
    if (!is.matrix(Sigma) || nrow(Sigma) != J || ncol(Sigma) != J) {
      stop("'Sigma' must be a J x J matrix, with J = ", J, ".", call. = FALSE)
    }
    if (!isTRUE(all.equal(Sigma, t(Sigma)))) {
      stop("'Sigma' must be symmetric.", call. = FALSE)
    }
    if (any(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
      stop("'Sigma' must be positive definite: all its eigenvalues must be ",
           "greater than zero.", call. = FALSE)
    }
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 ||
                         seed != round(seed))) {
    stop("'seed' must be a single whole number, or NULL.", call. = FALSE)
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
  # sample(x, n) muestrea de 1:x cuando x tiene longitud 1, no de x. Con
  # prop_contam_F1 alto puede quedar un solo lote disponible, y entonces se
  # contaminaria un lote distinto del elegido, sin aviso.
  outlier_batches_actual <- if (outlier_batches_F1 > 0 && length(available_for_outliers) > 0) {
    n_take <- min(outlier_batches_F1, length(available_for_outliers))
    available_for_outliers[sample.int(length(available_for_outliers), n_take)]
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
      status <- "Faulty"
    } else if (has_outliers) {
      X <- gen_batch_with_outliers(I, mu, Sigma, outlier_rate, outlier_shift)
      contam_type <- "Outliers"
      status <- "Faulty"
    } else {
      X <- gen_batch(I, mu, Sigma)
      contam_type <- "Clean"
      status <- "Fault-free"
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
        Status = if (is_ooc) "Faulty" else "Fault-free",
        # "Shifted", igual que en Fase 1: es la misma linea de codigo, mismo
        # desplazamiento de la media y misma covarianza. Lo que cambia es el
        # papel, absorber frente a detectar, y eso ya lo dice la columna Phase.
        ContaminationType = if (is_ooc) "Shifted" else "Clean",
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
  # Status es VERDAD de campo, no veredicto de la carta, asi que usa el
  # vocabulario de la verdad: Faulty / Fault-free. "Out of Control" era el
  # vocabulario del veredicto y confundia las dos ideas justo en la columna
  # que sirve para construir el argumento 'faulty'.
  sim_data$Status <- factor(sim_data$Status,
                            levels = c("Fault-free", "Faulty"))
  sim_data$ContaminationType <- factor(sim_data$ContaminationType,
                                       levels = c("Clean", "Outliers", "Shifted"))

  return(sim_data)
}
