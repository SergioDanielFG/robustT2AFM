#' Simulated Upper Control Limit for AFM-MCD (Phase 2)
#'
#' Computes the Upper Control Limit for the Phase 2 Hotelling T-squared
#' statistic by simulation, as an alternative to the analytic limit of
#' \code{\link{ucl_F_adjusted}}. The limit is the 1 - alpha quantile of the
#' T-squared of in-control batches when the whole method, calibration
#' included, is applied to synthetic Phase 1 samples generated from the
#' robust estimates of the actual calibration (a parametric bootstrap). It
#' therefore carries the sampling variability of the per-batch MCD, of the
#' AFM weights and of the robust center, which the F approximation ignores.
#'
#' @param calibration A list returned by \code{\link{calibrate_afm_mcd}}.
#' @param I Integer. Number of observations per batch in Phase 1 (e.g., 20).
#' @param alpha Numeric in (0, 1). Nominal false alarm rate. Default 0.001
#'   (in-control ARL0 = 1000), as in \code{\link{ucl_F_adjusted}}.
#' @param B Integer. Number of synthetic Phase 1 samples. Default 500.
#' @param n_new Integer. Number of new in-control batches evaluated against
#'   each synthetic calibration. Default 200. The limit is the quantile of
#'   the \code{B * n_new} simulated T-squared values; with the defaults and
#'   alpha = 0.001 about 100 of them lie beyond it.
#' @param seed Integer or \code{NULL}. With a number, the simulation runs
#'   under that seed and the user's random number stream is restored
#'   afterwards, so the same call always gives the same limit. With
#'   \code{NULL} (default) the current stream is used and advanced.
#'
#' @return A list containing:
#' \describe{
#'   \item{UCL}{Upper Control Limit value.}
#'   \item{method}{Character string: "Simulated (parametric bootstrap)".}
#'   \item{parameters}{Named list with J, K, I, alpha, B, n_new, n_T2 (the
#'     number of simulated T-squared values), seed and the calibration
#'     options used for the synthetic samples (mcd_alpha, scaling, center,
#'     center_alpha).}
#' }
#'
#' @details
#' The procedure repeats B times:
#' \enumerate{
#'   \item generate a synthetic Phase 1 of K batches of I observations from
#'         a multivariate normal distribution with mean \code{mu_r} and
#'         covariance \code{Sw} of the actual calibration;
#'   \item calibrate it with \code{\link{calibrate_afm_mcd}}, with the same
#'         options as the actual calibration;
#'   \item generate \code{n_new} in-control batch means from the normal
#'         distribution with mean \code{mu_r} and covariance \code{Sw / I},
#'         and compute their T-squared against the synthetic calibration.
#' }
#' The limit is the 1 - alpha quantile (type 8) of all the T-squared values.
#'
#' The robust estimates play the role of the in-control distribution, which
#' is unknown. This is what makes the limit usable with real data, and also
#' what sets its reach: it assumes in-control normality, and it cannot
#' remove what contamination left in \code{mu_r} and \code{Sw}. In the
#' simulation study of Frutos-Galarza et al. (2026) it brought the in-control
#' ARL0 close to the nominal 1000 when the Phase 1 sample was clean, also for
#' K = 100 batches, where the analytic limit is conservative; with a
#' contaminated Phase 1 it gave the same ARL0 as the analytic limit. Between
#' seeds, the limit varied with a standard deviation of about 0.2 for
#' K = 30, J = 4 and the default B and n_new.
#'
#' The computing time is that of B calibrations: about two minutes with the
#' defaults for K = 30 batches of 20 observations, and proportionally more
#' for larger K. For a quick look use a smaller B, but the limit then rests
#' on fewer simulated values beyond it; the function warns when there are
#' fewer than 10.
#'
#' @seealso \code{\link{ucl_F_adjusted}} for the analytic limit.
#'
#' @importFrom stats quantile
#' @export
#'
#' @examples
#' data(afm_phase1)
#' vars <- paste0("Var", 1:4)
#' cal  <- calibrate_afm_mcd(afm_phase1, vars)
#'
#' # Analytic limit, instantaneous:
#' ucl_F_adjusted(cal, I = 20)$UCL
#'
#' \donttest{
#' # Simulated limit. A small B keeps the example short; use the default
#' # B = 500 in practice.
#' us <- ucl_simulated(cal, I = 20, alpha = 0.01, B = 20, seed = 1)
#' us$UCL
#' }
ucl_simulated <- function(calibration, I, alpha = 0.001, B = 500,
                          n_new = 200, seed = NULL) {

  # --- Input validation ---
  if (!is.list(calibration) ||
      !all(c("mu_r", "Sw", "mcd_alpha", "weights") %in% names(calibration))) {
    stop("'calibration' must be a list from calibrate_afm_mcd().")
  }
  if (!is.numeric(I) || length(I) != 1 || I < 2 || I != round(I)) {
    stop("'I' must be a single positive integer >= 2.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value in (0, 1).")
  }
  if (!is.numeric(B) || length(B) != 1 || B < 1 || B != round(B)) {
    stop("'B' is the number of synthetic Phase 1 samples, so it must be a ",
         "single whole number >= 1. You provided: ", B, ".", call. = FALSE)
  }
  if (!is.numeric(n_new) || length(n_new) != 1 || n_new < 1 ||
      n_new != round(n_new)) {
    stop("'n_new' is the number of new in-control batches per synthetic ",
         "calibration, so it must be a single whole number >= 1. You ",
         "provided: ", n_new, ".", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 ||
                         seed != round(seed))) {
    stop("'seed' must be a single whole number, or NULL.", call. = FALSE)
  }

  J <- length(calibration$mu_r)        # number of variables
  K <- length(calibration$weights)     # number of Phase 1 batches
  if (I <= J) {
    stop("'I' (", I, ") must be larger than the number of variables (", J,
         "): the MCD of each synthetic batch needs more observations than ",
         "variables.", call. = FALSE)
  }

  # --- Defensive checks, as in ucl_F_adjusted() ---
  sizes <- calibration$batch_sizes
  if (!is.null(sizes) && length(unique(as.integer(sizes))) > 1L) {
    warning("Phase 1 batch sizes were not equal (", min(sizes), " to ",
            max(sizes), "); the synthetic batches all have I = ", I,
            " observations, so this limit is an approximation.",
            call. = FALSE)
  } else if (!is.null(calibration$I_phase1) && calibration$I_phase1 != I) {
    warning("'I' (", I, ") does not match the Phase 1 batch size recorded ",
            "during calibration (I_phase1 = ", calibration$I_phase1, "). ",
            "The synthetic batches will have I = ", I, " observations.")
  }
  # Con pocos valores simulados por encima del cuantil, el limite depende de
  # unos pocos lotes extremos y cambia mucho de una semilla a otra.
  n_beyond <- B * n_new * alpha
  if (n_beyond < 10) {
    warning("Only about ", round(n_beyond, 1), " simulated T2 values lie ",
            "beyond the limit (B * n_new * alpha), so it is imprecise. ",
            "Increase B or n_new.", call. = FALSE)
  }

  # --- Options of the actual calibration ---
  # Una calibracion de la version 0.2.0 no guarda las opciones: son las que
  # entonces eran las unicas (sin escala, centro = media de los centros).
  opt <- function(nm, default) {
    if (is.null(calibration[[nm]])) default else calibration[[nm]]
  }
  scaling      <- opt("scaling", "none")
  center       <- opt("center", "mean")
  center_alpha <- opt("center_alpha", 0.5)
  mu   <- calibration$mu_r
  Sw   <- calibration$Sw
  vars <- names(mu)
  if (is.null(vars)) vars <- paste0("Var", seq_len(J))

  # --- The simulation ---
  # Cada muestra sintetica pasa por la misma calibracion que la real, de modo
  # que el limite recoge la variabilidad del MCD por lote, de los pesos y del
  # centro robusto. Solo las medias de lote entran en T2, y la media de I
  # observaciones normales se genera directamente como N(mu, Sw / I).
  run <- function() {
    t2 <- numeric(0)
    for (b in seq_len(B)) {
      X  <- MASS::mvrnorm(K * I, mu, Sw)
      f1 <- data.frame(Batch = rep(sprintf("S%03d", seq_len(K)), each = I), X)
      names(f1)[-1] <- vars
      cb <- calibrate_afm_mcd(f1, vars, mcd_alpha = calibration$mcd_alpha,
                              scaling = scaling, center = center,
                              center_alpha = center_alpha)
      xb <- matrix(MASS::mvrnorm(n_new, mu, Sw / I), nrow = n_new)
      D  <- sweep(xb, 2, cb$mu_r)
      t2 <- c(t2, I * rowSums((D %*% solve(cb$Sw)) * D))
    }
    t2
  }
  t2 <- if (is.null(seed)) run() else with_local_seed(seed, run())

  UCL <- unname(stats::quantile(t2, 1 - alpha, type = 8))

  return(list(
    UCL = UCL,
    method = "Simulated (parametric bootstrap)",
    parameters = list(
      J = J,
      K = K,
      I = I,
      alpha = alpha,
      B = B,
      n_new = n_new,
      n_T2 = length(t2),
      seed = seed,
      mcd_alpha = calibration$mcd_alpha,
      scaling = scaling,
      center = center,
      center_alpha = center_alpha
    )
  ))
}
