# =====================================================================
# 23_REFERENCE_CENTER_AUDIT_v4.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   An audit of the reference center over six Phase 1 contamination
#   scenarios, answering two questions:
#     (1) Do the AFM weights, built from lambda_1 of S_k^MCD, identify
#         batches that are shifted as a whole when their covariance
#         does not change?
#     (2) Does weighting the center with those same weights improve on
#         the current center?
#             mu_S = (1/K) sum_k mu_k^MCD = cal$mu_r   (simple average)
#             mu_W = sum_k w_k mu_k^MCD               (weighted)
#   For each scenario it reports how many of the six evaluated batches
#   fall among the six lowest weights, the healthy / contaminated weight
#   ratio, the Mahalanobis distance of each center to the true center,
#   and the FAR and TPR each center produces, with standard errors and
#   95% confidence intervals for the differences.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - The published method remains (mu_S, Sw). mu_W is evaluated only as
#     an exploratory variant.
#   - Both variants use the SAME Sw and the SAME UCL, so any difference
#     is attributable to the center and not to the covariance or to the
#     limit.
#   - Scenarios 2 to 4 (pure shift of 1, 3 and 6 SD) share the same
#     seed, so the simulator picks the same six batches and generates
#     the same random cloud; only the translation changes. That is what
#     makes the three rows comparable to each other.
#   - Scenario 5 is built by hand because the simulator cannot give
#     2*Sigma to the contaminated batches only.
#   - In the clean scenario six batches are pseudo-labelled with no
#     contamination applied, which gives the chance level against which
#     the weight diagnostics are read; that level is the hypergeometric
#     expectation N_CONTA^2/K1.
#   - Each replicate seeds its own stream with SEED_BASE*1200 + idx, so
#     the result does not depend on how work is split between workers.
#
# ANCHORS AND THEIR EXPECTED VALUES
#   ANCHOR 0  MCD translation equivariance: a pure translation
#             X -> X + 3 must leave the MCD covariance and lambda_1
#             unchanged and move the center by exactly 3
#             (tolerance 1e-07 on all three). The same seed is forced
#             before both covMcd() calls so the algorithm uses the same
#             internal randomness.
#   ANCHOR 1  Base configuration of the manuscript:
#               contaminated batches detected = 6
#               calibrate_afm_mcd() returns weights, mu_r, Sw,
#                 mcd_centers, mcd_covariances and lambda1
#               UCL = 19.69285  (tolerance 1e-04)
#               max |mean(mcd_centers) - cal$mu_r| < 1e-08
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   Both anchors stop the script on failure; the campaign is not run
#   otherwise. If the mechanism is purely locational, n_in_lowest and
#   ratio_w must be practically the same at 1, 3 and 6 SD and close to
#   the chance control.
#
# OUTPUT
#   02_resultados/weighted_center_audit_v4.csv
#   02_resultados/center_audit_v4_detail_scenario_0..5.csv
# =====================================================================

library(parallel)
library(robustT2AFM)
library(MASS)
library(robustbase)

dir.create("02_resultados", showWarnings = FALSE)

# =====================================================================
# 1. PARAMETERS
# =====================================================================

N_REP     <- 2000
K1        <- 30
I         <- 20
J         <- 4
RHO       <- 0.6
H_MCD     <- 0.67
ALPHA     <- 0.001

N_CONTA   <- 6
OR        <- 0.20
OS        <- 4

N_F2_SANO <- 100
N_F2_OOC  <- 100
DELTA_F2  <- 1.0

SEED_BASE <- 2026
N_WORKERS <- 8

VARS <- paste0("Var", seq_len(J))
AZAR <- N_CONTA * N_CONTA / K1  # hypergeometric expectation = 1.2

make_equi <- function(p, rho) {
  S <- matrix(rho, p, p)
  diag(S) <- 1
  S
}

SIGMA     <- make_equi(J, RHO)
SIGINV    <- solve(SIGMA)
SIGMA_SD  <- sqrt(diag(SIGMA))

maha <- function(mu) {
  mu <- as.numeric(mu)
  as.numeric(sqrt(t(mu) %*% SIGINV %*% mu))
}

# =====================================================================
# 2. HELPER FUNCTIONS
# =====================================================================

extrae_ucl <- function(x) {
  if (is.numeric(x) && length(x) == 1L) return(as.numeric(x))

  if (is.list(x)) {
    for (nn in c("UCL", "ucl", "UCL_F", "ucl_F",
                 "limit", "UCL_robust", "value")) {
      v <- x[[nn]]
      if (!is.null(v) && is.numeric(v) && length(v) >= 1L)
        return(as.numeric(v[1]))
    }
  }

  stop("Could not extract the UCL from ucl_F_adjusted().")
}

lotes_contaminados <- function(f1) {
  if (!"ContaminationType" %in% names(f1))
    stop("The ContaminationType column is missing.")

  tipo <- as.character(f1$ContaminationType)
  unique(as.character(f1$Batch[tipo != "Clean"]))
}

extraer_centros_mcd <- function(cal, nombres_lotes) {
  if (is.null(cal$mcd_centers))
    stop("calibrate_afm_mcd() does not return cal$mcd_centers.")

  mc <- cal$mcd_centers

  if (is.list(mc)) {
    if (is.null(names(mc)))
      stop("cal$mcd_centers is an unnamed list.")

    if (!all(nombres_lotes %in% names(mc)))
      stop("Not every batch in weights appears in mcd_centers.")

    M <- do.call(rbind, mc[nombres_lotes])
    M <- as.matrix(M)
    rownames(M) <- nombres_lotes
    return(M)
  }

  M <- as.matrix(mc)

  if (!is.null(rownames(M))) {
    if (!all(nombres_lotes %in% rownames(M)))
      stop("Not every batch appears in the mcd_centers matrix.")
    M <- M[nombres_lotes, , drop = FALSE]
  } else {
    if (nrow(M) != length(nombres_lotes))
      stop("The order of mcd_centers cannot be determined.")
    rownames(M) <- nombres_lotes
  }

  M
}

# =====================================================================
# 3. ANCHOR 0: MCD TRANSLATION EQUIVARIANCE
# =====================================================================

cat("\n====================================================\n")
cat(" ANCHOR 0: MCD TRANSLATION EQUIVARIANCE\n")
cat("====================================================\n")

set.seed(111)
X0 <- MASS::mvrnorm(I, rep(0, J), SIGMA)
X3 <- sweep(X0, 2, rep(3, J), "+")

set.seed(222)
m0 <- robustbase::covMcd(X0, alpha = H_MCD)

set.seed(222)
m3 <- robustbase::covMcd(X3, alpha = H_MCD)

dif_cov <- max(abs(m0$cov - m3$cov))
dif_cen <- max(abs((m3$center - m0$center) - rep(3, J)))
lam0 <- max(eigen(m0$cov, symmetric = TRUE, only.values = TRUE)$values)
lam3 <- max(eigen(m3$cov, symmetric = TRUE, only.values = TRUE)$values)

cat(sprintf("max |Cov_MCD(X) - Cov_MCD(X+3)| = %.3e\n", dif_cov))
cat(sprintf("max |(center_X+3-center_X)-3|   = %.3e\n", dif_cen))
cat(sprintf("lambda1 original = %.10f | translated = %.10f\n", lam0, lam3))

if (dif_cov > 1e-7 || dif_cen > 1e-7 || abs(lam0 - lam3) > 1e-7) {
  stop("The translation anchor fails. Do not run the campaign.")
}

cat("ANCHOR 0 OK.\n")

# =====================================================================
# 4. PACKAGE ANCHORS IN THE BASE CONFIGURATION OF THE MANUSCRIPT
# =====================================================================

cat("\n====================================================\n")
cat(" ANCHOR 1: BASE CONFIGURATION OF THE MANUSCRIPT\n")
cat("====================================================\n")

set.seed(SEED_BASE * 1200 + 1)

prueba <- simulate_batch_process(
  K1 = K1, K2 = 1, I = I, J = J,
  rho = RHO, Sigma = SIGMA,
  outlier_batches_F1 = N_CONTA,
  outlier_rate = OR,
  outlier_shift = OS,
  prop_contam_F1 = 0,
  prop_ooc_F2 = 0
)

f1p <- subset(prueba, Phase == "Phase 1")
contp <- lotes_contaminados(f1p)

cat(sprintf("Contaminated batches detected: %d (expected %d)\n",
            length(contp), N_CONTA))
stopifnot(length(contp) == N_CONTA)

calp <- suppressMessages(
  calibrate_afm_mcd(f1p, VARS, mcd_alpha = H_MCD)
)

requeridos <- c("weights", "mu_r", "Sw",
                "mcd_centers", "mcd_covariances", "lambda1")
faltan <- setdiff(requeridos, names(calp))

if (length(faltan) > 0)
  stop("Missing components in calibrate_afm_mcd(): ",
       paste(faltan, collapse = ", "))

stopifnot(length(calp$weights) == K1)

ucl_p <- extrae_ucl(
  ucl_F_adjusted(calp, I = I, alpha = ALPHA)
)

cat(sprintf("UCL obtained = %.5f | expected = 19.69285\n", ucl_p))
if (abs(ucl_p - 19.69285) > 1e-4)
  stop("UCL ANCHOR FAILED.")

nm_p <- names(calp$weights)
centros_p <- extraer_centros_mcd(calp, nm_p)

mu_manual <- colMeans(centros_p)
mu_paquete <- as.numeric(calp$mu_r)
dif_mu <- max(abs(mu_manual - mu_paquete))

cat(sprintf("max |mean(mcd_centers) - cal$mu_r| = %.3e\n", dif_mu))
if (dif_mu > 1e-8)
  stop("cal$mu_r does not match the average of mcd_centers.")

cat("ANCHOR 1 OK.\n")

# =====================================================================
# 5. PHASE 1 GENERATORS
# =====================================================================
#
# 0. clean:
#    thirty clean batches are generated and six are pseudo-labelled
#    only to measure the chance level of n_in_lowest and ratio_w.
#
# 1. outliers:
#    exactly the mechanism of the manuscript: 6 batches with 20% of
#    their observations shifted 4 SD.
#
# 2-4. shift:
#    uses prop_contam_F1 of the package DIRECTLY. It is a pure shift:
#    it changes the batch mean and leaves Sigma untouched.
#
# 5. shift_cov2:
#    built by hand because simulate_batch_process() cannot give 2*Sigma
#    to the contaminated batches only.
# =====================================================================

genera_f1 <- function(esc, d = 0, factor_cov = 1) {

  if (esc == "limpia") {
    sim <- simulate_batch_process(
      K1 = K1, K2 = 0, I = I, J = J,
      rho = RHO, Sigma = SIGMA,
      outlier_batches_F1 = 0,
      prop_contam_F1 = 0,
      prop_ooc_F2 = 0
    )

    f1 <- subset(sim, Phase == "Phase 1")
    b <- unique(as.character(f1$Batch))
    pseudo <- b[seq_len(N_CONTA)]

    return(list(
      datos = f1[, c("Batch", VARS)],
      cont = pseudo,
      etiqueta_real = FALSE
    ))
  }

  if (esc == "atipicos") {
    sim <- simulate_batch_process(
      K1 = K1, K2 = 0, I = I, J = J,
      rho = RHO, Sigma = SIGMA,
      outlier_batches_F1 = N_CONTA,
      outlier_rate = OR,
      outlier_shift = OS,
      prop_contam_F1 = 0,
      prop_ooc_F2 = 0
    )

    f1 <- subset(sim, Phase == "Phase 1")
    return(list(
      datos = f1[, c("Batch", VARS)],
      cont = lotes_contaminados(f1),
      etiqueta_real = TRUE
    ))
  }

  if (esc == "shift") {
    sim <- simulate_batch_process(
      K1 = K1, K2 = 0, I = I, J = J,
      rho = RHO, Sigma = SIGMA,
      outlier_batches_F1 = 0,
      prop_contam_F1 = N_CONTA / K1,
      shift_contam = d,
      prop_ooc_F2 = 0
    )

    f1 <- subset(sim, Phase == "Phase 1")
    return(list(
      datos = f1[, c("Batch", VARS)],
      cont = lotes_contaminados(f1),
      etiqueta_real = TRUE
    ))
  }

  if (esc == "shift_cov2") {
    nb <- paste0("F1_B", sprintf("%02d", seq_len(K1)))
    malos <- nb[seq_len(N_CONTA)]
    filas <- vector("list", K1)

    for (k in seq_len(K1)) {
      es_malo <- nb[k] %in% malos
      mu_k <- if (es_malo) d * SIGMA_SD else rep(0, J)
      S_k  <- if (es_malo) factor_cov * SIGMA else SIGMA

      X <- MASS::mvrnorm(I, mu_k, S_k)
      colnames(X) <- VARS

      filas[[k]] <- data.frame(
        Batch = nb[k],
        X,
        stringsAsFactors = FALSE
      )
    }

    return(list(
      datos = do.call(rbind, filas),
      cont = malos,
      etiqueta_real = TRUE
    ))
  }

  stop("Unknown scenario: ", esc)
}

# =====================================================================
# 6. ONE REPLICATE
# =====================================================================

one_rep <- function(idx, esc, d, factor_cov) {

  # Same seed for d = 1, 3 and 6: the simulator picks the same six
  # batches and generates the same random cloud; only the translation
  # changes.
  set.seed(SEED_BASE * 1200 + idx)

  g <- genera_f1(esc, d, factor_cov)
  f1 <- g$datos
  cont <- as.character(g$cont)

  if (length(cont) != N_CONTA)
    stop("Expected ", N_CONTA, " batches in the evaluated group, found ",
         length(cont), ".")

  cal <- suppressMessages(
    calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD)
  )

  if (length(cal$weights) != K1)
    stop("The calibration did not retain the ", K1, " batches.")

  w <- cal$weights
  nm <- names(w)

  if (!all(cont %in% nm))
    stop("The labels of the evaluated group do not match names(weights).")

  # ---------------- WEIGHTS ----------------
  seis_menores <- nm[order(w)][seq_len(N_CONTA)]
  n_abajo <- sum(seis_menores %in% cont)

  es_cont <- nm %in% cont
  ratio_w <- mean(w[!es_cont]) / mean(w[es_cont])

  # ---------------- CENTERS ----------------
  centros <- extraer_centros_mcd(cal, nm)

  mu_S <- as.numeric(cal$mu_r)
  mu_check <- colMeans(centros)

  if (max(abs(mu_S - mu_check)) > 1e-8)
    stop("mu_r does not match the average of mcd_centers.")

  mu_W <- as.numeric(t(centros) %*% as.numeric(w[nm]))

  D_S <- maha(mu_S)
  D_W <- maha(mu_W)

  # ---------------- SAME Sw / SAME UCL ----------------
  SwI <- solve(cal$Sw)
  ucl <- extrae_ucl(ucl_F_adjusted(cal, I = I, alpha = ALPHA))

  t2_lote <- function(Xb, mu) {
    dd <- colMeans(Xb) - mu
    as.numeric(I * t(dd) %*% SwI %*% dd)
  }

  genera_lote_f2 <- function(delta) {
    X <- MASS::mvrnorm(
      n = I,
      mu = delta * SIGMA_SD,
      Sigma = SIGMA
    )
    colnames(X) <- VARS
    X
  }

  sanos <- lapply(seq_len(N_F2_SANO), function(z) genera_lote_f2(0))
  oocs  <- lapply(seq_len(N_F2_OOC),  function(z) genera_lote_f2(DELTA_F2))

  t_s_sanos <- vapply(sanos, t2_lote, numeric(1), mu = mu_S)
  t_w_sanos <- vapply(sanos, t2_lote, numeric(1), mu = mu_W)
  t_s_ooc   <- vapply(oocs,  t2_lote, numeric(1), mu = mu_S)
  t_w_ooc   <- vapply(oocs,  t2_lote, numeric(1), mu = mu_W)

  FAR_S <- mean(t_s_sanos > ucl)
  FAR_W <- mean(t_w_sanos > ucl)
  TPR_S <- mean(t_s_ooc > ucl)
  TPR_W <- mean(t_w_ooc > ucl)

  c(
    n_in_lowest = n_abajo,
    ratio_w = ratio_w,
    D_S = D_S,
    D_W = D_W,
    FAR_S = FAR_S,
    FAR_W = FAR_W,
    TPR_S = TPR_S,
    TPR_W = TPR_W,
    delta_D = D_S - D_W,
    delta_FAR = FAR_W - FAR_S,
    delta_TPR = TPR_W - TPR_S
  )
}

# =====================================================================
# 7. SCENARIOS
# =====================================================================

ESC <- list(
  list(nombre = "0. Limpia (control negativo)",
       esc = "limpia", shift = 0, factor_cov = 1),

  list(nombre = "1. Atipicos internos (paper)",
       esc = "atipicos", shift = 0, factor_cov = 1),

  list(nombre = "2. Desplazamiento puro 1 SD",
       esc = "shift", shift = 1, factor_cov = 1),

  list(nombre = "3. Desplazamiento puro 3 SD",
       esc = "shift", shift = 3, factor_cov = 1),

  list(nombre = "4. Desplazamiento puro 6 SD",
       esc = "shift", shift = 6, factor_cov = 1),

  list(nombre = "5. Desplaz. 3 SD + covarianza x2",
       esc = "shift_cov2", shift = 3, factor_cov = 2)
)

# =====================================================================
# 8. CLUSTER
# =====================================================================

cat("\n====================================================\n")
cat(" STARTING THE MONTE CARLO CAMPAIGN\n")
cat("====================================================\n")
cat(sprintf("N_REP = %d | F2 healthy = %d | F2 OOC = %d | workers = %d\n\n",
            N_REP, N_F2_SANO, N_F2_OOC, N_WORKERS))

cl <- makeCluster(N_WORKERS)
on.exit(stopCluster(cl), add = TRUE)

clusterEvalQ(cl, {
  library(robustT2AFM)
  library(MASS)
  library(robustbase)
  NULL
})

clusterExport(
  cl,
  c("K1", "I", "J", "RHO", "H_MCD", "ALPHA",
    "N_CONTA", "OR", "OS",
    "N_F2_SANO", "N_F2_OOC", "DELTA_F2",
    "SEED_BASE", "VARS",
    "SIGMA", "SIGINV", "SIGMA_SD",
    "maha", "extrae_ucl", "lotes_contaminados",
    "extraer_centros_mcd", "genera_f1", "one_rep")
)

# =====================================================================
# 9. SUMMARY
# =====================================================================

media <- function(x) mean(x, na.rm = TRUE)

se <- function(x) {
  x <- x[is.finite(x)]
  sd(x) / sqrt(length(x))
}

ic95 <- function(x) {
  mm <- media(x)
  ss <- se(x)
  c(mm - 1.96 * ss, mm + 1.96 * ss)
}

resumen <- data.frame()
detalle <- list()
t_ini <- Sys.time()

for (e in seq_along(ESC)) {

  E <- ESC[[e]]
  cat("\n----------------------------------------------------\n")
  cat(E$nombre, "\n")
  cat("----------------------------------------------------\n")

  t0 <- Sys.time()

  out <- parLapplyLB(
    cl,
    seq_len(N_REP),
    one_rep,
    esc = E$esc,
    d = E$shift,
    factor_cov = E$factor_cov
  )

  M <- as.data.frame(do.call(rbind, out))
  detalle[[e]] <- M

  icD   <- ic95(M$delta_D)
  icFAR <- ic95(M$delta_FAR)
  icTPR <- ic95(M$delta_TPR)

  fila <- data.frame(
    scenario = E$nombre,
    shift_SD = E$shift,
    factor_cov = E$factor_cov,

    n_in_lowest = media(M$n_in_lowest),
    SE_n_in_lowest = se(M$n_in_lowest),
    chance = AZAR,

    ratio_w = media(M$ratio_w),
    SE_ratio_w = se(M$ratio_w),

    D_simple = media(M$D_S),
    D_weighted = media(M$D_W),
    gain_D = media(M$delta_D),
    SE_gain_D = se(M$delta_D),
    CI95_D_low = icD[1],
    CI95_D_high = icD[2],

    FAR_simple = media(M$FAR_S),
    FAR_weighted = media(M$FAR_W),
    FAR_difference = media(M$delta_FAR),
    SE_FAR_difference = se(M$delta_FAR),
    CI95_FAR_low = icFAR[1],
    CI95_FAR_high = icFAR[2],

    TPR_simple = media(M$TPR_S),
    TPR_weighted = media(M$TPR_W),
    gain_TPR = media(M$delta_TPR),
    SE_gain_TPR = se(M$delta_TPR),
    CI95_TPR_low = icTPR[1],
    CI95_TPR_high = icTPR[2],

    stringsAsFactors = FALSE
  )

  resumen <- rbind(resumen, fila)

  cat(sprintf("n_in_lowest = %.3f of 6 | chance = %.1f | ratio_w = %.3f\n",
              fila$n_in_lowest, AZAR, fila$ratio_w))
  cat(sprintf("D: simple %.4f | weighted %.4f | gain %+.4f | 95%% CI [%+.4f,%+.4f]\n",
              fila$D_simple, fila$D_weighted, fila$gain_D,
              fila$CI95_D_low, fila$CI95_D_high))
  cat(sprintf("FAR: %.6f -> %.6f | delta %+.6f | 95%% CI [%+.6f,%+.6f]\n",
              fila$FAR_simple, fila$FAR_weighted, fila$FAR_difference,
              fila$CI95_FAR_low, fila$CI95_FAR_high))
  cat(sprintf("TPR: %.4f -> %.4f | gain %+.4f | 95%% CI [%+.4f,%+.4f]\n",
              fila$TPR_simple, fila$TPR_weighted, fila$gain_TPR,
              fila$CI95_TPR_low, fila$CI95_TPR_high))
  cat(sprintf("Scenario elapsed: %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# =====================================================================
# 10. KEY CHECK: d = 1, 3 and 6 MUST GIVE THE SAME WEIGHTS ON AVERAGE
# =====================================================================

cat("\n====================================================\n")
cat(" PURE-SHIFT CHECK\n")
cat("====================================================\n")

rp <- resumen[3:5, c("shift_SD", "n_in_lowest", "ratio_w")]
print(rp, row.names = FALSE)

cat("\nIf the mechanism is purely locational, n_in_lowest and ratio_w must be\n")
cat("practically the same at 1, 3 and 6 SD and close to the chance control.\n")

# =====================================================================
# 11. SAVING
# =====================================================================

ruta <- file.path("02_resultados", "weighted_center_audit_v4.csv")
write.csv(resumen, ruta, row.names = FALSE)

for (e in seq_along(detalle)) {
  write.csv(
    detalle[[e]],
    file.path("02_resultados",
              sprintf("center_audit_v4_detail_scenario_%d.csv", e - 1)),
    row.names = FALSE
  )
}

cat("\n====================================================\n")
cat(" FINAL RESULT\n")
cat("====================================================\n")
print(resumen, row.names = FALSE)

cat(sprintf("\nTotal elapsed: %.1f min\n",
            as.numeric(difftime(Sys.time(), t_ini, units = "mins"))))
cat("\nSaved to:\n", ruta, "\n")
cat("=== END ===\n")
