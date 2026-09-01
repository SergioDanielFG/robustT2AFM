# =====================================================================
# 26_arl0_spatial_median_v2.R
# ---------------------------------------------------------------------
# PREDICTION REGISTERED BEFORE RUNNING
#   The spatial median has a higher breakdown value but a lower
#   efficiency than the mean at the Gaussian model, so with a clean
#   Phase 1 the ARL0 is expected to fall somewhat below the 1600 of the
#   simple average. Under the within-batch contamination of Section 3.1
#   the MCD already removes the outliers before the centers are
#   combined, so both estimators are expected to behave similarly.
#
# OUTPUT
#   02_resultados/arl0_spatial_median_v2.csv
# =====================================================================

library(robustT2AFM)
library(parallel)

# ------------------------- Configuration -----------------------------
N_REP     <- 2000
K         <- 30
I         <- 20
J         <- 4
H_MCD     <- 0.67
ALPHA     <- 0.001
RHO       <- 0.6
N_F2      <- 100           # in-control Phase 2 batches per replicate
N_CONTA   <- 6             # contaminated batches in the dirty condition
OR        <- 0.20          # outlier rate inside each contaminated batch
OS        <- 4             # outlier shift, in standard deviations
SEED_BASE <- 2026
N_CORES   <- 10
DIR_OUT   <- "02_resultados"

VARS <- paste0("Var", 1:J)

# ------------------- Spatial median (Weiszfeld) ----------------------
spatial_median <- function(X, tol = 1e-10, maxit = 2000) {
  m <- colMeans(X)
  for (it in seq_len(maxit)) {
    d <- sqrt(rowSums((X - matrix(m, nrow(X), ncol(X), byrow = TRUE))^2))
    d[d < 1e-12] <- 1e-12
    w <- 1 / d
    m_new <- colSums(X * w) / sum(w)
    if (max(abs(m_new - m)) < tol) return(m_new)
    m <- m_new
  }
  warning("Weiszfeld did not converge")
  m
}

# ----------------------- Anchors A1 and A2 ---------------------------
ucl_analitico <- function(K, I, J, h, alpha) {
  m  <- round(I * h)
  d2 <- K * m - K - J + 1
  (J * (K + 1) * (m - 1) / d2) * qf(1 - alpha, J, d2)
}

UCL <- ucl_analitico(K, I, J, H_MCD, ALPHA)
cat("=== Anchor A1: deterministic UCL at K = 30 ===\n")
cat(sprintf("  computed = %.4f | anchor = 19.6929\n", UCL))
stopifnot(abs(UCL - 19.6929) < 5e-4)

sim_a2 <- simulate_batch_process(
  K1 = K, K2 = 1, I = I, J = J, rho = RHO,
  outlier_batches_F1 = 0, outlier_rate = OR, outlier_shift = OS,
  prop_contam_F1 = 0, prop_ooc_F2 = 0, seed = SEED_BASE)
cal_a2 <- suppressMessages(
  calibrate_afm_mcd(subset(sim_a2, Phase == "Phase 1"), VARS, mcd_alpha = H_MCD))
ucl_pkg <- ucl_F_adjusted(cal_a2, I = I, alpha = ALPHA)$UCL
cat("=== Anchor A2: agreement with the package ucl_F_adjusted() ===\n")
cat(sprintf("  package = %.6f | hand-written = %.6f\n", ucl_pkg, UCL))
stopifnot(abs(ucl_pkg - UCL) < 1e-6)
cat("  Anchors A1 and A2 verified.\n\n")

# --------------------- One Monte Carlo replicate ---------------------
una_replica <- function(r, K, I, J, N_F2, H_MCD, RHO, VARS, UCL,
                        SEED_BASE, N_CONTA, OR, OS, contaminada,
                        spatial_median) {
  sim <- simulate_batch_process(
    K1 = K, K2 = N_F2, I = I, J = J, rho = RHO,
    outlier_batches_F1 = if (contaminada) N_CONTA else 0,
    outlier_rate       = OR,
    outlier_shift      = OS,
    prop_contam_F1     = 0,
    prop_ooc_F2        = 0,
    seed = SEED_BASE + r)
  ph1 <- subset(sim, Phase == "Phase 1")
  ph2 <- subset(sim, Phase == "Phase 2")
  
  cal <- suppressMessages(calibrate_afm_mcd(ph1, VARS, mcd_alpha = H_MCD))
  
  # (a) published estimator: simple average of the K robust centers
  fa_mean <- sum(monitor_afm_mcd(ph2, cal, VARS)$T2 > UCL)
  
  # (b) spatial median of those same K centers, standardized
  C      <- do.call(rbind, cal$mcd_centers)[, VARS, drop = FALSE]
  sd_ph1 <- apply(ph1[, VARS], 2, sd)
  cal_sm <- cal
  cal_sm$mu_r <- spatial_median(sweep(C, 2, sd_ph1, "/")) * sd_ph1
  fa_sm <- sum(monitor_afm_mcd(ph2, cal_sm, VARS)$T2 > UCL)
  
  c(fa_mean, fa_sm)
}

# --------------------------- Both conditions -------------------------
cl <- makeCluster(N_CORES)
clusterEvalQ(cl, library(robustT2AFM))

resultados <- data.frame()
t_ini <- Sys.time()

for (contaminada in c(FALSE, TRUE)) {
  etiqueta <- if (contaminada) "contaminated" else "clean"
  cat(sprintf("--- Phase 1 %s (K = %d, UCL = %.4f, %d replicates) ---\n",
              etiqueta, K, UCL, N_REP))
  t0 <- Sys.time()
  
  out <- parSapply(
    cl, seq_len(N_REP), una_replica,
    K = K, I = I, J = J, N_F2 = N_F2, H_MCD = H_MCD, RHO = RHO,
    VARS = VARS, UCL = UCL, SEED_BASE = SEED_BASE,
    N_CONTA = N_CONTA, OR = OR, OS = OS, contaminada = contaminada,
    spatial_median = spatial_median)
  
  lotes_tot <- N_REP * N_F2
  
  for (k in 1:2) {
    est     <- c("mean", "spatial_median")[k]
    alarmas <- out[k, ]
    fa_tot  <- sum(alarmas)
    far     <- fa_tot / lotes_tot
    se_far  <- sd(alarmas / N_F2) / sqrt(N_REP)
    arl0    <- if (fa_tot > 0) 1 / far else NA_real_
    
    if (fa_tot > 0) {
      arl0_lo <- lotes_tot / (qchisq(0.975, 2 * fa_tot + 2) / 2)
      arl0_hi <- lotes_tot / (qchisq(0.025, 2 * fa_tot) / 2)
    } else {
      arl0_lo <- lotes_tot / (qchisq(0.975, 2) / 2); arl0_hi <- Inf
    }
    
    cat(sprintf("   %-15s false alarms: %5d of %d | ARL0 = %7.0f  CI95 = [%.0f, %.0f]\n",
                est, fa_tot, lotes_tot, arl0, arl0_lo, arl0_hi))
    
    resultados <- rbind(resultados, data.frame(
      phase1 = etiqueta, center_estimator = est, K = K, N_rep = N_REP,
      total_batches = lotes_tot, false_alarms = fa_tot,
      FAR = far, SE_FAR = se_far, ARL0 = arl0,
      ARL0_CI95_low = arl0_lo, ARL0_CI95_high = arl0_hi, UCL = UCL))
  }
  
  cat(sprintf("   elapsed: %.1f min\n\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

stopCluster(cl)

# ------------------------- Anchors A3 and A4 -------------------------
fa_ref <- resultados$false_alarms[resultados$phase1 == "clean" &
                                    resultados$center_estimator == "mean"]
cat("=== Anchor A3: reproduction of the K = 30 row of Table 2 ===\n")
cat(sprintf("  false alarms with the simple average, clean Phase 1: %d (expected 125)\n",
            fa_ref))
if (fa_ref != 125)
  cat("  *** DOES NOT REPRODUCE. Do not use these figures until understood. ***\n")
stopifnot(fa_ref == 125)
cat("  Anchor A3 verified.\n\n")

arl_cont <- resultados$ARL0[resultados$phase1 == "contaminated" &
                              resultados$center_estimator == "mean"]
cat("=== Anchor A4: contaminated Phase 1 with the simple average ===\n")
cat(sprintf("  ARL0 = %.0f (must fall inside [1211, 1698])\n", arl_cont))
if (!(arl_cont > 1211 && arl_cont < 1698))
  cat("  *** OUTSIDE THE PUBLISHED INTERVAL. The contamination mechanism does not match. ***\n")
stopifnot(arl_cont > 1211 && arl_cont < 1698)
cat("  Anchor A4 verified.\n\n")

# --------------------------- Results ---------------------------------
cat("=== SUMMARY ===\n")
print(resultados, row.names = FALSE)

cat("\n=== PAIRED COMPARISON ===\n")
for (cond in c("clean", "contaminated")) {
  s <- resultados[resultados$phase1 == cond, ]
  cat(sprintf("  Phase 1 %-12s : mean %.0f vs spatial median %.0f  (ratio %.2f)\n",
              cond, s$ARL0[s$center_estimator == "mean"],
              s$ARL0[s$center_estimator == "spatial_median"],
              s$ARL0[s$center_estimator == "spatial_median"] /
                s$ARL0[s$center_estimator == "mean"]))
}

cat(sprintf("\nTotal elapsed: %.1f min\n",
            as.numeric(difftime(Sys.time(), t_ini, units = "mins"))))

if (!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)
write.csv(resultados, file.path(DIR_OUT, "arl0_spatial_median_v2.csv"),
          row.names = FALSE)
cat("Saved: 02_resultados/arl0_spatial_median_v2.csv\n")
cat("\n=== END. Copy the whole output. ===\n\n")