# =====================================================================
# 15_arl0_convergence_K_v3.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The ARL0 that the analytic limit of Equation (8) actually produces
#   as a function of the number of calibration batches K, for
#   K = 30, 50, 100, 200 and 500. For each K it reports the UCL, the
#   number of false alarms, the FAR with its standard error, the ARL0
#   and its exact Poisson 95% confidence interval. It also tests
#   whether the ARL0 depends on K (homogeneity and trend) and shows how
#   the UCL approaches its asymptotic value chi2_J(1-alpha).
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - Phase 1 is CLEAN in every cell: the question is the behaviour of
#     the limit itself, not its resistance to contamination.
#   - All contamination arguments are passed EXPLICITLY, so the result
#     does not depend on package defaults that might change.
#   - The five rows share the same code and the same seeds
#     (SEED_BASE + r per replicate), so they differ only in K.
#   - The ARL0 confidence interval is the exact Poisson one on the
#     false-alarm count, which is the appropriate one for rare events.
#   - The standard error of the FAR uses the between-replicate
#     definition, the same one as Table 1.
#
# ANCHORS AND THEIR EXPECTED VALUES
#   A1  Deterministic UCL per K:
#         K=30  -> 19.6929      K=50  -> 19.1941
#         K=100 -> 18.8274      K=200 -> 18.6464
#         K=500 -> 18.5385      (tolerance 5e-04)
#   A2  The hand-written ucl_analitico() must equal the package
#       ucl_F_adjusted() at K = 100 (tolerance 1e-06).
#   A3  Reproduction of the earlier campaign: its ARL0 reference values
#       for K = 30, 50 and 100 must fall inside the new 95% Poisson
#       confidence intervals.
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   A1 and A2 halt the script if they fail. If A3 fails the results are
#   not used and the discrepancy is investigated first.
#
# OUTPUT
#   02_resultados/arl0_convergence_K_v3.csv
# =====================================================================

library(robustT2AFM)
library(parallel)

# ------------------------- Configuration -----------------------------
N_REP     <- 2000
K_GRID    <- c(30, 50, 100, 200, 500)
I         <- 20
J         <- 4
H_MCD     <- 0.67
ALPHA     <- 0.001
RHO       <- 0.6
N_F2      <- 100           # in-control Phase 2 batches per replicate
SEED_BASE <- 2026
N_CORES   <- 10            # 12 physical cores; two left free
DIR_OUT   <- "02_resultados"

VARS <- paste0("Var", 1:J)

# ----------------- Anchor A1: deterministic UCL per K ----------------
ucl_analitico <- function(K, I, J, h, alpha) {
  m  <- round(I * h)
  d2 <- K * m - K - J + 1
  (J * (K + 1) * (m - 1) / d2) * qf(1 - alpha, J, d2)
}

ANCLAS_UCL <- c(`30`  = 19.6929, `50`  = 19.1941, `100` = 18.8274,
                `200` = 18.6464, `500` = 18.5385)

cat("=== Anchor A1: deterministic UCL ===\n")
for (K in K_GRID) {
  ucl_k <- ucl_analitico(K, I, J, H_MCD, ALPHA)
  cat(sprintf("  K=%4d  computed = %.4f  | anchor = %.4f\n",
              K, ucl_k, ANCLAS_UCL[as.character(K)]))
  stopifnot(abs(ucl_k - ANCLAS_UCL[as.character(K)]) < 5e-4)
}
cat("  Anchor A1 verified: the five UCLs match.\n\n")

# --------- Anchor A2: the hand-written UCL == the package one --------
cat("=== Anchor A2: agreement with the package ucl_F_adjusted() ===\n")
sim_a2 <- simulate_batch_process(
  K1 = 100, K2 = 1, I = I, J = J, rho = RHO,
  outlier_batches_F1 = 0, prop_contam_F1 = 0, prop_ooc_F2 = 0,
  seed = SEED_BASE
)
cal_a2 <- suppressMessages(
  calibrate_afm_mcd(subset(sim_a2, Phase == "Phase 1"), VARS, mcd_alpha = H_MCD)
)
ucl_pkg  <- ucl_F_adjusted(cal_a2, I = I, alpha = ALPHA)$UCL
ucl_mano <- ucl_analitico(100, I, J, H_MCD, ALPHA)
cat(sprintf("  package = %.6f   hand-written = %.6f\n", ucl_pkg, ucl_mano))
stopifnot(abs(ucl_pkg - ucl_mano) < 1e-6)
cat("  Anchor A2 verified.\n\n")

# --------------------- One Monte Carlo replicate ---------------------
una_replica <- function(r, K, I, J, N_F2, H_MCD, RHO, VARS, UCL_K, SEED_BASE) {
  sim <- simulate_batch_process(
    K1 = K, K2 = N_F2, I = I, J = J, rho = RHO,
    outlier_batches_F1 = 0,     # Phase 1 clean
    prop_contam_F1     = 0,     # Phase 1 clean
    prop_ooc_F2        = 0,     # Phase 2 in control
    seed = SEED_BASE + r
  )
  ph1 <- subset(sim, Phase == "Phase 1")
  ph2 <- subset(sim, Phase == "Phase 2")

  cal <- suppressMessages(calibrate_afm_mcd(ph1, VARS, mcd_alpha = H_MCD))
  mon <- monitor_afm_mcd(ph2, cal, VARS)

  sum(mon$T2 > UCL_K)           # false alarms of this replicate
}

# ------------------------- Loop over K -------------------------------
cl <- makeCluster(N_CORES)
clusterEvalQ(cl, library(robustT2AFM))

resultados <- data.frame()
t_ini <- Sys.time()

for (K in K_GRID) {
  UCL_K <- ucl_analitico(K, I, J, H_MCD, ALPHA)
  cat(sprintf("--- K = %d  (UCL = %.4f, %d replicates) ---\n", K, UCL_K, N_REP))

  t0 <- Sys.time()
  alarmas <- parSapply(
    cl, seq_len(N_REP), una_replica,
    K = K, I = I, J = J, N_F2 = N_F2, H_MCD = H_MCD, RHO = RHO,
    VARS = VARS, UCL_K = UCL_K, SEED_BASE = SEED_BASE
  )
  t1 <- Sys.time()

  lotes_tot <- N_REP * N_F2
  fa_tot    <- sum(alarmas)
  far       <- fa_tot / lotes_tot
  se_far    <- sd(alarmas / N_F2) / sqrt(N_REP)   # same definition as Table 1
  arl0      <- if (fa_tot > 0) 1 / far else NA_real_

  # EXACT Poisson 95% CI on the false-alarm count
  if (fa_tot > 0) {
    lam_lo  <- qchisq(0.025, 2 * fa_tot) / 2
    lam_hi  <- qchisq(0.975, 2 * fa_tot + 2) / 2
    arl0_lo <- lotes_tot / lam_hi
    arl0_hi <- lotes_tot / lam_lo
  } else {
    arl0_lo <- lotes_tot / (qchisq(0.975, 2) / 2)
    arl0_hi <- Inf
  }

  cat(sprintf("   false alarms: %d of %d | FAR = %.6f (SE %.2e)\n",
              fa_tot, lotes_tot, far, se_far))
  cat(sprintf("   ARL0 = %.0f   Poisson 95%% CI = [%.0f, %.0f]\n",
              arl0, arl0_lo, arl0_hi))
  if (fa_tot < 30)
    cat("   WARNING: fewer than 30 alarms, the ARL0 has no resolution\n")
  cat(sprintf("   elapsed: %.1f min\n\n",
              as.numeric(difftime(t1, t0, units = "mins"))))

  resultados <- rbind(resultados, data.frame(
    K = K, N_rep = N_REP, total_batches = lotes_tot,
    false_alarms = fa_tot, FAR = far, SE_FAR = se_far,
    ARL0 = arl0, ARL0_CI95_low = arl0_lo, ARL0_CI95_high = arl0_hi,
    UCL = UCL_K
  ))
}

stopCluster(cl)
cat(sprintf("Total elapsed: %.1f min\n\n",
            as.numeric(difftime(Sys.time(), t_ini, units = "mins"))))

# --------------------------- Results ---------------------------------
cat("=== SUMMARY ===\n")
print(resultados, row.names = FALSE)

# ---- Anchor A3: does it reproduce the earlier campaign? --------------
cat("\n=== Anchor A3: does it reproduce the earlier campaign? ===\n")
if (N_REP < 500)
  cat("  WARNING: with few replicates this anchor does not discriminate.\n")

ANTERIOR <- list(`30` = 1439, `50` = 1626, `100` = 1709)
todos_ok <- TRUE
for (K in names(ANTERIOR)) {
  fila <- resultados[resultados$K == as.numeric(K), ]
  if (nrow(fila) != 1 || is.na(fila$ARL0)) next
  ref    <- ANTERIOR[[K]]
  dentro <- ref >= fila$ARL0_CI95_low && ref <= fila$ARL0_CI95_high
  todos_ok <- todos_ok && dentro
  cat(sprintf("  K=%3s  earlier = %4.0f   new = %5.0f  95%% CI [%5.0f, %6.0f]   %s\n",
              K, ref, fila$ARL0, fila$ARL0_CI95_low, fila$ARL0_CI95_high,
              if (dentro) "INSIDE" else "*** OUTSIDE ***"))
  if (fila$false_alarms < 30)
    cat("         (fewer than 30 alarms: no power to discriminate)\n")
}
cat(if (todos_ok)
  "\n  Anchor A3 verified. The paper table becomes five rows from a\n  single source. The three earlier CSV files are superseded and go\n  to 99_obsoleto.\n"
  else
  "\n  *** Anchor A3 FAILS. STOP and review before using these results. ***\n")

# ---- Does the ARL0 depend on K? -------------------------------------
cat("\n=== Does the ARL0 depend on K? ===\n")
tab <- cbind(resultados$false_alarms,
             resultados$total_batches - resultados$false_alarms)
ht  <- chisq.test(tab)
cat(sprintf("  Homogeneity of the %d measurements: X2 = %.3f, p = %.3f\n",
            nrow(resultados), ht$statistic, ht$p.value))

fa  <- resultados$false_alarms
lam <- mean(fa)
for (nm in c("K", "log(K)")) {
  sc <- if (nm == "K") resultados$K else log(resultados$K)
  z  <- sum((fa - lam) * (sc - mean(sc))) / sqrt(lam * sum((sc - mean(sc))^2))
  cat(sprintf("  Trend with %-7s: z = %+.3f, p = %.3f\n",
              nm, z, 2 * (1 - pnorm(abs(z)))))
}

far_pool <- sum(fa) / sum(resultados$total_batches)
cat(sprintf("\n  Pooling the %d measurements: ARL0 = %.0f\n",
            nrow(resultados), 1 / far_pool))
cat(sprintf("  UCL: from %.2f at K=%d to %.2f at K=%d; asymptotic chi2_%d(%.3f) = %.4f\n",
            resultados$UCL[1], resultados$K[1],
            resultados$UCL[nrow(resultados)], resultados$K[nrow(resultados)],
            J, 1 - ALPHA, qchisq(1 - ALPHA, J)))

# --- Rows ready to be copied into the paper table --------------------
cat("\n=== ROWS FOR THE PAPER TABLE (LaTeX) ===\n")
for (i in seq_len(nrow(resultados))) {
  r <- resultados[i, ]
  cat(sprintf("$%d$ & $%.2f$ & $%d$ & $%.0f$ & $[%.0f,\\ %.0f]$ \\\\\n",
              r$K, r$UCL, r$false_alarms, r$ARL0,
              r$ARL0_CI95_low, r$ARL0_CI95_high))
}

# ------------------------------ Output -------------------------------
if (!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)
ruta <- file.path(DIR_OUT, "arl0_convergence_K_v3.csv")
write.csv(resultados, ruta, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", ruta))
cat("The three earlier CSV files are NOT overwritten; move them to\n")
cat("99_obsoleto only if anchor A3 held.\n\n")
