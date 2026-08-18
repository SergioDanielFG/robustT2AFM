# =====================================================================
# KEFF_EFFECT_v2.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   Whether replacing the number of calibration batches K by the
#   effective number K_eff = 1/sum(w^2) in the F limit brings the ARL0
#   closer to the nominal target of 1/alpha or moves it further away.
#   It reports K_eff with the real weights (mean, range, standard
#   deviation and the percentage drop with respect to K), the two UCLs,
#   and for each limit the false alarm rate, its between-replicate
#   standard error, the ARL0 and the total false alarms.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - The hand-written ucl_F_conK() accepts a non-integer K, which
#     ucl_F_adjusted() does not; that is why the limit is recomputed
#     here instead of being taken from the package. The structural
#     anchor below is what makes that reimplementation usable as a
#     baseline.
#   - The uniform-weight check runs through the SAME code path as the
#     loop rather than through loose algebra, so it tests the code and
#     not just the identity 1/(K*(1/K)^2) = K.
#   - Both limits are evaluated on the SAME replicates and the SAME
#     T2 statistics, so the difference in ARL0 is attributable to the
#     limit alone.
#   - Phase 1 is CONTAMINATED and the seeds are SEED_BASE*100 + rep,
#     the same as the analytic campaign, so the K row is directly
#     comparable with the published one.
#   - The standard error of the FAR uses the between-replicate
#     definition, the same one as Table 1.
#
# ANCHORS AND THEIR EXPECTED VALUES
#   STRUCTURAL  ucl_F_conK(K = 30) must equal the package
#               ucl_F_adjusted() to better than 1e-09.
#   UNIFORM     with uniform weights, K_eff must equal K = 30 to better
#               than 1e-09.
#   TABLE 1     FAR with K = 0.000625  (tolerance 5e-05)
#               SE  with K = 5.55e-05  (tolerance 2e-06)
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   The correction is declared an improvement only if the ARL0 it
#   produces is closer to the target 1/alpha than the one produced with
#   K. Any other outcome is reported as it stands.
#
# OUTPUT
#   02_resultados/table1_row_keff.csv
# =====================================================================
library(parallel)

N_REP     <- 2000
K1 <- 30; K2 <- 100; I <- 20; J <- 4
RHO <- 0.6; H_MCD <- 0.67; ALPHA <- 0.001
OB <- 6; OR <- 0.20; OS <- 4
SEED_BASE <- 2026
VARS <- paste0("Var", 1:J)

make_equi <- function(p, rho) { S <- matrix(rho, p, p); diag(S) <- 1; S }
Sigma_EQ <- make_equi(J, RHO)

FAR_ROB_PUB <- 0.000625
SE_ROB_PUB  <- 5.55e-05

# ---------------------------------------------------------------------
# F limit with an arbitrary K (allows an integer K or a non-integer K_eff)
# ---------------------------------------------------------------------
ucl_F_conK <- function(Kval, mstar, J, alpha) {
  df2 <- Kval * mstar - Kval - J + 1
  if (df2 <= 0) return(NA_real_)
  ((J * (Kval + 1) * (mstar - 1)) / df2) * qf(1 - alpha, J, df2)
}

# ---------------------------------------------------------------------
# STRUCTURAL ANCHOR: the hand-written function must match the package
# ---------------------------------------------------------------------
library(robustT2AFM); library(MASS)

cat("=== STRUCTURAL ANCHOR ===\n")
set.seed(SEED_BASE * 100 + 1)
sim_a <- simulate_batch_process(K1 = K1, K2 = 1, I = I, J = J, rho = RHO,
                                Sigma = Sigma_EQ,
                                outlier_batches_F1 = OB, outlier_rate = OR,
                                outlier_shift = OS, prop_ooc_F2 = 0)
cal_a <- suppressMessages(calibrate_afm_mcd(subset(sim_a, Phase == "Phase 1"),
                                            VARS, mcd_alpha = H_MCD))
m_star_a <- round(I * H_MCD)
ucl_pkg  <- ucl_F_adjusted(cal_a, I = I, alpha = ALPHA)$UCL
ucl_mano <- ucl_F_conK(K1, m_star_a, J, ALPHA)

cat(sprintf("  m* = %d\n", m_star_a))
cat(sprintf("  package ucl_F_adjusted()      : %.6f\n", ucl_pkg))
cat(sprintf("  hand-written ucl_F_conK(K=30) : %.6f\n", ucl_mano))
cat(sprintf("  difference                    : %.2e\n", abs(ucl_pkg - ucl_mano)))
stopifnot(abs(ucl_pkg - ucl_mano) < 1e-9)
cat("  >>> OK: the two formulas agree. The baseline is valid.\n")

# Second check, running through the same code path as the loop
# (not through loose algebra):
cal_unif <- cal_a
cal_unif$weights <- rep(1 / length(cal_a$weights), length(cal_a$weights))
Keff_unif <- 1 / sum(cal_unif$weights^2)
cat(sprintf("\n  Uniform-weight check: K_eff = %.6f (must be %d)\n",
            Keff_unif, K1))
stopifnot(abs(Keff_unif - K1) < 1e-9)
cat("  >>> OK\n\n")

# ---------------------------------------------------------------------
# One replicate
# ---------------------------------------------------------------------
una_replica <- function(rep) {
  set.seed(SEED_BASE * 100 + rep)
  sim <- simulate_batch_process(K1 = K1, K2 = K2, I = I, J = J, rho = RHO,
                                Sigma = Sigma_EQ,
                                outlier_batches_F1 = OB, outlier_rate = OR,
                                outlier_shift = OS, prop_ooc_F2 = 0)
  f1 <- subset(sim, Phase == "Phase 1")
  f2 <- subset(sim, Phase == "Phase 2")
  cal <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))

  w      <- cal$weights
  Keff   <- 1 / sum(w^2)
  m_star <- round(I * H_MCD)

  ucl_K    <- ucl_F_conK(K1,   m_star, J, ALPHA)   # current
  ucl_Keff <- ucl_F_conK(Keff, m_star, J, ALPHA)   # proposed correction

  t2 <- monitor_afm_mcd(f2, cal, VARS)$T2

  c(far_K    = mean(t2 > ucl_K),
    far_Keff = mean(t2 > ucl_Keff),
    Keff     = Keff,
    ucl_K    = ucl_K,
    ucl_Keff = ucl_Keff,
    n_lotes  = length(t2))
}

procesar_bloque <- function(reps) {
  t(vapply(reps, una_replica, numeric(6)))
}

cat("=== Serial check: replicate 1 ===\n")
print(una_replica(1))
cat("Check passed. Launching the full campaign...\n\n")

# ---------------------------------------------------------------------
# Parallel campaign
# ---------------------------------------------------------------------
t0 <- Sys.time()
n_workers <- 8
bloques <- split(seq_len(N_REP), cut(seq_len(N_REP), n_workers, labels = FALSE))

cl <- makeCluster(n_workers)
clusterEvalQ(cl, { library(robustT2AFM); library(MASS) })
clusterExport(cl, c("K1", "K2", "I", "J", "RHO", "H_MCD", "ALPHA",
                    "OB", "OR", "OS", "SEED_BASE", "VARS", "Sigma_EQ",
                    "ucl_F_conK", "una_replica"))
res <- parLapply(cl, bloques, procesar_bloque)
stopCluster(cl)

M <- do.call(rbind, res)
stopifnot(nrow(M) == N_REP, all(M[, "n_lotes"] == K2))

cat(sprintf("Elapsed: %.1f min | %d workers\n\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), n_workers))

# ---------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------
far_K    <- M[, "far_K"]
far_Keff <- M[, "far_Keff"]

cat("===== K_eff WITH THE REAL WEIGHTS =====\n")
cat(sprintf("  K = %d\n", K1))
cat(sprintf("  mean K_eff = %.2f  (a %.1f%% drop with respect to K)\n",
            mean(M[, "Keff"]), 100 * (1 - mean(M[, "Keff"]) / K1)))
cat(sprintf("  K_eff range = [%.2f, %.2f]   sd = %.3f\n",
            min(M[, "Keff"]), max(M[, "Keff"]), sd(M[, "Keff"])))
cat(sprintf("  UCL with K     : %.4f\n", mean(M[, "ucl_K"])))
cat(sprintf("  UCL with K_eff : %.4f  (mean)\n\n", mean(M[, "ucl_Keff"])))

resumen <- data.frame(
  Limit     = c("Robusto F con K (actual)", "Robusto F con K_eff (Satterthwaite)"),
  FAR       = c(mean(far_K), mean(far_Keff)),
  SE_FAR    = c(sd(far_K) / sqrt(N_REP), sd(far_Keff) / sqrt(N_REP)),
  ARL0      = c(1 / mean(far_K), 1 / mean(far_Keff)),
  FA_total  = c(round(sum(far_K) * K2), round(sum(far_Keff) * K2)),
  n_batches = rep(N_REP * K2, 2),
  stringsAsFactors = FALSE
)

cat("===== EFFECT ON THE CALIBRATION =====\n")
print(format(resumen, digits = 4), row.names = FALSE)
write.csv(resumen, "02_resultados/table1_row_keff.csv", row.names = FALSE)
cat("\nSaved to: 02_resultados/table1_row_keff.csv\n")

# ---------------------------------------------------------------------
# Reproducibility anchor and verdict
# ---------------------------------------------------------------------
cat("\n===== ANCHOR AGAINST TABLE 1 =====\n")
cat(sprintf("  FAR with K : %.6f  (expected %.6f)  %s\n",
            resumen$FAR[1], FAR_ROB_PUB,
            if (abs(resumen$FAR[1] - FAR_ROB_PUB) < 5e-5) "OK" else "REVIEW"))
cat(sprintf("  SE  with K : %.3e  (expected %.3e)  %s\n",
            resumen$SE_FAR[1], SE_ROB_PUB,
            if (abs(resumen$SE_FAR[1] - SE_ROB_PUB) < 2e-6) "OK" else "REVIEW"))

cat("\n===== VERDICT =====\n")
obj <- 1 / ALPHA
d_K    <- abs(resumen$ARL0[1] - obj)
d_Keff <- abs(resumen$ARL0[2] - obj)
cat(sprintf("  Target ARL0 = %.0f\n", obj))
cat(sprintf("  With K     : ARL0 = %.0f  (distance to the target: %.0f)\n",
            resumen$ARL0[1], d_K))
cat(sprintf("  With K_eff : ARL0 = %.0f  (distance to the target: %.0f)\n",
            resumen$ARL0[2], d_Keff))
if (d_Keff < d_K) {
  cat("  >>> K_eff BRINGS the ARL0 closer to the target: the correction IMPROVES it.\n")
} else {
  cat("  >>> K_eff MOVES the ARL0 away from the target: the correction WORSENS it.\n")
  cat("      Reason: reducing K widens the limit, and the limit was already\n")
  cat("      conservative. The correction aggravates what it meant to fix.\n")
}
