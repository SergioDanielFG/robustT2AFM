# =====================================================================
# CENTERED_BOOTSTRAP_v5.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The three resampling alternatives to the analytic control limit, and
#   compares them with it. For each of the four limits (analytic F,
#   bootstrap of the MCD subset, bootstrap of the full batch, and
#   bootstrap of the full batch weighted by w_k) it reports the false
#   alarm rate, its between-replicate standard error, the ARL0 and the
#   total number of false alarms, plus a resolution diagnostic showing
#   which rows were estimated from too few events to support a point
#   ARL0.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - SCALE CORRECTION OF THE MCD RESIDUALS. The residuals of the MCD
#     subset are the m* most central observations of each batch, already
#     centred. Keeping only the central ones shrinks the cloud: their
#     covariance is not Sigma but c^-1 * Sigma, with c the consistency
#     factor of the MCD. Those residuals are evaluated against
#     Si = solve(cal$Sw), and Sw is built from covMcd()$cov, which DOES
#     carry the consistency correction. Without the correction the
#     numerator and the denominator sit on different scales, the
#     bootstrap MCD quantile comes out low, the limit is permissive and
#     the ARL0 is artificially close to the target. Multiplying the
#     residuals by sqrt(c) returns them to the Sigma scale.
#     Consistency factor reference: Hubert et al. (2018).
#   - The bootstrap resamples batches with replacement and, within the
#     drawn batch, observations with replacement, so it reproduces both
#     sources of variation.
#   - The weighted variant reuses the SAME bootstrap statistic as the
#     full-batch variant and only changes the quantile, so the two
#     differ solely in the weighting.
#   - Phase 1 is CONTAMINATED, matching the design of the manuscript.
#   - Each replicate seeds its own stream with SEED_BASE*100 + rep, the
#     same seeds as the analytic campaign, so the analytic row is
#     directly comparable with the published one.
#   - The serial test replicate runs a scale diagnostic and stops the
#     script if the corrected residuals fall outside the band
#     0.90 - 1.10.
#
# ANCHORS AND THEIR EXPECTED VALUES
#   The correction only touches the MCD residuals, so the analytic F row
#   must come out unchanged:
#     FAR analytic F = 0.000625   (must be identical to 1e-09)
#     SE  analytic F = 5.55e-05   (tolerance 2e-06)
#   Scale diagnostic on the test replicate: trace(Sw^-1 cov(res)) / J
#   must be 1.00 for both the corrected MCD residuals and the
#   full-batch control residuals.
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   A resampling limit is declared acceptable if its ARL0 falls inside
#   the band [ARL0_MIN, ARL0_MAX] defined below. The band, the seeds and
#   B are not modified after seeing the result.
#
# OUTPUT
#   02_resultados/table1_rows_bootstrap_v5.csv
# ---------------------------------------------------------------------
# PRE-REGISTERED SUCCESS CRITERION — DO NOT MODIFY
ARL0_MIN <- 900
ARL0_MAX <- 1200
# ---------------------------------------------------------------------
library(parallel)
N_REP  <- 2000
B_BOOT <- 10000
K1 <- 30; K2 <- 100; I <- 20; J <- 4
RHO <- 0.6; H_MCD <- 0.67; ALPHA <- 0.001
OB <- 6; OR <- 0.20; OS <- 4
SEED_BASE <- 2026
VARS <- paste0("Var", 1:J)
make_equi <- function(p, rho) { S <- matrix(rho, p, p); diag(S) <- 1; S }
Sigma_EQ <- make_equi(J, RHO)
FAR_ROB_PUB <- 0.000625
SE_ROB_PUB  <- 5.55e-05
# --- Consistency factor of the trimming ------------------------------
# The covariance of the m* most central observations of a multivariate
# normal equals (1/c) * Sigma, with
#     c = alpha_ret / F_{chi2_{J+2}}( q_{chi2_J}(alpha_ret) )
# Multiplying the residuals by sqrt(c) returns them to the Sigma scale.
m_star    <- round(I * H_MCD)
ALPHA_RET <- m_star / I
C_CONS    <- ALPHA_RET / pchisq(qchisq(ALPHA_RET, J), J + 2)
F_ESCALA  <- sqrt(C_CONS)
cat(sprintf("Effective retention  : %d of %d = %.4f\n", m_star, I, ALPHA_RET))
cat(sprintf("Consistency factor   : c = %.4f\n", C_CONS))
cat(sprintf("Factor on residuals  : sqrt(c) = %.4f\n\n", F_ESCALA))
# --- Weighted quantile -----------------------------------------------
wtd_quantile <- function(x, w, p) {
o <- order(x); x <- x[o]; w <- w[o]
cw <- cumsum(w) / sum(w)
x[which(cw >= p)[1]]
}
# ---------------------------------------------------------------------
# One replicate
# ---------------------------------------------------------------------
una_replica <- function(rep, verbose = FALSE) {
set.seed(SEED_BASE * 100 + rep)
sim <- simulate_batch_process(K1 = K1, K2 = K2, I = I, J = J, rho = RHO,
Sigma = Sigma_EQ,
outlier_batches_F1 = OB, outlier_rate = OR,
outlier_shift = OS, prop_ooc_F2 = 0)
f1 <- subset(sim, Phase == "Phase 1")
f2 <- subset(sim, Phase == "Phase 2")
f1$Batch <- as.character(f1$Batch)
cal <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
Si  <- solve(cal$Sw)
lotes <- lapply(split(f1[, VARS], f1$Batch), as.matrix)
lotes <- lotes[vapply(lotes, function(X) nrow(X) == I, logical(1))]
nm <- names(lotes)
stopifnot(all(nm %in% names(cal$weights)))
stopifnot(all(nm %in% names(cal$mcd_centers)))
w <- as.numeric(cal$weights[match(nm, names(cal$weights))])
stopifnot(!anyNA(w))
ucl_F <- ucl_F_adjusted(cal, I = I, alpha = ALPHA)$UCL
res_full <- lapply(lotes, function(X) sweep(X, 2, colMeans(X)))
# --- Residuals of the MCD subset, RESCALED ------------------------
res_mcd <- lapply(nm, function(k) {
X  <- lotes[[k]]
d  <- mahalanobis(X, cal$mcd_centers[[k]], cal$mcd_covariances[[k]])
Xm <- X[order(d)[seq_len(m_star)], , drop = FALSE]
sweep(Xm, 2, colMeans(Xm)) * F_ESCALA
})
names(res_mcd) <- nm
# --- Scale diagnostic (serial test replicate only) ----------------
if (verbose) {
e_m <- sum(diag(Si %*% cov(do.call(rbind, res_mcd))))  / J
e_f <- sum(diag(Si %*% cov(do.call(rbind, res_full)))) / J
cat(sprintf("  batches used: %d | weight min-max: %.4f - %.4f\n",
length(lotes), min(w), max(w)))
cat(sprintf("  SCALE res_mcd  after correction: %.4f  (target 1.00)\n", e_m))
cat(sprintf("  SCALE res_full (control)       : %.4f  (target 1.00)\n", e_f))
if (e_m < 0.90 || e_m > 1.10)
stop("The scale of res_mcd is still outside the 0.90-1.10 band. STOP.")
cat("  Scale within the band. Continue.\n")
}
# --- Bootstrap ----------------------------------------------------
idx <- sample.int(length(lotes), B_BOOT, replace = TRUE)
t2_mcd  <- numeric(B_BOOT)
t2_full <- numeric(B_BOOT)
w_boot  <- w[idx]
for (b in seq_len(B_BOOT)) {
k  <- idx[b]
Rf <- res_full[[k]]
Rm <- res_mcd[[k]]
rb  <- colMeans(Rf[sample.int(I, I, replace = TRUE), , drop = FALSE])
t2_full[b] <- I * as.numeric(t(rb) %*% Si %*% rb)
rb2 <- colMeans(Rm[sample.int(m_star, I, replace = TRUE), , drop = FALSE])
t2_mcd[b]  <- I * as.numeric(t(rb2) %*% Si %*% rb2)
}
ucl_bMCD  <- quantile(t2_mcd,  1 - ALPHA, names = FALSE)
ucl_bFull <- quantile(t2_full, 1 - ALPHA, names = FALSE)
ucl_bPond <- wtd_quantile(t2_full, w_boot, 1 - ALPHA)
if (verbose) {
cat(sprintf("  UCL: F=%.2f | bMCD=%.2f | bFull=%.2f | bPond=%.2f\n",
ucl_F, ucl_bMCD, ucl_bFull, ucl_bPond))
cat(sprintf("  (without the rescaling bMCD would be around %.2f)\n",
ucl_bMCD * 0.5808))
cat(sprintf("  (reference: 0.999 quantile of chi2_%d = %.2f)\n",
J, qchisq(1 - ALPHA, J)))
}
t2 <- monitor_afm_mcd(f2, cal, VARS)$T2
c(far_F     = mean(t2 > ucl_F),
far_bMCD  = mean(t2 > ucl_bMCD),
far_bFull = mean(t2 > ucl_bFull),
far_bPond = mean(t2 > ucl_bPond),
n_lotes   = length(t2))
}
procesar_bloque <- function(reps) t(vapply(reps, una_replica, numeric(5)))
# ---------------------------------------------------------------------
# SERIAL TEST (stops by itself if the scale was not corrected)
# ---------------------------------------------------------------------
library(robustT2AFM); library(MASS)
cat("=== Serial check: replicate 1 (scale corrected) ===\n")
print(una_replica(1, verbose = TRUE))
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
"B_BOOT", "wtd_quantile", "una_replica",
"m_star", "F_ESCALA"))
res <- parLapply(cl, bloques, procesar_bloque)
stopCluster(cl)
M <- do.call(rbind, res)
stopifnot(nrow(M) == N_REP, all(M[, "n_lotes"] == K2))
cat(sprintf("Elapsed: %.1f min | %d workers\n\n",
as.numeric(difftime(Sys.time(), t0, units = "mins")), n_workers))
# ---------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------
cols <- c("far_F", "far_bMCD", "far_bFull", "far_bPond")
etiq <- c("F analitico (actual)",
"Bootstrap subconjunto MCD (centrado, escala corregida)",
"Bootstrap lote completo (centrado)",
"Bootstrap ponderado por w_k (centrado)")
resumen <- data.frame(
Limit     = etiq,
FAR       = vapply(cols, function(c) mean(M[, c]), numeric(1)),
SE_FAR    = vapply(cols, function(c) sd(M[, c]) / sqrt(N_REP), numeric(1)),
ARL0      = vapply(cols, function(c) 1 / mean(M[, c]), numeric(1)),
FA_total  = vapply(cols, function(c) round(sum(M[, c]) * K2), numeric(1)),
n_batches = rep(N_REP * K2, 4),
row.names = NULL, stringsAsFactors = FALSE
)
cat("===== COMPARISON OF LIMITS =====\n")
print(format(resumen, digits = 4), row.names = FALSE)
write.csv(resumen, "02_resultados/table1_rows_bootstrap_v5.csv", row.names = FALSE)
cat("\nSaved to: 02_resultados/table1_rows_bootstrap_v5.csv\n")
cat("(the earlier version file is NOT overwritten)\n")
# ---------------------------------------------------------------------
# Resolution of each row
# ---------------------------------------------------------------------
cat("\n===== RESOLUTION OF EACH ROW =====\n")
cat("  When SE/FAR approaches 1, the FAR was estimated from one or two\n")
cat("  events and is NOT distinguishable from zero. Those rows should be\n")
cat("  reported as a FAR with its SE, avoiding a point ARL0.\n\n")
for (i in seq_len(nrow(resumen))) {
r <- resumen$SE_FAR[i] / resumen$FAR[i]
cat(sprintf("  %-55s FA=%-5.0f SE/FAR=%.3f %s\n",
resumen$Limit[i], resumen$FA_total[i], r,
if (r > 0.5) "<-- no resolution" else ""))
}
# ---------------------------------------------------------------------
# ANCHOR: the analytic F row must not have moved
# ---------------------------------------------------------------------
cat("\n===== ANCHOR AGAINST TABLE 1 =====\n")
cat("  The correction only touches res_mcd. The analytic F row must come\n")
cat("  out IDENTICAL to the earlier one. If it changes, something else moved.\n\n")
cat(sprintf("  FAR analytic F : %.6f  (expected %.6f)  %s\n",
resumen$FAR[1], FAR_ROB_PUB,
if (abs(resumen$FAR[1] - FAR_ROB_PUB) < 1e-9) "IDENTICAL" else "REVIEW"))
cat(sprintf("  SE  analytic F : %.3e  (expected %.3e)  %s\n",
resumen$SE_FAR[1], SE_ROB_PUB,
if (abs(resumen$SE_FAR[1] - SE_ROB_PUB) < 2e-6) "OK" else "REVIEW"))
cat("\n===== VERDICT (PRE-REGISTERED criterion, not modified) =====\n")
cat(sprintf("  Band fixed in advance: ARL0 within [%d, %d]\n", ARL0_MIN, ARL0_MAX))
for (i in 2:4) {
a <- resumen$ARL0[i]
cumple <- is.finite(a) && a >= ARL0_MIN && a <= ARL0_MAX
cat(sprintf("  %-55s ARL0=%-10s %s\n", resumen$Limit[i],
if (is.finite(a)) sprintf("%.0f", a) else "Inf",
if (cumple) "MEETS IT" else "does not meet it"))
}
cat("\n  Without the scale correction the bootstrap MCD limit met the band.\n")
cat("  If it now moves away from the target, that agreement came from the\n")
cat("  scale error and not from a better calibration.\n\n")
