# =====================================================================
# 20_arl0_clean_vs_contaminated_v2.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   Whether the false-alarm rate of the proposed analytic limit changes
#   when Phase 1 is contaminated. Ten cells are evaluated: five
#   disjoint seed sets crossed with Phase 1 clean and Phase 1
#   contaminated with six batches. For each cell it reports the false
#   alarms out of 200 000 in-control batches, the ARL0 and its exact
#   Poisson 95% CI; it then pools by condition, tests homogeneity
#   across the ten cells and clean against contaminated, and reports
#   the spread across seed sets at fixed K.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - The five seed bases are separated by 20 000, far above N_REP.
#     Since the seed of a replicate is SEED_BASE + r with r up to
#     N_REP, bases closer together than N_REP would make the seed
#     ranges overlap and the campaigns would silently be the same one.
#     The script CHECKS that separation before starting and stops if it
#     is not met.
#   - K is fixed at 30 in every cell, so the spread observed across
#     seed sets is attributable to the estimation procedure and not to
#     the number of calibration batches.
#   - The UCL is the deterministic analytic limit, identical in all ten
#     cells, so the comparison isolates the effect of contamination on
#     the statistic and not on the limit.
#
# ANCHORS AND THEIR EXPECTED VALUES
#   A1  Deterministic UCL = 19.6929  (tolerance 5e-04)
#   A2  The cell (seed base 2026, CLEAN) must give EXACTLY 125 false
#       alarms, reproducing the K = 30 row of arl0_convergence_K_v3.csv
#
# SUCCESS CRITERIA FIXED IN ADVANCE
#   (1) homogeneity across the ten cells is not rejected (p > 0.05)
#   (2) the ratio ARL0_clean / ARL0_contaminated lies in [0.85, 1.18]
#   Both       -> INSENSITIVE, published as a property of the method.
#   Only (1)   -> NOT CONCLUSIVE, modest wording.
#   Neither    -> SENSITIVE, modest wording.
#
# OUTPUT
#   02_resultados/arl0_clean_vs_contaminated_v2.csv
# =====================================================================

library(robustT2AFM)
library(parallel)

N_REP      <- 2000
K          <- 30
I          <- 20
J          <- 4
H_MCD      <- 0.67
ALPHA      <- 0.001
RHO        <- 0.6
N_F2       <- 100
OB         <- 6
OR         <- 0.20
OS         <- 4
SEED_BASES <- c(2026, 20000, 40000, 60000, 80000)
N_CORES    <- 10
DIR_OUT    <- "02_resultados"

VARS <- paste0("Var", 1:J)

# ---- The seed bases must be farther apart than N_REP, or the seed
# ---- ranges of the campaigns overlap.
sep_min <- min(diff(sort(SEED_BASES)))
cat("=== Check that the seed sets are independent ===\n")
cat(sprintf("  minimum separation between bases = %d | required > %d\n",
            sep_min, N_REP))
if (sep_min <= N_REP) {
  stop("The seed bases overlap. Separate them by more than N_REP.")
}
cat("  OK. The campaigns use disjoint seed sets.\n\n")

ucl_analitico <- function(K, I, J, h, alpha) {
  m  <- round(I * h)
  d2 <- K * m - K - J + 1
  (J * (K + 1) * (m - 1) / d2) * qf(1 - alpha, J, d2)
}
UCL <- ucl_analitico(K, I, J, H_MCD, ALPHA)

cat("=== Anchor A1: deterministic UCL ===\n")
cat(sprintf("  UCL = %.4f   | expected 19.6929\n", UCL))
stopifnot(abs(UCL - 19.6929) < 5e-4)
cat("  OK.\n\n")

una_replica <- function(r, K, I, J, N_F2, H_MCD, RHO, VARS, UCL,
                        SEED_BASE, OB_local, OR, OS) {
  sim <- simulate_batch_process(
    K1 = K, K2 = N_F2, I = I, J = J, rho = RHO,
    outlier_batches_F1 = OB_local,
    outlier_rate       = OR,
    outlier_shift      = OS,
    prop_contam_F1     = 0,
    prop_ooc_F2        = 0,
    seed = SEED_BASE + r
  )
  cal <- suppressMessages(
    calibrate_afm_mcd(subset(sim, Phase == "Phase 1"), VARS, mcd_alpha = H_MCD))
  mon <- monitor_afm_mcd(subset(sim, Phase == "Phase 2"), cal, VARS)
  sum(mon$T2 > UCL)
}

condiciones <- expand.grid(seed_base = SEED_BASES,
                           contam = c(0, OB), KEEP.OUT.ATTRS = FALSE)

cl <- makeCluster(N_CORES)
clusterEvalQ(cl, library(robustT2AFM))

resultados <- data.frame()
for (i in seq_len(nrow(condiciones))) {
  sb  <- condiciones$seed_base[i]
  obl <- condiciones$contam[i]
  etq <- if (obl == 0) "limpia" else "contaminada"
  cat(sprintf("--- base %d | Phase 1 %s ---\n", sb,
              if (obl == 0) "clean" else "contaminated"))
  t0 <- Sys.time()
  alarmas <- parSapply(cl, seq_len(N_REP), una_replica,
    K = K, I = I, J = J, N_F2 = N_F2, H_MCD = H_MCD, RHO = RHO,
    VARS = VARS, UCL = UCL, SEED_BASE = sb, OB_local = obl, OR = OR, OS = OS)
  lt <- N_REP * N_F2; fa <- sum(alarmas)
  a0 <- if (fa > 0) lt / fa else NA_real_
  lo <- if (fa > 0) lt / (qchisq(0.975, 2*fa + 2)/2) else NA_real_
  hi <- if (fa > 0) lt / (qchisq(0.025, 2*fa)/2)     else NA_real_
  cat(sprintf("   %d of %d | ARL0 = %.0f  95%% CI [%.0f, %.0f] | %.1f min\n\n",
      fa, lt, a0, lo, hi, as.numeric(difftime(Sys.time(), t0, units="mins"))))
  resultados <- rbind(resultados, data.frame(
    seed_base = sb, condition = etq, contaminated_batches = obl, UCL = UCL,
    total_batches = lt, false_alarms = fa, FAR = fa/lt, ARL0 = a0,
    ARL0_CI95_low = lo, ARL0_CI95_high = hi))
}
stopCluster(cl)

cat("=== Anchor A2 ===\n")
anc <- resultados[resultados$seed_base == 2026 & resultados$condition == "limpia", ]
cat(sprintf("  false alarms = %d (expected 125) %s | ARL0 = %.0f\n\n",
            anc$false_alarms,
            if (anc$false_alarms == 125) "OK" else "*** FAILS ***", anc$ARL0))
if (anc$false_alarms != 125)
  stop("Anchor A2 failed. STOP: the comparisons are not valid.")

cat("=== FULL TABLE ===\n")
print(resultados[, c("seed_base","condition","false_alarms","ARL0",
                     "ARL0_CI95_low","ARL0_CI95_high")], row.names = FALSE)

lim <- resultados[resultados$condition == "limpia", ]
con <- resultados[resultados$condition == "contaminada", ]
fl <- sum(lim$false_alarms)/sum(lim$total_batches)
fc <- sum(con$false_alarms)/sum(con$total_batches)
cociente <- (1/fl)/(1/fc)

cat("\n=== POOLED ===\n")
cat(sprintf("  CLEAN        : %d of %d -> ARL0 = %.0f\n",
            sum(lim$false_alarms), sum(lim$total_batches), 1/fl))
cat(sprintf("  CONTAMINATED : %d of %d -> ARL0 = %.0f\n",
            sum(con$false_alarms), sum(con$total_batches), 1/fc))
cat(sprintf("  ratio = %.4f\n", cociente))

cat("\n=== SPREAD ACROSS SEED SETS, K FIXED AT 30 (Phase 1 clean) ===\n")
cat(sprintf("  minimum ARL0 = %.0f | maximum = %.0f | range = %.0f\n",
            min(lim$ARL0), max(lim$ARL0), max(lim$ARL0) - min(lim$ARL0)))
cat("  Compare with the range obtained varying K from 30 to 500: 1504 to 1802.\n")

cat("\n=== TESTS ===\n")
tabN <- cbind(resultados$false_alarms,
              resultados$total_batches - resultados$false_alarms)
hN <- suppressWarnings(chisq.test(tabN))
cat(sprintf("  (1) homogeneity of the %d cells: X2 = %.3f, p = %.4f\n",
            nrow(resultados), hN$statistic, hN$p.value))
pt <- prop.test(c(sum(lim$false_alarms), sum(con$false_alarms)),
                c(sum(lim$total_batches), sum(con$total_batches)), correct = FALSE)
cat(sprintf("      clean vs contaminated: p = %.4f\n", pt$p.value))
cat(sprintf("  (2) ratio = %.4f within [0.85, 1.18]: %s\n",
            cociente, if (cociente >= 0.85 && cociente <= 1.18) "YES" else "NO"))

c1 <- hN$p.value > 0.05; c2 <- cociente >= 0.85 && cociente <= 1.18
cat("\n=========================================================\n")
if (c1 && c2) {
  cat(" VERDICT: INSENSITIVE. Published as a property of the method.\n")
} else if (c1 && !c2) {
  cat(" VERDICT: NOT CONCLUSIVE. Modest wording.\n")
} else {
  cat(" VERDICT: SENSITIVE. Modest wording.\n")
}
cat("=========================================================\n")

if (!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)
write.csv(resultados, file.path(DIR_OUT, "arl0_clean_vs_contaminated_v2.csv"),
          row.names = FALSE)
cat("\nSaved. Move the v1 file to 99_obsoleto.\n")
cat("=== END. Copy the whole output. ===\n\n")
