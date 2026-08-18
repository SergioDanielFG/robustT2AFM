# =====================================================================
# ZTEST_TABLE1_FAITHFUL_v2.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The two analytic rows of Table 1 (false-alarm calibration under
#   control): false alarm rate (FAR), its standard error, the ARL0 and
#   the false alarm count, for the proposed robust F limit with m* and
#   for the classical Montgomery limit. It also runs the two-proportion
#   test and Fisher's exact test between both limits.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - The FAR of EVERY replicate is stored, not just the total count.
#     The standard error is therefore the between-replicate one,
#     sd(far)/sqrt(N_REP), which is the criterion under which the two
#     published standard errors of Table 1 were obtained.
#   - Phase 1 is CONTAMINATED (6 batches, 20% of observations shifted
#     4 SD), matching the design of the manuscript.
#   - Each replicate seeds its own stream with SEED_BASE*100 + rep, so
#     the result does not depend on how replicates are split between
#     workers.
#   - The design is balanced (K2 constant across replicates), so
#     mean(fa_rep/K2) = sum(fa_rep)/(N_REP*K2); the point FAR is the
#     same under either criterion and only the SE differs.
#   - Each method keeps its OWN analytic limit; the two are not brought
#     to a common false-alarm rate here.
#
# ANCHORS AND THEIR EXPECTED VALUES
#   FAR robust    = 0.000625   (tolerance 5e-05)
#   SE  robust    = 5.55e-05   (tolerance 2e-06)   <-- decisive anchor
#   FAR classical = 0.000730   (tolerance 5e-05)
#   SE  classical = 6.35e-05   (tolerance 2e-06)
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   The two SE anchors must report OK. That establishes that the
#   between-replicate criterion is the one of the original campaign and
#   that further rows of Table 1 can be computed the same way.
#
# OUTPUT
#   02_resultados/table1_ztest_faithful.csv
#   02_resultados/far_per_replicate.csv
# =====================================================================
library(parallel)

N_REP     <- 2000
K1 <- 30; K2 <- 100; I <- 20; J <- 4
RHO <- 0.6; H_MCD <- 0.67; ALPHA <- 0.001
OB <- 6; OR <- 0.20; OS <- 4        # Phase 1 CONTAMINATED, as in the campaign
SEED_BASE <- 2026
VARS <- paste0("Var", 1:J)

make_equi <- function(p, rho) { S <- matrix(rho, p, p); diag(S) <- 1; S }
Sigma_EQ <- make_equi(J, RHO)

# --- Published values, used by the anchors ---------------------------
FAR_ROB_PUB <- 0.000625
SE_ROB_PUB  <- 5.55e-05
FAR_CLS_PUB <- 0.000730
SE_CLS_PUB  <- 6.35e-05

# ---------------------------------------------------------------------
# One replicate: returns the FAR of that replicate for each method
# ---------------------------------------------------------------------
una_replica <- function(rep) {
  set.seed(SEED_BASE * 100 + rep)
  sim <- simulate_batch_process(K1 = K1, K2 = K2, I = I, J = J, rho = RHO,
                                Sigma = Sigma_EQ,
                                outlier_batches_F1 = OB, outlier_rate = OR,
                                outlier_shift = OS, prop_ooc_F2 = 0)
  f1 <- subset(sim, Phase == "Phase 1")
  f2 <- subset(sim, Phase == "Phase 2")

  calR <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
  calH <- hotelling_classical_calibrate(f1, VARS)

  uR <- ucl_F_adjusted(calR, I = I, alpha = ALPHA)$UCL
  uH <- hotelling_classical_ucl(K = calH$n_batches, I = I, J = J,
                                alpha = ALPHA, phase = "II")$UCL

  t2R <- monitor_afm_mcd(f2, calR, VARS)$T2
  t2H <- hotelling_classical_monitor(f2, calH, VARS)$T2

  c(far_rob = mean(t2R > uR),
    far_cls = mean(t2H > uH),
    n_lotes = length(t2R))
}

procesar_bloque <- function(reps) {
  t(vapply(reps, una_replica, numeric(3)))
}

# ---------------------------------------------------------------------
# SERIAL CHECK: one replicate before opening the cluster, so that any
# error surfaces with its real message instead of "8 nodes produced
# errors"
# ---------------------------------------------------------------------
library(robustT2AFM); library(MASS)
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
                    "una_replica"))
resultados <- parLapply(cl, bloques, procesar_bloque)
stopCluster(cl)

M <- do.call(rbind, resultados)
far_rob <- M[, "far_rob"]
far_cls <- M[, "far_cls"]

# Integrity checks before computing anything
stopifnot(length(far_rob) == N_REP)
stopifnot(all(M[, "n_lotes"] == K2))   # balanced design: mean(far) = total/n

cat(sprintf("Elapsed: %.1f min | %d workers\n\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), n_workers))

# ---------------------------------------------------------------------
# Summary under the BETWEEN-REPLICATE criterion
# ---------------------------------------------------------------------
n_eval <- N_REP * K2
fa_rob <- round(sum(far_rob) * K2)
fa_cls <- round(sum(far_cls) * K2)

resumen <- data.frame(
  Limit     = c("Robusto F (m*)", "Clasico (Montgomery)"),
  FAR       = c(mean(far_rob), mean(far_cls)),
  SE_FAR    = c(sd(far_rob) / sqrt(N_REP), sd(far_cls) / sqrt(N_REP)),
  ARL0      = c(1 / mean(far_rob), 1 / mean(far_cls)),
  FA_total  = c(fa_rob, fa_cls),
  n_batches = c(n_eval, n_eval),
  stringsAsFactors = FALSE
)

cat("===== TABLE 1 (between-replicate criterion) =====\n")
print(format(resumen, digits = 4), row.names = FALSE)
write.csv(resumen, "02_resultados/table1_ztest_faithful.csv", row.names = FALSE)
write.csv(data.frame(replicate = seq_len(N_REP),
                     far_rob = far_rob,
                     far_cls = far_cls),
          "02_resultados/far_per_replicate.csv", row.names = FALSE)

cat("\nSaved to: 02_resultados/table1_ztest_faithful.csv\n")


# ---------------------------------------------------------------------
# Overdispersion diagnostic: explains why the classical standard error
# exceeds the binomial one and the robust one does not
# ---------------------------------------------------------------------
cat("\n===== DIAGNOSTIC: variance / mean of the false alarms =====\n")
cat("  (=1 -> Poisson behaviour; >1 -> the calibration varies between replicates)\n")
for (nm in c("rob", "cls")) {
  fa_rep <- (if (nm == "rob") far_rob else far_cls) * K2
  se_bin <- sqrt(mean(fa_rep / K2) * (1 - mean(fa_rep / K2)) / n_eval)
  se_ent <- sd(fa_rep / K2) / sqrt(N_REP)
  cat(sprintf("  %s : mean=%.5f  var=%.5f  var/mean=%.4f | binomial SE=%.3e  between-replicate SE=%.3e  ratio=%.4f\n",
              toupper(nm), mean(fa_rep), var(fa_rep), var(fa_rep) / mean(fa_rep),
              se_bin, se_ent, se_ent / se_bin))
}

# ---------------------------------------------------------------------
# Formal two-proportion test, on the counts
# ---------------------------------------------------------------------
cat("\n===== FORMAL TWO-PROPORTION TEST =====\n")
pt <- prop.test(x = c(fa_rob, fa_cls), n = c(n_eval, n_eval), correct = FALSE)
ft <- fisher.test(matrix(c(fa_rob, n_eval - fa_rob,
                           fa_cls, n_eval - fa_cls), nrow = 2))
cat(sprintf("  z-test: X-squared=%.4f, p-value=%.4f\n", pt$statistic, pt$p.value))
cat(sprintf("  Fisher: p-value=%.4f\n", ft$p.value))

# ---------------------------------------------------------------------
# ANCHORS
# ---------------------------------------------------------------------
chk <- function(obt, esp, tol) if (abs(obt - esp) < tol) "OK" else "REVIEW"

cat("\n===== ANCHORS AGAINST THE PUBLISHED TABLE 1 =====\n")
cat(sprintf("  FAR robust    : %.6f  (expected %.6f)  %s\n",
            resumen$FAR[1], FAR_ROB_PUB, chk(resumen$FAR[1], FAR_ROB_PUB, 5e-5)))
cat(sprintf("  SE  robust    : %.3e  (expected %.3e)  %s   <-- decisive anchor\n",
            resumen$SE_FAR[1], SE_ROB_PUB, chk(resumen$SE_FAR[1], SE_ROB_PUB, 2e-6)))
cat(sprintf("  FAR classical : %.6f  (expected %.6f)  %s\n",
            resumen$FAR[2], FAR_CLS_PUB, chk(resumen$FAR[2], FAR_CLS_PUB, 5e-5)))
cat(sprintf("  SE  classical : %.3e  (expected %.3e)  %s\n",
            resumen$SE_FAR[2], SE_CLS_PUB, chk(resumen$SE_FAR[2], SE_CLS_PUB, 2e-6)))

cat("\nIf the two SE anchors report OK, it is established that the\n")
cat("between-replicate criterion is the one of the original campaign,\n")
cat("and the further rows of Table 1 can be computed the same way.\n")
