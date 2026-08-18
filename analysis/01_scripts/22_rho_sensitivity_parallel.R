# =====================================================================
# 22_rho_sensitivity_parallel.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The sensitivity of the detection advantage to the correlation
#   between variables. It evaluates rho = 0.3, 0.6 and 0.9 under the
#   most demanding scenario of the sensitivity table (6 contaminated
#   Phase 1 batches, mean shift delta = 1.0), reporting for each rho the
#   true positive rate of the classical and of the proposed method with
#   their standard errors, and the advantage of the proposed method.
#
#   The equicorrelation structure fixes Sigma_jj' = rho, Sigma_jj = 1,
#   whose eigenvalues have a closed form:
#       lambda_1 = 1 + (J-1)*rho      (once)
#       lambda_j = 1 - rho            (multiplicity J-1)
#   The AFM weighting is defined through lambda_1, so varying rho varies
#   the spectral separation the mechanism relies on: with J = 4 the
#   ratio lambda_1 / lambda_rest goes from 2.7 at rho = 0.3 to 37.0 at
#   rho = 0.9.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - Everything other than rho stays at the base configuration
#     (K = 30, I = 20, J = 4, h = 0.67, alpha = 0.001, N = 2000
#     replicates), so the observed effect is attributable to rho alone.
#   - The two methods are equalised to the same ARL0 through the
#     empirical 1-alpha quantile of their own statistic over an
#     independent set of N_CAL in-control batches, the same set for both
#     methods within each replicate.
#   - The seeds are deliberately the ones of the power campaign
#     (SEED_BASE*200 / *250 / *300). This makes the rho = 0.6 row
#     reproduce the base row exactly, so the "coincides by
#     construction" footnote is literal, and it makes the three values
#     of rho share common random numbers, which sharpens the comparison
#     between them.
#   - Each replicate seeds its own stream with SEED_BASE*m + idx, so the
#     result does not depend on how work is split between workers and is
#     bit-for-bit identical to the serial loop. clusterSetRNGStream and
#     L'Ecuyer-CMRG are deliberately NOT used: they would break the
#     correspondence with the anchor.
#
# ANCHOR AND ITS EXPECTED VALUES
#   The rho = 0.6 row must reproduce the base row of the sensitivity
#   table: TPR_Hot = 0.2700, TPR_Rob = 0.8629, Advantage = 0.5929,
#   with tolerance 1e-04. rho = 0.6 runs first; if the anchor fails the
#   script stops and rho = 0.3 and 0.9 are not run.
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   The advantage of the proposed method stays above +0.40 at all three
#   values of rho. If it drops below that at some rho the result is
#   reported all the same: it would be a measured limitation, not a
#   failure of the script.
#
# OUTPUT
#   02_resultados/table_3_4c_rho.csv
# =====================================================================

library(parallel)
library(robustT2AFM)
library(MASS)

dir.create("02_resultados", showWarnings = FALSE)

# --- Parameters: identical to the main campaign ----------------------
N_REP     <- 2000
N_CAL     <- 5000
K1        <- 30
K2        <- 100
I         <- 20
J         <- 4
H_MCD     <- 0.67
ALPHA     <- 0.001
OB        <- 6
OR        <- 0.20
OS        <- 4
SHIFT     <- 1.0
SEED_BASE <- 2026
VARS      <- paste0("Var", 1:J)

RHO_GRID  <- c(0.6, 0.3, 0.9)   # 0.6 first: it is the anchor
N_WORKERS <- 8                  # 8 of 12 cores

# --- Anchor and tolerance --------------------------------------------
ANCLA <- c(TPR_Hot = 0.2700, TPR_Rob = 0.8629, Ventaja = 0.5929)
TOL   <- 1e-4

# --- Utilities (identical to the campaign ones) ----------------------
make_equi <- function(p, rho) { S <- matrix(rho, p, p); diag(S) <- 1; S }

t2_manual <- function(Xb, mu, Si, n) {
  d <- colMeans(Xb) - mu
  as.numeric(n * t(d) %*% Si %*% d)
}

se <- function(x) sd(x) / sqrt(length(x))

# --- One full replicate ----------------------------------------------
# 'idx' is the replicate index. The name 'rep' is avoided so as not to
# shadow the base function rep().
one_rep <- function(idx, rho) {

  Sg <- make_equi(J, rho)

  # Phase 1: calibration with 6 contaminated batches
  set.seed(SEED_BASE * 200 + idx)
  sim <- simulate_batch_process(
    K1 = K1, K2 = 1, I = I, J = J, rho = rho, Sigma = Sg,
    outlier_batches_F1 = OB, outlier_rate = OR,
    outlier_shift = OS, prop_ooc_F2 = 0
  )
  f1 <- subset(sim, Phase == "Phase 1")

  calR <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
  calH <- hotelling_classical_calibrate(f1, VARS)

  SiR <- solve(calR$Sw); muR <- calR$mu_r
  SiH <- solve(calH$Sp); muH <- calH$mu_global

  # ARL0 equalisation by empirical quantile, same set for both methods
  set.seed(SEED_BASE * 250 + idx)
  f1c <- lapply(seq_len(N_CAL), function(k) {
    X <- mvrnorm(I, rep(0, J), Sg); colnames(X) <- VARS; X
  })
  qR <- quantile(sapply(f1c, t2_manual, muR, SiR, I), 1 - ALPHA)
  qH <- quantile(sapply(f1c, t2_manual, muH, SiH, I), 1 - ALPHA)

  # Phase 2: batches shifted by delta = 1.0
  set.seed(SEED_BASE * 300 + idx)
  f2 <- lapply(seq_len(K2), function(k) {
    X <- mvrnorm(I, rep(SHIFT, J), Sg); colnames(X) <- VARS; X
  })

  c(
    tprR = mean(sapply(f2, t2_manual, muR, SiR, I) > qR),
    tprH = mean(sapply(f2, t2_manual, muH, SiH, I) > qH)
  )
}

# --- Cluster start-up -------------------------------------------------
cat(sprintf("Starting %d workers out of %d available cores\n",
            N_WORKERS, detectCores()))

cl <- makeCluster(N_WORKERS)
on.exit(stopCluster(cl), add = TRUE)

clusterEvalQ(cl, {
  library(robustT2AFM)
  library(MASS)
  NULL
})

clusterExport(cl, c(
  "N_CAL", "K1", "K2", "I", "J", "H_MCD", "ALPHA",
  "OB", "OR", "OS", "SHIFT", "SEED_BASE", "VARS",
  "make_equi", "t2_manual", "one_rep"
))

# --- Loop over rho ----------------------------------------------------
res <- data.frame()
t_ini <- Sys.time()

for (rho in RHO_GRID) {

  lam1  <- 1 + (J - 1) * rho
  lam_r <- 1 - rho

  cat(sprintf("\n== rho = %.1f  (lambda_1 = %.2f against %.2f, ratio %.1f) ==\n",
              rho, lam1, lam_r, lam1 / lam_r))

  t0  <- Sys.time()
  out <- parLapplyLB(cl, seq_len(N_REP), one_rep, rho = rho)
  dt  <- difftime(Sys.time(), t0, units = "mins")

  tprR <- vapply(out, function(z) z[["tprR"]], numeric(1))
  tprH <- vapply(out, function(z) z[["tprH"]], numeric(1))

  fila <- data.frame(
    config       = sprintf("rho=%.1f", rho),
    rho          = rho,
    lambda1      = round(lam1, 4),
    lambda_rest  = round(lam_r, 4),
    lambda_ratio = round(lam1 / lam_r, 4),
    TPR_Hot      = round(mean(tprH), 4),
    SE_Hot       = round(se(tprH), 4),
    TPR_Rob      = round(mean(tprR), 4),
    SE_Rob       = round(se(tprR), 4),
    Advantage    = round(mean(tprR) - mean(tprH), 4),
    stringsAsFactors = FALSE
  )
  res <- rbind(res, fila)

  cat(sprintf("   TPR_Hot=%.4f  TPR_Rob=%.4f  Advantage=%+.4f   [%.1f min]\n",
              fila$TPR_Hot, fila$TPR_Rob, fila$Advantage, as.numeric(dt)))

  # --- Anchor check, only at rho = 0.6 -------------------------------
  if (abs(rho - 0.6) < 1e-12) {
    d1 <- abs(fila$TPR_Hot - ANCLA[["TPR_Hot"]])
    d2 <- abs(fila$TPR_Rob - ANCLA[["TPR_Rob"]])
    d3 <- abs(fila$Advantage - ANCLA[["Ventaja"]])
    cat("\n   --- Anchor check (base row of the sensitivity table) ---\n")
    cat(sprintf("   TPR_Hot    expected %.4f  obtained %.4f  diff %.6f\n",
                ANCLA[["TPR_Hot"]], fila$TPR_Hot, d1))
    cat(sprintf("   TPR_Rob    expected %.4f  obtained %.4f  diff %.6f\n",
                ANCLA[["TPR_Rob"]], fila$TPR_Rob, d2))
    cat(sprintf("   Advantage  expected %.4f  obtained %.4f  diff %.6f\n",
                ANCLA[["Ventaja"]], fila$Advantage, d3))

    if (max(d1, d2, d3) > TOL) {
      cat("\n   ANCHOR NOT REPRODUCED. The script stops here.\n")
      cat("   rho = 0.3 and rho = 0.9 are not run.\n")
      cat("   Keep this whole output before changing anything.\n")
      write.csv(res, "02_resultados/table_3_4c_rho_ANCHOR_FAILED.csv",
                row.names = FALSE)
      stop("Anchor not reproduced at rho = 0.6.")
    }
    cat("   ANCHOR REPRODUCED. Continuing with rho = 0.3 and rho = 0.9.\n")

    eta <- as.numeric(dt) * (length(RHO_GRID) - 1)
    cat(sprintf("   Estimated time remaining: %.0f min\n", eta))
  }
}

# --- Output -----------------------------------------------------------
res <- res[order(res$rho), ]

cat("\n================ RESULT ================\n")
print(res, row.names = FALSE)

cat(sprintf("\nTotal elapsed: %.1f min\n",
            as.numeric(difftime(Sys.time(), t_ini, units = "mins"))))

cat("\n--- Success criterion fixed in advance: Advantage > 0.40 ---\n")
for (k in seq_len(nrow(res))) {
  cat(sprintf("   %-9s Advantage %+.4f   %s\n",
              res$config[k], res$Advantage[k],
              ifelse(res$Advantage[k] > 0.40, "MEETS IT", "DOES NOT MEET IT (report anyway)")))
}

write.csv(res, "02_resultados/table_3_4c_rho.csv", row.names = FALSE)
cat("\nSaved to 02_resultados/table_3_4c_rho.csv\n")
