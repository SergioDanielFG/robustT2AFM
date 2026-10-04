# =====================================================================
# SCRIPT_M2M3_diagnostics.R
# ---------------------------------------------------------------------
# Two short checks on results already obtained. Nothing new is
# simulated for the tables of the paper.
#
# PART 1  Paired standard errors for M2 (Reviewer 1, comments 10 and 11)
#   M2 reported each chart with its own standard error. The charts share
#   the same data in every replicate, so the right test of a difference
#   is paired: NEW - other, replicate by replicate. Same reading rule as
#   M3, fixed in advance: a difference is declared only when |z| > 3.
#   Uses the saved replicates in 03_T2_crudos/M2/.
#   Check: the means recomputed from the files must reproduce
#   table_M2_all.csv (tolerance 5e-04) for all 27 configurations.
#
# PART 2  Why the pooled charts collapse in block W of M3
#   Block W: 6 of the 30 Phase 1 batches shifted as a whole by 1.5 or 3
#   sd in every variable. Pooled RMCD fell to TPR 0.12-0.16 in direction
#   D1. The Phase 1 data of every replicate are regenerated with the
#   same seeds and the same call as M3 (SEED_BASE*1400 + rep), and the
#   estimators are fitted in the same order, so the fits are the ones
#   that produced the M3 results. For each replicate it records:
#     - share of the shifted observations kept by the raw MCD subset of
#       the pooled RMCD, and by its reweighting step
#     - offset of each reference center along D1, in sd (the mean of its
#       four coordinates; every variable has variance 1)
#     - first eigenvalue of each covariance (true value 2.8)
#   All 2000 replicates of each shift, in parallel.
#   Check: the classical center (average of the batch means) must be
#   offset by 0.2 x shift (6 of 30 batches), i.e. 0.30 and 0.60.
#
# OUTPUT (05_revision_R1/02_resultados/)
#   table_M2_paired.csv, table_M3_W_diagnostic.csv
#   Run from C:/temp_paper. Expected time: a few minutes.
# =====================================================================
library(parallel)
library(robustT2AFM)
library(robustbase)
library(rrcov)
library(MASS)
stopifnot(packageVersion("robustT2AFM") >= "0.3.0")

DIR_OUT <- "05_revision_R1/02_resultados"
DIR_M2  <- "05_revision_R1/03_T2_crudos/M2"
N_REP <- 2000; BLOCK <- 200; SEED_BASE <- 2026; OR <- 0.20; N_WORKERS <- 10
CHARTS <- c("a_classical", "b_MCD_unif", "c_AFM_cls", "d_V7", "e_NEW", "f_MRCD", "g_RMCD_pooled")

# =====================================================================
# PART 1: paired standard errors for M2
# =====================================================================
cat("=== PART 1: M2 paired differences ===\n")
m2 <- read.csv(file.path(DIR_OUT, "table_M2_all.csv"))
ok1 <- TRUE
PAIR2 <- do.call(rbind, lapply(m2$id, function(i) {
  f <- file.path(DIR_M2, sprintf("M2_%s_%04d.rds", i, seq(1, N_REP, BLOCK)))
  if (!all(file.exists(f))) stop("Missing M2 replicate files for ", i)
  M <- do.call(rbind, lapply(f, readRDS)); stopifnot(nrow(M) == N_REP)
  got <- colMeans(M[, CHARTS]); ref <- unlist(m2[m2$id == i, CHARTS])
  if (any(abs(got - ref) >= 5e-4)) { ok1 <<- FALSE; cat("  MISMATCH with table_M2_all.csv in", i, "\n") }
  r <- m2[m2$id == i, ]
  do.call(rbind, lapply(c("a_classical", "d_V7", "f_MRCD", "g_RMCD_pooled"), function(o) {
    dd <- M[, "e_NEW"] - M[, o]; se <- sd(dd) / sqrt(length(dd)); z <- if (se > 0) mean(dd) / se else NA
    data.frame(id = i, block = r$block, dir = r$dir, versus = o,
               NEW = round(mean(M[, "e_NEW"]), 4), other = round(mean(M[, o]), 4),
               diff = round(mean(dd), 4), SE_paired = round(se, 4), z = round(z, 1),
               verdict = if (is.na(z) || abs(z) <= 3) "no difference" else if (z > 0) "NEW better" else "NEW worse")
  }))
}))
cat(if (ok1) "  Check OK: the saved replicates reproduce table_M2_all.csv.\n\n"
    else "  CHECK FAILED: the files do not reproduce table_M2_all.csv. Do NOT use Part 1.\n\n")
write.csv(PAIR2, file.path(DIR_OUT, "table_M2_paired.csv"), row.names = FALSE)
print(PAIR2[, c("id", "versus", "NEW", "other", "diff", "z", "verdict")], row.names = FALSE)
cat("\n  Verdicts of NEW against each chart:\n")
print(table(PAIR2$versus, PAIR2$verdict))

# =====================================================================
# PART 2: pooled estimators in block W of M3
# =====================================================================
cat("\n=== PART 2: block W of M3, why the pooled charts collapse ===\n")
Sg <- matrix(0.6, 4, 4); diag(Sg) <- 1
diag_rep <- function(sc, rep) {
  vv   <- paste0("Var", 1:4)
  lam1 <- function(S) max(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
  set.seed(SEED_BASE * 1400 + rep)                         # same seed and call as M3, block W
  sim <- simulate_batch_process(K1 = 30, K2 = 1, I = 20, J = 4, rho = 0.6, Sigma = Sg,
                                outlier_batches_F1 = 0, outlier_rate = OR,
                                outlier_shift = 0, prop_contam_F1 = 0.2,
                                shift_contam = sc, prop_ooc_F2 = 0)
  f1 <- subset(sim, Phase == "Phase 1")
  shifted <- f1$ContaminationType == "Shifted"
  cal <- suppressMessages(calibrate_afm_mcd(f1, vv, mcd_alpha = 0.67))   # same order of fits as M3
  X   <- as.matrix(f1[, vv])
  mr  <- suppressWarnings(rrcov::CovMrcd(X, alpha = 0.67))
  rm_ <- robustbase::covMcd(X, alpha = 0.67)
  stopifnot(!is.null(rm_$best), !is.null(rm_$mcd.wt))
  C      <- do.call(rbind, cal$mcd_centers)[, vv, drop = FALSE]
  mu_cls <- colMeans(do.call(rbind, lapply(unique(f1$Batch), function(bb) colMeans(f1[f1$Batch == bb, vv]))))
  c(shift = sc,
    RMCD_raw_keeps_shifted   = sum(shifted[rm_$best]) / sum(shifted),
    RMCD_rew_keeps_shifted   = sum(rm_$mcd.wt[shifted] > 0) / sum(shifted),
    RMCD_rew_share_shifted   = sum(rm_$mcd.wt[shifted] > 0) / sum(rm_$mcd.wt > 0),
    offset_classical         = mean(mu_cls),
    offset_V7_mean_centers   = mean(colMeans(C)),
    offset_NEW_MCD_centers   = mean(cal$mu_r),
    offset_RMCD_raw          = mean(rm_$raw.center),
    offset_RMCD_reweighted   = mean(rm_$center),
    offset_MRCD              = mean(rrcov::getCenter(mr)),
    lambda1_NEW_Sw           = lam1(cal$Sw),
    lambda1_RMCD_reweighted  = lam1(rm_$cov),
    lambda1_MRCD             = lam1(rrcov::getCov(mr)))
}
jobs <- split(expand.grid(sc = c(1.5, 3), rep = seq_len(N_REP)), rep(seq_len(N_WORKERS * 4), length.out = 2 * N_REP))
t0 <- Sys.time()
clu <- makeCluster(N_WORKERS)
clusterEvalQ(clu, { library(robustT2AFM); library(robustbase); library(rrcov); library(MASS); NULL })
clusterExport(clu, c("diag_rep", "Sg", "SEED_BASE", "OR"))
res <- parLapplyLB(clu, jobs, function(jb) t(mapply(diag_rep, jb$sc, jb$rep)))
stopCluster(clu)
res <- as.data.frame(do.call(rbind, res))
cat(sprintf("  %d replicates in %.1f min\n", nrow(res), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

WD <- aggregate(. ~ shift, data = res, FUN = mean)
WD$lambda1_true <- 2.8
WD[, -1] <- round(WD[, -1], 3)
write.csv(WD, file.path(DIR_OUT, "table_M3_W_diagnostic.csv"), row.names = FALSE)
print(t(WD))

ok2 <- all(abs(WD$offset_classical - 0.2 * WD$shift) < 0.02)
cat(if (ok2) "\n  Check OK: the classical center is offset by 0.2 x shift, as it must be.\n"
    else "\n  CHECK FAILED: the classical offset is not 0.2 x shift; the data are not those of M3.\n")
cat(if (ok1 && ok2) "\nALL CHECKS OK.\n" else "\nCHECK FAILURES: do NOT use these results.\n")
