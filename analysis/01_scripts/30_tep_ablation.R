# =====================================================================
# 30_tep_ablation.R
#
# WHAT IT COMPUTES
#   The four rows of the ablation block of Table 10: the 2 x 2
#   decomposition of the method on the Tennessee Eastman application,
#   with the same Phase 1 (24 healthy + 6 IDV 7) and Phase 2 (10 healthy
#   + 20 IDV 7) as PHASE1_PIPELINE_STEP30.R.
#
#     (a) classical covariance, uniform weights   -> classical T2
#     (b) MCD covariance,       uniform weights   -> MCD only
#     (c) classical covariance, AFM weights       -> AFM only
#     (d) MCD covariance,       AFM weights       -> AFM-MCD (proposed)
#
#   Each variant uses its own analytic limit: (a) and (c) the classical
#   Phase 2 limit with K, I, J; (b) and (d) the adjusted limit with m*.
#   The reference center is the average of the MCD centers for (b) and
#   (d) and the average of the batch means for (a) and (c), exactly as
#   in SCRIPT_FINAL_ABLATIONS_2x2.R on simulated data.
#
# DEPENDENCIES
#   Sources PHASE1_PIPELINE_STEP30.R, which leaves f1, f2, cal_rob,
#   cal_cls, tr, tc, ef, ucl_rob and ucl_cls in memory. Nothing from the
#   pipeline is recomputed: variants (a) and (d) must reproduce tc and
#   tr, which is the internal check of this script.
#
# ANCHORS (tep_ablation.csv, values of the manuscript)
#   (a) UCL 19.46 | detected  2/20 | false alarms 0/10 | T2 median 11.5
#   (b) UCL 19.69 | detected 20/20 | false alarms 0/10 | T2 median 56.5
#   (c) UCL 19.46 | detected 20/20 | false alarms 0/10 | T2 median 54.4
#   (d) UCL 19.69 | detected 20/20 | false alarms 0/10 | T2 median 88.7
#   T2 median is the median of the statistic over the 20 faulty batches.
#   covMcd() explores subsets at random, so the medians of (b) and (d)
#   may move by a few tenths between runs; the counts must not.

# OUTPUT
#   02_resultados/tep_ablation.csv
# =====================================================================
library(robustT2AFM)
if (!exists("cal_rob") || !exists("tr") || !exists("tc"))
  source("01_scripts/PHASE1_PIPELINE_STEP30.R", echo = FALSE)
stopifnot(exists("f1"), exists("f2"), exists("cal_rob"), exists("cal_cls"),
          exists("tr"), exists("tc"), exists("ef"),
          exists("ucl_rob"), exists("ucl_cls"), exists("VARS"), exists("I_LOTE"))

DIR_OUT <- "02_resultados"
dir.create(DIR_OUT, showWarnings = FALSE)

# --- 1. Building blocks from Phase 1 ---------------------------------
batches1 <- unique(f1$Batch)
S_mcd  <- cal_rob$mcd_covariances
mu_mcd <- cal_rob$mu_r
S_cls  <- lapply(batches1, function(b) cov(f1[f1$Batch == b, VARS]))
mu_cls <- colMeans(do.call(rbind, lapply(batches1, function(b)
  colMeans(f1[f1$Batch == b, VARS]))))

lam1_cls  <- sapply(S_cls, function(S)
  max(eigen(S, symmetric = TRUE, only.values = TRUE)$values))
w_afm_cls <- (1 / lam1_cls) / sum(1 / lam1_cls)
unif      <- rep(1 / length(batches1), length(batches1))
comb      <- function(Sl, w) Reduce(`+`, Map(function(S, wi) wi * S, Sl, w))

refs <- list(
  a = list(label = "(a) Classical", cov = "Classical", wts = "Uniform",
           mu = mu_cls, Si = solve(comb(S_cls, unif)), ucl = ucl_cls),
  b = list(label = "(b) MCD only",  cov = "MCD",       wts = "Uniform",
           mu = mu_mcd, Si = solve(comb(S_mcd, unif)), ucl = ucl_rob),
  c = list(label = "(c) AFM only",  cov = "Classical", wts = "AFM",
           mu = mu_cls, Si = solve(comb(S_cls, w_afm_cls)), ucl = ucl_cls),
  d = list(label = "(d) AFM-MCD",   cov = "MCD",       wts = "AFM",
           mu = mu_mcd, Si = solve(cal_rob$Sw), ucl = ucl_rob)
)

# --- 2. Phase 2 statistic for each variant ---------------------------
t2_manual <- function(Xb, mu, Si, n) {
  d <- colMeans(Xb) - mu
  as.numeric(n * t(d) %*% Si %*% d)
}
batches2 <- unique(f2$Batch)
stopifnot(identical(ef, grepl("fallo", batches2)))
t2_of <- function(R) sapply(batches2, function(b)
  t2_manual(as.matrix(f2[f2$Batch == b, VARS]), R$mu, R$Si, I_LOTE))

res <- do.call(rbind, lapply(refs, function(R) {
  t2 <- t2_of(R)
  data.frame(Variant = R$label, Covariance = R$cov, Weights = R$wts,
             UCL = round(R$ucl, 2),
             Detected = sum(t2[ef] > R$ucl),
             FalseAlarms = sum(t2[!ef] > R$ucl),
             T2_median = round(median(t2[ef]), 1),
             stringsAsFactors = FALSE)
}))
rownames(res) <- NULL

# --- 3. Internal check: (a) and (d) must reproduce the pipeline ------
cat("\n===== INTERNAL CHECK: variants (a) and (d) against the pipeline =====\n")
t2_a <- t2_of(refs$a); t2_d <- t2_of(refs$d)
cat(sprintf("  (a) max |T2 - tc| = %.2e   %s\n", max(abs(t2_a - tc)),
            ifelse(max(abs(t2_a - tc)) < 1e-6, "OK", "*** FAILS ***")))
cat(sprintf("  (d) max |T2 - tr| = %.2e   %s\n", max(abs(t2_d - tr)),
            ifelse(max(abs(t2_d - tr)) < 1e-6, "OK", "*** FAILS ***")))

# --- 4. Anchors of the manuscript ------------------------------------
anchor <- data.frame(Detected = c(2, 20, 20, 20), FalseAlarms = c(0, 0, 0, 0),
                     T2_median = c(11.5, 56.5, 54.4, 88.7))
cat("\n===== ANCHORS (tep_ablation.csv) =====\n")
for (i in 1:4) {
  ok_cnt <- res$Detected[i] == anchor$Detected[i] &&
    res$FalseAlarms[i] == anchor$FalseAlarms[i]
  ok_med <- abs(res$T2_median[i] - anchor$T2_median[i]) <= 1.0
  cat(sprintf("  %-14s UCL %.2f | det %2d/20 | FA %d/10 | T2 med %6.1f (exp %5.1f)  counts:%s median:%s\n",
              res$Variant[i], res$UCL[i], res$Detected[i], res$FalseAlarms[i],
              res$T2_median[i], anchor$T2_median[i],
              ifelse(ok_cnt, "OK", "*** FAILS ***"),
              ifelse(ok_med, "OK", "REVIEW")))
}
print(res)

write.csv(res, file.path(DIR_OUT, "tep_ablation.csv"), row.names = FALSE)
cat(sprintf("\nSaved to %s/tep_ablation.csv\n", DIR_OUT))