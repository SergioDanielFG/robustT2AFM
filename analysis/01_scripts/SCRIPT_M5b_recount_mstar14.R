# =====================================================================
# SCRIPT_M5b_recount_mstar14.R
# ---------------------------------------------------------------------
# Revision round 1, module M5, second part (Reviewer 1, comments 2 and 3).
#
# WHY
#   SCRIPT_M5_simulated_limit.R (28-sep, 08:39 UTC) compared the simulated
#   limit with Eq. (8) computed with m* = 13 (round(I*h)), because the
#   m_star argument of ucl_F_adjusted() did not exist yet. The paper now
#   uses m* = 14 (MCD subset size). M5 saved, per outer replicate, the
#   limits and the false-alarm counts, not the T2 values, so the Eq. (8)
#   column cannot be recounted directly.
#
# WHAT IT DOES
#   For each of the 600 outer replicates of M5 (3 scenarios x 200) it
#   rebuilds, with the same seeds, the real Phase 1, its calibration and
#   the 5000 in-control evaluation batches, and counts false alarms of:
#     Eq. (8) with m* = 13 (anchor: must equal M5),
#     the simulated limit saved by M5 (anchor: must equal M5),
#     Eq. (8) with m* = 14 (new column).
#   The simulated limit itself is NOT recomputed (it is read from M5), so
#   the run takes minutes, not hours.
#
# ANCHORS
#   UCL of Eq. (8) with m* = 13 and both false-alarm counts must match the
#   M5 files exactly in every replicate. The script stops otherwise.
#
# OUTPUT
#   05_revision_R1/02_resultados/table_M5b_eq8_mstar14.csv
#   Run from C:/temp_paper.
# =====================================================================
library(robustT2AFM)
library(MASS)
stopifnot(packageVersion("robustT2AFM") == "0.3.0")

DIR_OUT  <- "05_revision_R1/02_resultados"
DIR_TASK <- "05_revision_R1/03_T2_crudos/M5"
R_OUT <- 200; N_EVAL <- 5000
I <- 20; J <- 4; ALPHA <- 0.001; RHO <- 0.6; OB <- 6; OR <- 0.20; OS <- 4
SEED_OUT <- 500000; SEED_EVAL <- 700000          # as in M5
VARS <- paste0("Var", 1:J)
Sigma_EQ <- matrix(RHO, J, J); diag(Sigma_EQ) <- 1
SCEN <- data.frame(scen = c("K100_clean", "K30_clean", "K30_contam"),
                   K = c(100, 30, 30), ob = c(0, 0, OB))

t2_means <- function(Xbar, mu, Sw, I) { D <- sweep(Xbar, 2, mu); I * rowSums((D %*% solve(Sw)) * D) }
arl_ci <- function(fa, n) c(ARL0 = if (fa > 0) round(n / fa) else NA,
  CI95_low = round(n / (qchisq(0.975, 2 * fa + 2) / 2)),
  CI95_high = if (fa > 0) round(n / (qchisq(0.025, 2 * fa) / 2)) else Inf)

rows <- list(); t0 <- Sys.time()
for (s in seq_len(nrow(SCEN))) {
  sc <- SCEN[s, ]; K <- sc$K
  cat(sprintf("--- %s ---\n", sc$scen))
  for (r in seq_len(R_OUT)) {
    saved <- readRDS(file.path(DIR_TASK, sprintf("M5_%s_r%03d.rds", sc$scen, r)))
    sim <- simulate_batch_process(K1 = K, K2 = 1, I = I, J = J, rho = RHO,
                                  outlier_batches_F1 = sc$ob, outlier_rate = OR,
                                  outlier_shift = OS, prop_contam_F1 = 0,
                                  prop_ooc_F2 = 0, seed = SEED_OUT + 1000 * K + r)
    ph1 <- subset(sim, Phase == "Phase 1")
    cal <- calibrate_afm_mcd(ph1, VARS)
    u13 <- ucl_F_adjusted(cal, I = I, alpha = ALPHA, m_star = "nominal")$UCL
    u14 <- ucl_F_adjusted(cal, I = I, alpha = ALPHA, m_star = "subset")$UCL
    set.seed(SEED_EVAL + 1000 * K + r)
    xe  <- MASS::mvrnorm(N_EVAL, rep(0, J), Sigma_EQ / I)
    t2n <- t2_means(xe, cal$mu_r, cal$Sw, I)
    rows[[length(rows) + 1]] <- data.frame(
      scen = sc$scen, r = r,
      UCL13 = u13, UCL13_M5 = saved$UCL_eq8, fa13 = sum(t2n > u13), fa13_M5 = saved$fa_eq8,
      fa_sim = sum(t2n > saved$UCL_sim), fa_sim_M5 = saved$fa_sim,
      UCL14 = u14, fa14 = sum(t2n > u14))
  }
}
res <- do.call(rbind, rows)

# --- Anchors -----------------------------------------------------------
ok <- isTRUE(all.equal(res$UCL13, res$UCL13_M5, tolerance = 1e-10)) &&
      all(res$fa13 == res$fa13_M5) && all(res$fa_sim == res$fa_sim_M5)
if (!ok) stop("ANCHOR FAILED: the rebuilt replicates do not match M5. Do not use.")
cat("\nAnchors OK: Eq. (8) m* = 13 and the simulated limit reproduce M5 exactly.\n")

# --- Summary -----------------------------------------------------------
out <- do.call(rbind, lapply(split(res, res$scen), function(d) {
  n <- nrow(d) * N_EVAL
  do.call(rbind, list(
    data.frame(scen = d$scen[1], limit = "simulated (B = 500)", UCL_mean = NA,
               false_alarms = sum(d$fa_sim), batches = n, t(arl_ci(sum(d$fa_sim), n))),
    data.frame(scen = d$scen[1], limit = "Eq. (8), m* = 13", UCL_mean = mean(d$UCL13),
               false_alarms = sum(d$fa13), batches = n, t(arl_ci(sum(d$fa13), n))),
    data.frame(scen = d$scen[1], limit = "Eq. (8), m* = 14", UCL_mean = mean(d$UCL14),
               false_alarms = sum(d$fa14), batches = n, t(arl_ci(sum(d$fa14), n)))))
}))
rownames(out) <- NULL
print(out, row.names = FALSE)
write.csv(out, file.path(DIR_OUT, "table_M5b_eq8_mstar14.csv"), row.names = FALSE)
cat(sprintf("\nSaved table_M5b_eq8_mstar14.csv | %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
