# =====================================================================
# SCRIPT_V_ucl_simulated_vs_M5.R
# ---------------------------------------------------------------------
# Check that the package function ucl_simulated() is the same procedure
# that was validated in module M5 (SCRIPT_M5_simulated_limit.R).
#
# It rebuilds three outer replicates of M5 (same seeds, same calls, in
# the same order) and computes their simulated limit with the package
# function. It must reproduce the UCL_sim saved by M5 in
# 03_T2_crudos/M5/ exactly (tolerance 1e-8): same random numbers, same
# calibration, same quantile.
#
#   K30_clean  r = 1 and r = 2
#   K30_contam r = 1
#
# Needs the package installed with ucl_simulated() (after commit 2).
# Run from C:/temp_paper. Expected time: about 6 minutes (3 x 500
# calibrations).
# =====================================================================
library(robustT2AFM)
stopifnot(exists("ucl_simulated", where = asNamespace("robustT2AFM")))

DIR_M5 <- "05_revision_R1/03_T2_crudos/M5"
I <- 20; J <- 4; ALPHA <- 0.001; RHO <- 0.6; OR <- 0.20; OS <- 4
B <- 500; M_IN <- 200
SEED_OUT <- 500000; SEED_IN <- 600000
VARS <- paste0("Var", 1:J)

check <- function(scen, K, ob, r) {
  saved <- readRDS(file.path(DIR_M5, sprintf("M5_%s_r%03d.rds", scen, r)))
  # Same Phase 1 and same calibration as outer_rep() in M5
  sim <- simulate_batch_process(K1 = K, K2 = 1, I = I, J = J, rho = RHO,
                                outlier_batches_F1 = ob, outlier_rate = OR,
                                outlier_shift = OS, prop_contam_F1 = 0,
                                prop_ooc_F2 = 0, seed = SEED_OUT + 1000 * K + r)
  ph1 <- subset(sim, Phase == "Phase 1")
  cal <- calibrate_afm_mcd(ph1, VARS)
  t0  <- Sys.time()
  us  <- ucl_simulated(cal, I = I, alpha = ALPHA, B = B, n_new = M_IN,
                       seed = SEED_IN + 1000 * K + r)$UCL
  ok  <- abs(us - saved$UCL_sim) < 1e-8
  cat(sprintf("  %-11s r = %d   package %.6f   M5 %.6f   %s   (%.1f min)\n",
              scen, r, us, saved$UCL_sim, if (ok) "OK" else "*** FAIL ***",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  ok
}

cat("=== ucl_simulated() against the saved M5 replicates ===\n")
res <- c(check("K30_clean",  30, 0, 1),
         check("K30_clean",  30, 0, 2),
         check("K30_contam", 30, 6, 1))
cat(if (all(res)) "\nALL CHECKS OK: the package function is the procedure validated in M5.\n"
    else "\nCHECK FAILURES: the package function does not reproduce M5. Do NOT commit.\n")
