# =====================================================================
# SCRIPT_R3_02_computation_time_v030.R
# ---------------------------------------------------------------------
# Revision round 1, Reviewer 3, note 2: "a representative indication of
# computational time for Phase I calibration".
#
# WHAT IT COMPUTES
#   Wall-clock time of the Phase 1 calibration, calibrate_afm_mcd(), and
#   of its two cheap companions, ucl_F_adjusted() and the Phase 2
#   statistic for one batch, measured on one core:
#     (1) the simulation design of Section 3 (K = 30, I = 20, J = 4) and
#         four variations of it (K = 100, 300; I = 50; J = 8);
#     (2) the Tennessee Eastman Phase 1 of Section 4 (30 batches of 20
#         observations, 4 variables), built by PHASE1_PIPELINE_STEP30.R.
#   For reference it also times the classical calibration and the MCD
#   of the K batch centers (the robust center under study for comment 7
#   of Reviewer 1).
#
# VERSION v030 (29-sep-2026): same design as
# SCRIPT_R3_02_computation_time.R, which ran on 27-sep with robustT2AFM
# 0.2.0, whose calibration had neither the MAD scaling nor the MCD of
# the batch centers. This copy runs with 0.3.0, so calibrate_afm_mcd()
# times the final method (the MCD of the centers is now inside it; its
# separate row is kept as a component). It also times ucl_simulated()
# with its defaults (B = 500, n_new = 200) in the base design and on the
# Tennessee Eastman, N_TIME_SIM repetitions each, since the paper offers
# it as the alternative limit. Output: *_v030.csv / *_v030.txt.
#
# DESIGN DECISIONS
#   - Serial, one core: the figure a practitioner gets on a laptop, not
#     the throughput of the parallel Monte Carlo study.
#   - Each timing is repeated N_TIME times on fresh data; the median and
#     the 10th and 90th percentiles are reported, in milliseconds. Data
#     generation is outside the timed block.
#   - Run it with the machine otherwise idle (no other simulation in
#     the background), since timings depend on the load.
#
# OUTPUT
#   05_revision_R1/02_resultados/table_R3_02_computation_time_v030.csv
#   05_revision_R1/02_resultados/sessionInfo_R3_02_v030.txt
# =====================================================================
library(robustT2AFM)
library(MASS)
library(robustbase)
stopifnot(packageVersion("robustT2AFM") == "0.3.0")
OUT_DIR <- "05_revision_R1/02_resultados"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_TIME    <- 100
RHO       <- 0.6
H_MCD     <- 0.67
ALPHA     <- 0.001
SEED_BASE <- 2026
N_TIME_SIM <- 3

make_equi <- function(p, rho) { S <- matrix(rho, p, p); diag(S) <- 1; S }

# Median and 10th / 90th percentiles of the elapsed time, in ms.
time_it <- function(expr_fun, data_fun, n = N_TIME) {
  t <- numeric(n)
  for (i in seq_len(n)) {
    d <- data_fun(i)
    t[i] <- system.time(expr_fun(d))[["elapsed"]] * 1000
  }
  c(median_ms = median(t), p10_ms = unname(quantile(t, 0.10)),
    p90_ms = unname(quantile(t, 0.90)))
}

# Phase 1 generator for the simulation design (clean Phase 1).
gen_f1 <- function(K, I, J, seed) {
  set.seed(seed)
  sim <- simulate_batch_process(K1 = K, K2 = 0, I = I, J = J, rho = RHO,
                                Sigma = make_equi(J, RHO),
                                outlier_batches_F1 = 0, prop_ooc_F2 = 0)
  subset(sim, Phase == "Phase 1")
}

rows <- list()
add_row <- function(setting, K, I, J, step, tm)
  rows[[length(rows) + 1]] <<- data.frame(setting = setting, K = K, I = I, J = J,
                                          step = step, t(round(tm, 1)))

# =====================================================================
# (1) Simulation design and variations
# =====================================================================
DESIGNS <- data.frame(
  setting = c("Base (Section 3)", "K = 100", "K = 300", "I = 50", "J = 8"),
  K = c(30, 100, 300, 30, 30),
  I = c(20,  20,  20, 50, 20),
  J = c( 4,   4,   4,  4,  8))

cat("== Simulation designs ==\n")
for (d in seq_len(nrow(DESIGNS))) {
  K <- DESIGNS$K[d]; I <- DESIGNS$I[d]; J <- DESIGNS$J[d]
  VARS <- paste0("Var", seq_len(J))
  data_fun <- function(i) gen_f1(K, I, J, SEED_BASE * 400 + i)

  tm_cal <- time_it(function(f1) suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD)), data_fun)
  tm_cls <- time_it(function(f1) hotelling_classical_calibrate(f1, VARS), data_fun)

  f1  <- data_fun(1)
  cal <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
  C   <- do.call(rbind, cal$mcd_centers)
  tm_ctr <- time_it(function(x) robustbase::covMcd(x, alpha = 0.5), function(i) C)
  tm_ucl <- time_it(function(x) ucl_F_adjusted(x, I = I, alpha = ALPHA), function(i) cal)
  one_batch <- f1[f1$Batch == unique(f1$Batch)[1], ]
  tm_mon <- time_it(function(x) monitor_afm_mcd(x, cal, VARS), function(i) one_batch)
  if (d == 1) tm_sim <- time_it(function(x) ucl_simulated(x, I = I, alpha = ALPHA, seed = SEED_BASE),
                                function(i) cal, n = N_TIME_SIM)

  add_row(DESIGNS$setting[d], K, I, J, "calibrate_afm_mcd (Phase 1)", tm_cal)
  add_row(DESIGNS$setting[d], K, I, J, "classical calibration (reference)", tm_cls)
  add_row(DESIGNS$setting[d], K, I, J, "MCD of the K centers (comment 7)", tm_ctr)
  add_row(DESIGNS$setting[d], K, I, J, "ucl_F_adjusted", tm_ucl)
  add_row(DESIGNS$setting[d], K, I, J, "monitor_afm_mcd, one batch", tm_mon)
  if (d == 1) add_row(DESIGNS$setting[d], K, I, J, "ucl_simulated (B = 500, n_new = 200)", tm_sim)
  cat(sprintf("  %-18s calibration median %.1f ms [p10 %.1f, p90 %.1f]\n",
              DESIGNS$setting[d], tm_cal[1], tm_cal[2], tm_cal[3]))
}

# =====================================================================
# (2) Tennessee Eastman Phase 1 (Section 4)
# =====================================================================
cat("\n== Tennessee Eastman ==\n")
if (file.exists("01_scripts/PHASE1_PIPELINE_STEP30.R") &&
    file.exists("04_datos/TEP_FaultFree_Testing.csv")) {
  source("01_scripts/PHASE1_PIPELINE_STEP30.R")      # leaves f1, VARS, I_LOTE in memory
  VARS_TEP <- VARS
  tm_tep <- time_it(function(x) suppressMessages(calibrate_afm_mcd(x, VARS_TEP, mcd_alpha = H_MCD)),
                    function(i) f1)
  add_row("Tennessee Eastman (Section 4)", length(unique(f1$Batch)), I_LOTE,
          length(VARS_TEP), "calibrate_afm_mcd (Phase 1)", tm_tep)
  cal_tep <- suppressMessages(calibrate_afm_mcd(f1, VARS_TEP, mcd_alpha = H_MCD))
  tm_tep_sim <- time_it(function(x) ucl_simulated(x, I = I_LOTE, alpha = ALPHA, seed = SEED_BASE),
                        function(i) cal_tep, n = N_TIME_SIM)
  add_row("Tennessee Eastman (Section 4)", length(unique(f1$Batch)), I_LOTE,
          length(VARS_TEP), "ucl_simulated (B = 500, n_new = 200)", tm_tep_sim)
  cat(sprintf("  TEP calibration median %.1f ms [p10 %.1f, p90 %.1f]\n", tm_tep[1], tm_tep[2], tm_tep[3]))
} else {
  cat("  TEP files not found from the working directory: TEP timing skipped.\n")
}

# =====================================================================
# Output
# =====================================================================
res <- do.call(rbind, rows); rownames(res) <- NULL
cat("\n== RESULT (milliseconds, one core) ==\n")
print(res, row.names = FALSE)

write.csv(res, file.path(OUT_DIR, "table_R3_02_computation_time_v030.csv"), row.names = FALSE)
writeLines(c(capture.output(sessionInfo()), "",
             paste("CPU:", Sys.getenv("PROCESSOR_IDENTIFIER")),
             paste("Cores detected:", parallel::detectCores())),
           file.path(OUT_DIR, "sessionInfo_R3_02_v030.txt"))
cat("\nSaved to", file.path(OUT_DIR, "table_R3_02_computation_time_v030.csv"), "\n")
