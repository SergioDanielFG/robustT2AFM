# =====================================================================
# SCRIPT_M5_simulated_limit.R
# ---------------------------------------------------------------------
# Revision round 1, module M5: simulation-based control limit
# (Reviewer 1, comment 2: "whether a simple numerical correction or
# simulation-based calibration could bring the operational ARL0 closer
# to its nominal target").
#
# THE PROCEDURE BEING EVALUATED (what a user would do with real data)
#   Given a Phase 1 calibration (mu_r, Sw) of the final method
#   (robustT2AFM 0.3.0 defaults):
#     1. generate a healthy synthetic Phase 1 of K batches x I
#        observations from N(mu_r, Sw);
#     2. calibrate it with the same method;
#     3. generate M_IN new healthy batch means from N(mu_r, Sw / I) and
#        compute their T2 against that synthetic calibration;
#     4. repeat B times and take the 1 - alpha quantile of all the T2.
#   Every synthetic Phase 1 goes through the same per-batch MCD, so the
#   small-sample inflation of the MCD covariance found in M1b is carried
#   into the limit automatically.
#
# PART 1: DOES IT GIVE ARL0 ~ 1000?  (nested simulation)
#   Scenarios: K = 30 clean, K = 30 contaminated (6 batches, 20% of
#   observations shifted 4 SD), K = 100 clean. For each scenario,
#   R_OUT = 200 "real" Phase 1 samples; for each one, the simulated
#   limit with B = 500 synthetic Phase 1 (100 000 T2), and the false
#   alarms of three limits on N_EVAL = 5000 in-control batches of the
#   true process: simulated limit, Eq. (8) with m* = 13, and the
#   classical chart with its own limit. 1 000 000 batches per scenario.
#   The outer replicates are 200, not 2000, because each one needs 500
#   calibrations; precision of the ARL0 is kept by evaluating 5000
#   batches per replicate (the same 1 000 000 batches as in M1).
#
# PART 2: TENNESSEE EASTMAN
#   The simulated limit for the TEP Phase 1 (6 IDV 7 batches), with
#   B = 1000, repeated with 100 seeds: spread of the limit across seeds
#   and detection of the 20 faulty / 10 healthy Phase 2 batches.
#
# DESIGN DECISIONS
#   - New batches are generated as batch means, x_bar ~ N(mu, Sigma / I),
#     which is exactly the distribution of the mean of I normal
#     observations; only the means enter T2.
#   - Each task (one outer replicate, or one TEP seed) sets its own
#     seeds, so results do not depend on the worker, and saves its own
#     file. If the run is interrupted, launching it again skips the
#     tasks already done (resumable overnight run).
#   - Parallel: 700 tasks over 10 workers with load balancing, largest
#     first (K = 100, then TEP, then K = 30).
#
# CHECKS
#   Eq. (8) limit 19.6929 (K = 30) and 18.8274 (K = 100). A serial test
#   with B = 20 runs before the campaign.
#
# OUTPUT
#   05_revision_R1/02_resultados/table_M5_simulated_limit.csv
#   05_revision_R1/02_resultados/table_M5_limit_spread.csv
#   05_revision_R1/02_resultados/table_M5_TEP.csv
#   05_revision_R1/03_T2_crudos/M5/  (one small file per task)
#   Run from C:/temp_paper. Expected time: 5 to 7 hours.
# =====================================================================
library(parallel)
library(robustT2AFM)
library(MASS)
stopifnot(packageVersion("robustT2AFM") >= "0.3.0")

DIR_OUT  <- "05_revision_R1/02_resultados"
DIR_TASK <- "05_revision_R1/03_T2_crudos/M5"
dir.create(DIR_TASK, recursive = TRUE, showWarnings = FALSE)

R_OUT <- 200; B <- 500; M_IN <- 200; N_EVAL <- 5000
B_TEP <- 1000; N_SEEDS_TEP <- 100; N_WORKERS <- 10
I <- 20; J <- 4; ALPHA <- 0.001; RHO <- 0.6; OB <- 6; OR <- 0.20; OS <- 4
SEED_OUT <- 500000; SEED_IN <- 600000; SEED_EVAL <- 700000; SEED_TEP <- 800000
VARS <- paste0("Var", 1:J)
Sigma_EQ <- matrix(RHO, J, J); diag(Sigma_EQ) <- 1
SCEN <- data.frame(scen = c("K100_clean", "K30_clean", "K30_contam"),
                   K = c(100, 30, 30), ob = c(0, 0, OB))

# --- The simulated limit ------------------------------------------------
t2_means <- function(Xbar, mu, Sw, I) { D <- sweep(Xbar, 2, mu); I * rowSums((D %*% solve(Sw)) * D) }
ucl_simulated <- function(mu, Sw, K, I, B, M_IN, alpha, vars) {
  J <- length(mu); t2 <- numeric(0)
  for (b in seq_len(B)) {
    X  <- MASS::mvrnorm(K * I, mu, Sw)
    f1 <- data.frame(Batch = rep(sprintf("S%03d", seq_len(K)), each = I), X)
    names(f1)[-1] <- vars
    cb <- calibrate_afm_mcd(f1, vars)
    xb <- MASS::mvrnorm(M_IN, mu, Sw / I)
    t2 <- c(t2, t2_means(xb, cb$mu_r, cb$Sw, I))
  }
  unname(quantile(t2, 1 - alpha, type = 8))
}

# --- One outer replicate of Part 1 ----------------------------------------
outer_rep <- function(scen, K, ob, r, B_use = B) {
  f <- file.path(DIR_TASK, sprintf("M5_%s_r%03d.rds", scen, r))
  if (B_use == B && file.exists(f)) return(readRDS(f))
  sim <- simulate_batch_process(K1 = K, K2 = 1, I = I, J = J, rho = RHO,
                                outlier_batches_F1 = ob, outlier_rate = OR,
                                outlier_shift = OS, prop_contam_F1 = 0,
                                prop_ooc_F2 = 0, seed = SEED_OUT + 1000 * K + r)
  ph1  <- subset(sim, Phase == "Phase 1")
  cal  <- calibrate_afm_mcd(ph1, VARS)
  u8   <- ucl_F_adjusted(cal, I = I, alpha = ALPHA)$UCL
  calH <- hotelling_classical_calibrate(ph1, VARS)
  uH   <- hotelling_classical_ucl(K = calH$n_batches, I = I, J = J, alpha = ALPHA, phase = "II")$UCL
  set.seed(SEED_IN + 1000 * K + r)
  us   <- ucl_simulated(cal$mu_r, cal$Sw, K, I, B_use, M_IN, ALPHA, VARS)
  set.seed(SEED_EVAL + 1000 * K + r)
  xe   <- MASS::mvrnorm(N_EVAL, rep(0, J), Sigma_EQ / I)
  t2n  <- t2_means(xe, cal$mu_r, cal$Sw, I)
  t2c  <- t2_means(xe, calH$mu_global, calH$Sp, I)
  out  <- data.frame(scen = scen, K = K, r = r, UCL_sim = us, UCL_eq8 = u8, UCL_cls = uH,
                     fa_sim = sum(t2n > us), fa_eq8 = sum(t2n > u8), fa_cls = sum(t2c > uH))
  if (B_use == B) saveRDS(out, f)
  out
}

# --- One TEP seed of Part 2 -----------------------------------------------
tep_seed <- function(s, mu, Sw) {
  f <- file.path(DIR_TASK, sprintf("M5_TEP_s%03d.rds", s))
  if (file.exists(f)) return(readRDS(f))
  set.seed(SEED_TEP + s)
  out <- data.frame(seed = s, UCL_sim = ucl_simulated(mu, Sw, 30, I, B_TEP, M_IN, ALPHA, names(mu)))
  saveRDS(out, f); out
}

# --- TEP calibration (master) ---------------------------------------------
cat("=== TEP: calibration of the real Phase 1 ===\n")
TVARS <- c("xmv_9", "xmv_8", "xmeas_19", "xmeas_17")
if (!exists("ff_test")) ff_test <- read.csv("04_datos/TEP_FaultFree_Testing.csv")
if (!exists("fy"))      fy      <- read.csv("04_datos/TEP_Faulty_Testing.csv")
lote <- function(dat, run, fault, desde, nom) {
  sel <- dat$simulationRun == run & dat$faultNumber == fault & dat$sample %in% seq(desde, by = 30, length.out = I)
  s <- dat[sel, TVARS]; if (nrow(s) != I) return(NULL); data.frame(Batch = nom, s, row.names = NULL) }
tf1 <- do.call(rbind, c(lapply(1:24, function(r) lote(ff_test, r, 0, 1,   sprintf("F1_sano_%02d", r))),
                        lapply(1:6,  function(r) lote(fy,      r, 7, 161, sprintf("F1_fallo_%02d", r)))))
tf2 <- do.call(rbind, c(lapply(25:34, function(r) lote(ff_test, r, 0, 1,   sprintf("F2_sano_%02d", r))),
                        lapply(7:26,  function(r) lote(fy,      r, 7, 161, sprintf("F2_fallo_%03d", r)))))
set.seed(2026)
tcal <- calibrate_afm_mcd(tf1, TVARS)
tmon <- monitor_afm_mcd(tf2, tcal, TVARS)
tucl8 <- ucl_F_adjusted(tcal, I = I, alpha = ALPHA)$UCL
cat(sprintf("  Eq. (8) UCL %.4f | %d/20 detected, %d/10 false alarms\n", tucl8,
            sum(tmon$T2[grepl("fallo", tmon$Batch)] > tucl8), sum(tmon$T2[grepl("sano", tmon$Batch)] > tucl8)))
rm(ff_test, fy); invisible(gc())                         # free memory before the cluster

# --- Checks ---------------------------------------------------------------
ucl_m <- function(K, m = 13) { d2 <- K * m - K - J + 1; (J * (K + 1) * (m - 1) / d2) * qf(1 - ALPHA, J, d2) }
stopifnot(abs(ucl_m(30) - 19.6929) < 5e-4, abs(ucl_m(100) - 18.8274) < 5e-4)
cat("\n=== Serial test: one outer replicate with B = 20 ===\n")
t0 <- Sys.time(); print(outer_rep("K30_contam", 30, OB, 1, B_use = 20))
cat(sprintf("  OK (%.1f s). Launching the campaign...\n\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# --- Campaign --------------------------------------------------------------
tasks <- c(
  lapply(seq_len(R_OUT), function(r) list(type = "out", scen = "K100_clean", K = 100, ob = 0, r = r)),
  lapply(seq_len(N_SEEDS_TEP), function(s) list(type = "tep", s = s)),
  unlist(lapply(2:3, function(i) lapply(seq_len(R_OUT), function(r)
    list(type = "out", scen = SCEN$scen[i], K = SCEN$K[i], ob = SCEN$ob[i], r = r))), recursive = FALSE))
run_task <- function(tk) if (tk$type == "out") outer_rep(tk$scen, tk$K, tk$ob, tk$r) else tep_seed(tk$s, TMU, TSW)
TMU <- tcal$mu_r; TSW <- tcal$Sw

t0 <- Sys.time()
clu <- makeCluster(N_WORKERS)
clusterEvalQ(clu, { library(robustT2AFM); library(MASS); NULL })
clusterExport(clu, c("t2_means", "ucl_simulated", "outer_rep", "tep_seed", "run_task",
                     "DIR_TASK", "B", "M_IN", "N_EVAL", "B_TEP", "I", "J", "ALPHA", "RHO",
                     "OR", "OS", "SEED_OUT", "SEED_IN", "SEED_EVAL", "SEED_TEP",
                     "VARS", "Sigma_EQ", "TMU", "TSW"))
res <- parLapplyLB(clu, tasks, run_task)
stopCluster(clu)
cat(sprintf("Campaign finished: %.1f h, %d workers, %d tasks\n\n",
            as.numeric(difftime(Sys.time(), t0, units = "hours")), N_WORKERS, length(tasks)))

# --- Part 1: summary ---------------------------------------------------------
P1 <- do.call(rbind, res[vapply(tasks, function(t) t$type == "out", logical(1))])
arl <- function(fa, n) c(ARL0 = round(n / fa), CI95_low = round(n / (qchisq(0.975, 2 * fa + 2) / 2)),
                         CI95_high = round(n / (qchisq(0.025, 2 * fa) / 2)))
tab <- do.call(rbind, lapply(split(P1, P1$scen), function(d) {
  n <- nrow(d) * N_EVAL
  rbind(data.frame(scen = d$scen[1], limit = "simulated (B = 500)", false_alarms = sum(d$fa_sim), batches = n, t(arl(sum(d$fa_sim), n))),
        data.frame(scen = d$scen[1], limit = "Eq. (8), m* = 13",     false_alarms = sum(d$fa_eq8), batches = n, t(arl(sum(d$fa_eq8), n))),
        data.frame(scen = d$scen[1], limit = "classical",            false_alarms = sum(d$fa_cls), batches = n, t(arl(sum(d$fa_cls), n))))
}))
spread <- do.call(rbind, lapply(split(P1, P1$scen), function(d) data.frame(
  scen = d$scen[1], replicates = nrow(d), UCL_eq8 = round(d$UCL_eq8[1], 4),
  UCL_sim_mean = round(mean(d$UCL_sim), 4), UCL_sim_sd = round(sd(d$UCL_sim), 4),
  UCL_sim_p05 = round(quantile(d$UCL_sim, 0.05), 4), UCL_sim_p95 = round(quantile(d$UCL_sim, 0.95), 4))))

# --- Part 2: TEP -------------------------------------------------------------
P2 <- do.call(rbind, res[vapply(tasks, function(t) t$type == "tep", logical(1))])
fl <- grepl("fallo", tmon$Batch)
P2$detected_of_20 <- vapply(P2$UCL_sim, function(u) sum(tmon$T2[fl] > u), numeric(1))
P2$false_alarms_of_10 <- vapply(P2$UCL_sim, function(u) sum(tmon$T2[!fl] > u), numeric(1))
tep <- data.frame(seeds = nrow(P2), UCL_eq8 = round(tucl8, 4),
                  UCL_sim_mean = round(mean(P2$UCL_sim), 4), UCL_sim_sd = round(sd(P2$UCL_sim), 4),
                  UCL_sim_min = round(min(P2$UCL_sim), 4), UCL_sim_max = round(max(P2$UCL_sim), 4),
                  detected_min = min(P2$detected_of_20), detected_max = max(P2$detected_of_20),
                  false_alarms_max = max(P2$false_alarms_of_10))

write.csv(tab,    file.path(DIR_OUT, "table_M5_simulated_limit.csv"), row.names = FALSE)
write.csv(spread, file.path(DIR_OUT, "table_M5_limit_spread.csv"),    row.names = FALSE)
write.csv(tep,    file.path(DIR_OUT, "table_M5_TEP.csv"),             row.names = FALSE)
cat("=== Part 1: ARL0 of each limit (1 000 000 batches per scenario) ===\n"); print(tab, row.names = FALSE)
cat("\n=== Spread of the simulated limit across Phase 1 samples ===\n");    print(spread, row.names = FALSE)
cat("\n=== Part 2: TEP, simulated limit over 100 seeds ===\n");             print(tep, row.names = FALSE)
cat(sprintf("\nReference from M1 (Eq. 8, NEW): K30 clean 1027, K30 contaminated 1258, K100 clean 1527.\n"))
