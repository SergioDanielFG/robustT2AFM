# =====================================================================
# SCRIPT_M3_adverse_scenarios.R
# ---------------------------------------------------------------------
# Revision round 1, module M3: detection power under less favourable
# Phase 1 conditions, with the two robust competitors of M2 included in
# every scenario (Reviewer 1, comment 9; Reviewer 2, comment 5).
#
# WHAT THE REVIEWERS ASKED
#   R1-9: "smaller contamination magnitudes, heavy-tailed distributions,
#          skewed distributions and whole batches shifted mainly in
#          location but not in covariance".
#   R2-5: multivariate skew-normal and skew-t distributions.
#   Block V (batches with inflated dispersion) is OUR addition, not a
#   reviewer request: it is the situation the AFM weighting was designed
#   for, and after M2 it is needed to answer "why not pool all the
#   observations and use a single robust estimator?".
#
# CHARTS (identical code to M2; all brought to the same ARL0 through the
#   empirical 1 - alpha quantile over N_CAL in-control batches drawn from
#   the SAME in-control distribution as the scenario)
#   a_classical, b_MCD_unif, c_AFM_cls, d_V7, e_NEW, f_MRCD, g_RMCD_pooled
#
# BLOCKS (K = 30, I = 20, J = 4, h = 0.67, rho = 0.6, delta = 1,
#         Phase 2 shift in D1 = all variables equal, and in D2 = contrast
#         (1,-1,0,0) at the same Mahalanobis distance, as in M2)
#   A  Anchor: M2 configuration P_d1.0_c6 with the same seeds
#      (*200/250/300). Must reproduce the 7 charts of table_M2_all.csv.
#   H  In-control distribution t5, t3, skew-normal (SN) and skew-t with
#      5 df (ST5), Phase 1 clean or with 6 batches carrying 20% outlying
#      observations shifted 4 sd (the base contamination of the paper).
#      Every distribution is standardised to mean 0 and covariance Sigma,
#      so the shift is at the same Mahalanobis distance as in the normal
#      case. SN and ST5: sn::rmsn / sn::rmst with Omega = Sigma and
#      shape alpha = (8, 0, 0, 0), i.e. one strongly skewed variable
#      (marginal skewness about 0.93 for SN) and three mildly skewed.
#      Seeds *1000/1050/1100.
#   C  Smaller contamination magnitude: 6 batches, 20% of observations
#      shifted 1.5 or 2.5 sd (the paper uses 4). Seeds *1200/1250/1300.
#   W  Whole batches shifted in location, not in covariance: 6 of the 30
#      Phase 1 batches (20%) shifted 1.5 or 3 sd in every variable, no
#      outlying observations (simulate_batch_process, prop_contam_F1).
#      Seeds *1400/1450/1500.
#   V  Batches with inflated dispersion: 6 batches with covariance
#      kappa * Sigma, kappa = 2.25 or 4 (sd x1.5 or x2), same mean.
#      Seeds *1600/1650/1700.
#   29 configurations x 2000 replicates, the same number as the paper.
#
# EXPECTATIONS WRITTEN BEFORE RUNNING (reported whatever the outcome)
#   H  Heavy tails: the classical sample covariance is inefficient under
#      t3, so the robust charts should not lose to it; skewness is less
#      predictable, because MCD trims the long tail as if it were
#      contamination. With contamination, the D1/D2 pattern of M2 should
#      repeat.
#   C  The advantage of every robust chart over the classical one should
#      shrink as the shift decreases: at 1.5 sd the outliers are neither
#      harmful nor separable.
#   W  The within-batch covariance of NEW ignores a batch shifted as a
#      whole and its center is robust, whereas f and g fold the
#      between-batch shift into their covariance. NEW is expected to beat
#      f and g here. This is the scenario that tests the batch structure.
#   V  The AFM weights should down-weight the inflated batches, so d and
#      e are expected to beat b (uniform weights) and the classical chart.
#
# READING RULE FIXED IN ADVANCE
#   Each difference NEW - other is computed replicate by replicate (the
#   charts share the same data), with its paired standard error. A
#   difference is declared only when |z| > 3; otherwise "no difference".
#
# ANCHOR (tolerance 5e-04)
#   A_anchor = row P_d1.0_c6 of 02_resultados/table_M2_all.csv, 7 charts.
#
# GENERATOR CHECK (before launching)
#   1e6 in-control observations of t5, t3, SN and ST5: max |mean| < 0.01
#   and max |cov - Sigma| < 0.03 (t3 only reported: its fourth moment is
#   infinite and the sample covariance converges too slowly to test).
#
# PARALLEL AND RESUMABLE
#   29 configurations x 10 blocks of 200 replicates = 290 tasks over 10
#   workers. Each task saves its file in 03_T2_crudos/M3/; launching the
#   script again skips the finished tasks.
#
# OUTPUT (05_revision_R1/02_resultados/)
#   table_M3_all.csv, table_M3_paired.csv, table_M3_generator_check.csv
#   Run from C:/temp_paper. Needs package sn. Expected time: 3 to 5 hours.
# =====================================================================
library(parallel)
library(robustT2AFM)
library(robustbase)
library(rrcov)
library(MASS)
stopifnot(packageVersion("robustT2AFM") >= "0.3.0")
if (!requireNamespace("sn", quietly = TRUE))
  stop("Package 'sn' is needed for the skew-normal and skew-t scenarios: install.packages(\"sn\")")
cat("sn version:", as.character(packageVersion("sn")), "\n")

DIR_OUT  <- "05_revision_R1/02_resultados"
DIR_TASK <- "05_revision_R1/03_T2_crudos/M3"
dir.create(DIR_TASK, recursive = TRUE, showWarnings = FALSE)

N_REP <- 2000; N_CAL <- 5000; BLOCK <- 200; N_WORKERS <- 10
K2 <- 100; ALPHA <- 0.001; OR <- 0.20; SEED_BASE <- 2026
SKEW_ALPHA1 <- 8
CHARTS <- c("a_classical", "b_MCD_unif", "c_AFM_cls", "d_V7", "e_NEW", "f_MRCD", "g_RMCD_pooled")
make_equi <- function(p, rho) { S <- matrix(rho, p, p); diag(S) <- 1; S }

# --- Shift directions at the Mahalanobis distance of rep(delta, J) (as M2) --
shift_vec <- function(dir, delta, J, rho) {
  Si <- solve(make_equi(J, rho)); m2 <- function(v) as.numeric(t(v) %*% Si %*% v)
  raw <- switch(dir, D1 = rep(1, J), D2 = c(1, -1, 0, 0), D3 = c(1, 1, -1, -1), D4 = c(1, 0, 0, 0),
                mD1 = -rep(1, J), mD4 = c(-1, 0, 0, 0))
  if (dir == "D1") return(rep(delta, J))
  if (dir == "mD1") return(-rep(delta, J))
  raw * sqrt(m2(rep(delta, J)) / m2(raw))
}

# --- In-control distributions, standardised to mean 0 and covariance Sigma --
# SN / ST: Azzalini parametrisation (xi = 0, Omega = Sigma, alpha). Their mean
# and covariance are known in closed form (Azzalini & Capitanio 2003):
#   mu = omega * b * delta,  b = sqrt(2/pi) (SN) or
#        sqrt(nu/pi) Gamma((nu-1)/2) / Gamma(nu/2) (ST)
#   Cov = c * Omega - mu mu',  c = 1 (SN) or nu/(nu-2) (ST)
# and the linear map X = A (Y - mu), A = chol(Sigma) chol(Cov)^-1, gives
# covariance Sigma while keeping the skewness (it is affine invariant).
dist_std <- function(dist, J = 4, rho = 0.6) {
  S <- make_equi(J, rho)
  if (dist %in% c("normal", "t5", "t3")) return(list(dist = dist, S = S))
  nu <- if (dist == "SN") Inf else 5
  a  <- c(SKEW_ALPHA1, rep(0, J - 1))
  w  <- sqrt(diag(S)); Ob <- S / tcrossprod(w)
  d  <- as.vector(Ob %*% a) / sqrt(1 + as.numeric(t(a) %*% Ob %*% a))
  b  <- if (is.infinite(nu)) sqrt(2 / pi) else sqrt(nu / pi) * gamma((nu - 1) / 2) / gamma(nu / 2)
  mu <- w * b * d
  Cv <- (if (is.infinite(nu)) 1 else nu / (nu - 2)) * S - tcrossprod(mu)
  A  <- t(chol(S)) %*% solve(t(chol(Cv)))
  list(dist = dist, S = S, Omega = S, alpha = a, nu = nu, mu = mu, A = A)
}
DSTD <- lapply(setNames(nm = c("normal", "t5", "t3", "SN", "ST5")), dist_std)

r_obs <- function(n, D) {
  J <- nrow(D$S)
  switch(D$dist,
    normal = MASS::mvrnorm(n, rep(0, J), D$S),
    t5 = , t3 = {
      nu <- if (D$dist == "t5") 5 else 3
      MASS::mvrnorm(n, rep(0, J), D$S * (nu - 2) / nu) / sqrt(rchisq(n, nu) / nu)
    },
    SN  = sweep(sn::rmsn(n, xi = rep(0, J), Omega = D$Omega, alpha = D$alpha), 2, D$mu) %*% t(D$A),
    ST5 = sweep(sn::rmst(n, xi = rep(0, J), Omega = D$Omega, alpha = D$alpha, nu = D$nu), 2, D$mu) %*% t(D$A))
}
batch_means <- function(nb, I, D) rowsum(r_obs(nb * I, D), rep(seq_len(nb), each = I)) / I

# --- Phase 1 generated here (blocks H and V) --------------------------------
phase1_custom <- function(cf, D) {
  K <- cf$K; I <- cf$I; J <- cf$J
  bad <- if (cf$ob > 0) sort(sample.int(K, cf$ob)) else integer(0)
  Xl <- lapply(seq_len(K), function(k) {
    X <- r_obs(I, D)
    if (k %in% bad) {
      if (cf$kappa > 1) X <- X * sqrt(cf$kappa)
      if (cf$os > 0) {
        idx <- sample.int(I, max(1, round(I * OR)))
        X[idx, ] <- X[idx, ] + cf$os
      }
    }
    X
  })
  X <- do.call(rbind, Xl); colnames(X) <- paste0("Var", seq_len(J))
  data.frame(Batch = rep(sprintf("F1_B%02d", seq_len(K)), each = I), X, stringsAsFactors = FALSE)
}

# --- Configurations ------------------------------------------------------
cfg <- function(id, block, gen, dist, dir, ob = 0, os = 0, pc = 0, sc = 0, kappa = 1, s,
                delta = 1, K = 30, I = 20, J = 4, h = 0.67, rho = 0.6)
  data.frame(id = id, block = block, gen = gen, dist = dist, dir = dir, delta = delta,
             ob = ob, os = os, pc = pc, sc = sc, kappa = kappa, K = K, I = I, J = J, h = h, rho = rho,
             s1 = s[1], s2 = s[2], s3 = s[3], stringsAsFactors = FALSE)
DIRS <- c("D1", "D2")
CFG <- rbind(
  cfg("A_anchor", "A", "pkg", "normal", "D1", ob = 6, os = 4, s = c(200, 250, 300)),
  do.call(rbind, lapply(c("t5", "t3", "SN", "ST5"), function(ds) do.call(rbind, lapply(c(0, 6), function(o)
    do.call(rbind, lapply(DIRS, function(dr)
      cfg(sprintf("H_%s_c%d_%s", ds, o, dr), "H", "custom", ds, dr, ob = o, os = if (o > 0) 4 else 0,
          s = c(1000, 1050, 1100)))))))),
  do.call(rbind, lapply(c(1.5, 2.5), function(o) do.call(rbind, lapply(DIRS, function(dr)
    cfg(sprintf("C_os%.1f_%s", o, dr), "C", "pkg", "normal", dr, ob = 6, os = o, s = c(1200, 1250, 1300)))))),
  do.call(rbind, lapply(c(1.5, 3), function(sh) do.call(rbind, lapply(DIRS, function(dr)
    cfg(sprintf("W_sh%.1f_%s", sh, dr), "W", "pkg", "normal", dr, pc = 0.2, sc = sh, s = c(1400, 1450, 1500)))))),
  do.call(rbind, lapply(c(2.25, 4), function(kp) do.call(rbind, lapply(DIRS, function(dr)
    cfg(sprintf("V_k%.2f_%s", kp, dr), "V", "custom", "normal", dr, ob = 6, kappa = kp, s = c(1600, 1650, 1700)))))))
rownames(CFG) <- NULL
stopifnot(nrow(CFG) == 29, !anyDuplicated(CFG$id))

# --- One replicate (chart code identical to M2) ------------------------------
one_rep <- function(cf, rep) {
  J <- cf$J; I <- cf$I; vv <- paste0("Var", 1:J); Sg <- make_equi(J, cf$rho)
  v <- shift_vec(cf$dir, cf$delta, J, cf$rho)
  D <- DSTD[[cf$dist]]
  set.seed(SEED_BASE * cf$s1 + rep)
  if (cf$gen == "pkg") {
    sim <- simulate_batch_process(K1 = cf$K, K2 = 1, I = I, J = J, rho = cf$rho, Sigma = Sg,
                                  outlier_batches_F1 = cf$ob, outlier_rate = OR,
                                  outlier_shift = cf$os, prop_contam_F1 = cf$pc,
                                  shift_contam = cf$sc, prop_ooc_F2 = 0)
    f1 <- subset(sim, Phase == "Phase 1")
  } else {
    f1 <- phase1_custom(cf, D)
  }
  b <- unique(f1$Batch)
  cal   <- suppressMessages(calibrate_afm_mcd(f1, vv, mcd_alpha = cf$h))
  S_mcd <- cal$mcd_covariances; C <- do.call(rbind, cal$mcd_centers)[, vv, drop = FALSE]
  lam   <- function(S) max(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
  comb  <- function(Sl, w) Reduce(`+`, Map(function(S, wi) wi * S, Sl, w))
  wr    <- function(Sl) { l <- sapply(Sl, lam); (1 / l) / sum(1 / l) }
  S_cls <- lapply(b, function(bb) cov(f1[f1$Batch == bb, vv]))
  mu_cls <- colMeans(do.call(rbind, lapply(b, function(bb) colMeans(f1[f1$Batch == bb, vv]))))
  unif  <- rep(1 / length(b), length(b))
  X     <- as.matrix(f1[, vv])
  mr    <- suppressWarnings(rrcov::CovMrcd(X, alpha = cf$h))
  rm_   <- robustbase::covMcd(X, alpha = cf$h)
  R <- list(
    a_classical   = list(mu = mu_cls,          S = comb(S_cls, unif)),
    b_MCD_unif    = list(mu = colMeans(C),     S = comb(S_mcd, unif)),
    c_AFM_cls     = list(mu = mu_cls,          S = comb(S_cls, wr(S_cls))),
    d_V7          = list(mu = colMeans(C),     S = comb(S_mcd, wr(S_mcd))),
    e_NEW         = list(mu = cal$mu_r,        S = cal$Sw),
    f_MRCD        = list(mu = rrcov::getCenter(mr), S = rrcov::getCov(mr)),
    g_RMCD_pooled = list(mu = rm_$center,      S = rm_$cov))
  t2 <- function(M, r) { Dm <- sweep(M, 2, r$mu); I * rowSums((Dm %*% solve(r$S)) * Dm) }

  set.seed(SEED_BASE * cf$s2 + rep)
  Mc <- if (cf$dist == "normal") {
    do.call(rbind, lapply(seq_len(N_CAL), function(k) colMeans(MASS::mvrnorm(I, rep(0, J), Sg))))
  } else batch_means(N_CAL, I, D)
  q  <- sapply(R, function(r) quantile(t2(Mc, r), 1 - ALPHA, names = FALSE))
  set.seed(SEED_BASE * cf$s3 + rep)
  Mo <- if (cf$dist == "normal") {
    do.call(rbind, lapply(seq_len(K2), function(k) colMeans(MASS::mvrnorm(I, v, Sg))))
  } else sweep(batch_means(K2, I, D), 2, v, "+")
  sapply(names(R), function(nm) mean(t2(Mo, R[[nm]]) > q[nm]))
}
run_task <- function(tk) {
  f <- file.path(DIR_TASK, sprintf("M3_%s_%04d.rds", tk$id, tk$from))
  if (file.exists(f)) return(readRDS(f))
  cf <- CFG[CFG$id == tk$id, ]
  out <- t(vapply(tk$from:tk$to, function(r) one_rep(cf, r), numeric(length(CHARTS))))
  saveRDS(out, f); out
}

# --- Generator check -----------------------------------------------------
cat("=== Generator check: 1e6 in-control observations per distribution ===\n")
set.seed(SEED_BASE * 990)
GCHK <- do.call(rbind, lapply(c("t5", "t3", "SN", "ST5"), function(ds) {
  X  <- r_obs(1e6, DSTD[[ds]])
  sk <- apply(X, 2, function(x) mean((x - mean(x))^3) / sd(x)^3)
  data.frame(dist = ds, max_abs_mean = round(max(abs(colMeans(X))), 4),
             max_abs_cov_dev = round(max(abs(cov(X) - DSTD[[ds]]$S)), 4),
             skew_Var1 = round(sk[1], 3), skew_max_other = round(max(abs(sk[-1])), 3))
}))
print(GCHK, row.names = FALSE)
write.csv(GCHK, file.path(DIR_OUT, "table_M3_generator_check.csv"), row.names = FALSE)
if (!(all(GCHK$max_abs_mean < 0.01) && all(GCHK$max_abs_cov_dev[GCHK$dist != "t3"] < 0.03)))
  stop("Generator check FAILED: the in-control distributions do not have mean 0 and covariance Sigma. Nothing launched.")
cat("  Generator check OK.\n\n")

# --- Serial check: one replicate of each kind of generator --------------------
cat("=== Serial check ===\n")
for (i in c("A_anchor", "H_ST5_c6_D1", "V_k4.00_D2")) {
  t0 <- Sys.time(); r1 <- one_rep(CFG[CFG$id == i, ], 1)
  cat(sprintf("  %-12s %s  (%.1f s)\n", i, paste(sprintf("%.2f", r1), collapse = " "),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
cat("  OK. Launching...\n\n")

# --- Run ------------------------------------------------------------------------
tasks <- do.call(rbind, lapply(CFG$id, function(i) data.frame(id = i, from = seq(1, N_REP, BLOCK), to = seq(BLOCK, N_REP, BLOCK))))
tasks <- split(tasks, seq_len(nrow(tasks)))
t0 <- Sys.time()
clu <- makeCluster(N_WORKERS)
clusterEvalQ(clu, { library(robustT2AFM); library(robustbase); library(rrcov); library(MASS); NULL })
clusterExport(clu, c("CFG", "one_rep", "run_task", "shift_vec", "make_equi", "DSTD", "r_obs", "batch_means",
                     "phase1_custom", "DIR_TASK", "N_CAL", "K2", "ALPHA", "OR", "SEED_BASE", "CHARTS"))
out <- parLapplyLB(clu, tasks, run_task)
stopCluster(clu)
cat(sprintf("Run finished: %.1f h, %d workers, %d tasks\n\n",
            as.numeric(difftime(Sys.time(), t0, units = "hours")), N_WORKERS, length(tasks)))

# --- Assemble --------------------------------------------------------------------
ids <- vapply(tasks, function(t) t$id, character(1))
MM  <- lapply(setNames(nm = CFG$id), function(i) { M <- do.call(rbind, out[ids == i]); stopifnot(nrow(M) == N_REP); M })
ALL <- do.call(rbind, lapply(CFG$id, function(i) {
  M <- MM[[i]]; cf <- CFG[CFG$id == i, ]
  data.frame(cf[, c("id", "block", "dist", "dir", "ob", "os", "pc", "sc", "kappa")],
             t(round(colMeans(M[, CHARTS]), 4)),
             t(setNames(round(apply(M[, CHARTS], 2, function(x) sd(x) / sqrt(length(x))), 4), paste0("SE_", CHARTS))),
             check.names = FALSE, row.names = NULL)
}))
PAIR <- do.call(rbind, lapply(CFG$id, function(i) {
  M <- MM[[i]]; cf <- CFG[CFG$id == i, ]
  do.call(rbind, lapply(c("a_classical", "d_V7", "f_MRCD", "g_RMCD_pooled"), function(o) {
    dd <- M[, "e_NEW"] - M[, o]; se <- sd(dd) / sqrt(length(dd)); z <- if (se > 0) mean(dd) / se else NA
    data.frame(id = i, block = cf$block, dir = cf$dir, versus = o,
               NEW = round(mean(M[, "e_NEW"]), 4), other = round(mean(M[, o]), 4),
               diff = round(mean(dd), 4), SE_paired = round(se, 4), z = round(z, 1),
               verdict = if (is.na(z) || abs(z) <= 3) "no difference" else if (z > 0) "NEW better" else "NEW worse")
  }))
}))

# --- Anchor ----------------------------------------------------------------------
cat("=== Anchor: A_anchor must reproduce M2 P_d1.0_c6 ===\n")
m2  <- read.csv(file.path(DIR_OUT, "table_M2_all.csv"))
ref <- unlist(m2[m2$id == "P_d1.0_c6", CHARTS]); got <- unlist(ALL[ALL$id == "A_anchor", CHARTS])
ok_anchor <- all(abs(got - ref) < 5e-4)
for (k in seq_along(CHARTS)) cat(sprintf("  %-14s got %.4f expected %.4f %s\n", CHARTS[k], got[k], ref[k],
                                         if (abs(got[k] - ref[k]) < 5e-4) "OK" else "*** FAIL ***"))

write.csv(ALL,  file.path(DIR_OUT, "table_M3_all.csv"),    row.names = FALSE)
write.csv(PAIR, file.path(DIR_OUT, "table_M3_paired.csv"), row.names = FALSE)

cat("\n=== TPR by configuration (delta = 1) ===\n")
print(ALL[, c("id", CHARTS)], row.names = FALSE)
cat("\n=== NEW minus each chart, paired (|z| > 3 to declare a difference) ===\n")
print(PAIR[, c("id", "versus", "NEW", "other", "diff", "z", "verdict")], row.names = FALSE)
cat(if (ok_anchor) "\nANCHOR OK: valid run.\n" else "\nANCHOR FAILURE: do NOT use these results.\n")
