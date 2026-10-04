# =====================================================================
# SCRIPT_M2b_sensitivity_sameseeds.R
# ---------------------------------------------------------------------
# Revision round 1, P10 (Section 3.5, Table 6).
#
# WHY
#   In SCRIPT_M2_power.R the rows of Table 6 came from three seed sets:
#   base and rho rows with the seeds of Table 4 (SEED_BASE*200/250/300),
#   h, J and I rows with *400/450/500, and 30% / 40% contamination with
#   *800/850/900 (inherited from the published scripts). The printed base
#   row is the one of Table 4, so the h, J and I rows were compared with a
#   base computed on other seeds (classical 0.2657 vs 0.2700). This script
#   recomputes those six rows with the seeds of Table 4, so that every row
#   of Table 6 shares one seed set and the base row is the reference of all
#   of them by construction.
#
# WHAT IT DOES
#   Same design, code and charts as SCRIPT_M2_power.R (one_rep copied
#   verbatim; only the seeds of the six rows change):
#     S2_h0.75, S2_J2, S2_J6, S2_I15  (6 contaminated batches, delta = 1)
#     S2_cont30 (9 batches), S2_cont40 (12 batches)
#   all with s = c(200, 250, 300). 2000 replicates each, 10 workers.
#   The base row (P_d1.0_c6) and the rho rows (S_rho0.3, S_rho0.9) already
#   use these seeds: they are read from 03_T2_crudos/M2/, not re-simulated.
#
# ANCHORS (the script stops if any fails)
#   B1  S2_base, replicates 1-200 only, rebuilt with this code, must equal
#       M2_P_d1.0_c6_0001.rds of SCRIPT_M2_power.R exactly (all 7 charts).
#       This proves the copied code and seeds reproduce M2.
#   B2  Means recomputed from the M2 files of P_d1.0_c6, S_rho0.3 and
#       S_rho0.9 must reproduce table_M2_all.csv (tolerance 5e-04).
#
# OUTPUT (05_revision_R1/02_resultados/) — M2 outputs are NOT overwritten
#   table_M2b_all.csv            the six new rows, 7 charts, means and SE
#   table_M2b_T6_sensitivity.csv complete Table 6 (9 rows) on one seed set
#   table_M2b_paired.csv         NEW - other, paired, |z| > 3 rule (as M2M3)
#   Means are written with 6 decimals to avoid rounding ties.
#   Raw replicates: 05_revision_R1/03_T2_crudos/M2b/ (resumable).
#   Run from C:/temp_paper. The serial check prints the time per replicate.
# =====================================================================
library(parallel)
library(robustT2AFM)
library(robustbase)
library(rrcov)
library(MASS)
stopifnot(packageVersion("robustT2AFM") >= "0.3.0")

DIR_OUT  <- "05_revision_R1/02_resultados"
DIR_M2   <- "05_revision_R1/03_T2_crudos/M2"
DIR_TASK <- "05_revision_R1/03_T2_crudos/M2b"
dir.create(DIR_TASK, recursive = TRUE, showWarnings = FALSE)

N_REP <- 2000; N_CAL <- 5000; N_CAL2 <- 50000; BLOCK <- 200; N_WORKERS <- 10
K2 <- 100; ALPHA <- 0.001; OR <- 0.20; OS <- 4; SEED_BASE <- 2026
CHARTS <- c("a_classical", "b_MCD_unif", "c_AFM_cls", "d_V7", "e_NEW", "f_MRCD", "g_RMCD_pooled")
make_equi <- function(p, rho) { S <- matrix(rho, p, p); diag(S) <- 1; S }

shift_vec <- function(dir, delta, J, rho) {
  Si <- solve(make_equi(J, rho)); m2 <- function(v) as.numeric(t(v) %*% Si %*% v)
  raw <- switch(dir, D1 = rep(1, J), D2 = c(1, -1, 0, 0), D3 = c(1, 1, -1, -1), D4 = c(1, 0, 0, 0),
                mD1 = -rep(1, J), mD4 = c(-1, 0, 0, 0))
  if (dir == "D1") return(rep(delta, J))
  if (dir == "mD1") return(-rep(delta, J))
  raw * sqrt(m2(rep(delta, J)) / m2(raw))
}

cfg <- function(id, block, dir, delta, ob, K = 30, I = 20, J = 4, h = 0.67, rho = 0.6, s = c(200, 250, 300), ncal = FALSE)
  data.frame(id = id, block = block, dir = dir, delta = delta, ob = ob, K = K, I = I, J = J, h = h, rho = rho,
             s1 = s[1], s2 = s[2], s3 = s[3], ncal = ncal, stringsAsFactors = FALSE)
CFG <- rbind(
  cfg("S2_base",   "S", "D1", 1, 6),               # anchor only (replicates 1-200)
  cfg("S2_h0.75",  "S", "D1", 1, 6, h = 0.75),
  cfg("S2_J2",     "S", "D1", 1, 6, J = 2),
  cfg("S2_J6",     "S", "D1", 1, 6, J = 6),
  cfg("S2_I15",    "S", "D1", 1, 6, I = 15),
  cfg("S2_cont30", "S", "D1", 1, 9),
  cfg("S2_cont40", "S", "D1", 1, 12))

# --- One replicate (verbatim from SCRIPT_M2_power.R) ------------------------
one_rep <- function(c, rep) {
  J <- c$J; I <- c$I; vv <- paste0("Var", 1:J); Sg <- make_equi(J, c$rho)
  v <- shift_vec(c$dir, c$delta, J, c$rho)
  set.seed(SEED_BASE * c$s1 + rep)
  sim <- simulate_batch_process(K1 = c$K, K2 = 1, I = I, J = J, rho = c$rho, Sigma = Sg,
                                outlier_batches_F1 = c$ob, outlier_rate = OR,
                                outlier_shift = OS, prop_ooc_F2 = 0)
  f1 <- subset(sim, Phase == "Phase 1"); b <- unique(f1$Batch)
  cal   <- suppressMessages(calibrate_afm_mcd(f1, vv, mcd_alpha = c$h))
  S_mcd <- cal$mcd_covariances; C <- do.call(rbind, cal$mcd_centers)[, vv, drop = FALSE]
  lam   <- function(S) max(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
  comb  <- function(Sl, w) Reduce(`+`, Map(function(S, wi) wi * S, Sl, w))
  wr    <- function(Sl) { l <- sapply(Sl, lam); (1 / l) / sum(1 / l) }
  S_cls <- lapply(b, function(bb) cov(f1[f1$Batch == bb, vv]))
  mu_cls <- colMeans(do.call(rbind, lapply(b, function(bb) colMeans(f1[f1$Batch == bb, vv]))))
  unif  <- rep(1 / length(b), length(b))
  X     <- as.matrix(f1[, vv])
  mr    <- suppressWarnings(rrcov::CovMrcd(X, alpha = c$h))
  rm_   <- robustbase::covMcd(X, alpha = c$h)
  R <- list(
    a_classical   = list(mu = mu_cls,          S = comb(S_cls, unif)),
    b_MCD_unif    = list(mu = colMeans(C),     S = comb(S_mcd, unif)),
    c_AFM_cls     = list(mu = mu_cls,          S = comb(S_cls, wr(S_cls))),
    d_V7          = list(mu = colMeans(C),     S = comb(S_mcd, wr(S_mcd))),
    e_NEW         = list(mu = cal$mu_r,        S = cal$Sw),
    f_MRCD        = list(mu = rrcov::getCenter(mr), S = rrcov::getCov(mr)),
    g_RMCD_pooled = list(mu = rm_$center,      S = rm_$cov))
  t2 <- function(M, r) { D <- sweep(M, 2, r$mu); I * rowSums((D %*% solve(r$S)) * D) }

  set.seed(SEED_BASE * c$s2 + rep)
  Mc <- do.call(rbind, lapply(seq_len(N_CAL), function(k) colMeans(MASS::mvrnorm(I, rep(0, J), Sg))))
  q  <- sapply(R, function(r) quantile(t2(Mc, r), 1 - ALPHA, names = FALSE))
  if (c$ncal) {
    set.seed(SEED_BASE * 260 + rep)
    Mc2 <- rbind(Mc, MASS::mvrnorm(N_CAL2 - N_CAL, rep(0, J), Sg / I))
    q2  <- sapply(R, function(r) quantile(t2(Mc2, r), 1 - ALPHA, names = FALSE))
  }
  set.seed(SEED_BASE * c$s3 + rep)
  M2 <- do.call(rbind, lapply(seq_len(K2), function(k) colMeans(MASS::mvrnorm(I, v, Sg))))
  tpr  <- sapply(names(R), function(nm) mean(t2(M2, R[[nm]]) > q[nm]))
  tpr2 <- if (c$ncal) sapply(names(R), function(nm) mean(t2(M2, R[[nm]]) > q2[nm])) else rep(NA, length(R))
  c(tpr, setNames(tpr2, paste0(names(R), "_ncal50k")))
}
run_task <- function(tk) {
  f <- file.path(DIR_TASK, sprintf("M2b_%s_%04d.rds", tk$id, tk$from))
  if (file.exists(f)) return(readRDS(f))
  c <- CFG[CFG$id == tk$id, ]
  out <- t(vapply(tk$from:tk$to, function(r) one_rep(c, r), numeric(2 * length(CHARTS))))
  saveRDS(out, f); out
}
read_m2 <- function(id) {
  f <- list.files(DIR_M2, paste0("^M2_", id, "_[0-9]{4}\\.rds$"), full.names = TRUE)
  M <- do.call(rbind, lapply(sort(f), readRDS)); stopifnot(nrow(M) == N_REP); M
}

# --- Serial check ------------------------------------------------------------
cat("=== Serial check: one replicate of S2_base ===\n")
t0 <- Sys.time(); print(round(one_rep(CFG[CFG$id == "S2_base", ], 1)[CHARTS], 4))
cat(sprintf("  OK (%.1f s per replicate). Launching...\n\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# --- Campaign: S2_base only block 1 (anchor); the six rows complete ----------
tasks <- rbind(data.frame(id = "S2_base", from = 1, to = BLOCK),
               do.call(rbind, lapply(CFG$id[-1], function(i)
                 data.frame(id = i, from = seq(1, N_REP, BLOCK), to = seq(BLOCK, N_REP, BLOCK)))))
tasks <- split(tasks, seq_len(nrow(tasks)))
t0 <- Sys.time()
clu <- makeCluster(N_WORKERS)
clusterEvalQ(clu, { library(robustT2AFM); library(robustbase); library(rrcov); library(MASS); NULL })
clusterExport(clu, c("CFG", "one_rep", "run_task", "shift_vec", "make_equi", "DIR_TASK", "N_CAL", "N_CAL2",
                     "K2", "ALPHA", "OR", "OS", "SEED_BASE", "CHARTS"))
out <- parLapplyLB(clu, tasks, run_task)
stopCluster(clu)
cat(sprintf("Campaign finished: %.1f h, %d workers, %d tasks\n\n",
            as.numeric(difftime(Sys.time(), t0, units = "hours")), N_WORKERS, length(tasks)))
ids <- vapply(tasks, function(t) t$id, character(1))

# --- Anchors -------------------------------------------------------------------
ok_all <- TRUE
anc <- out[[which(ids == "S2_base")]][, CHARTS]
ref <- readRDS(file.path(DIR_M2, "M2_P_d1.0_c6_0001.rds"))[, CHARTS]
b1 <- isTRUE(all.equal(unname(anc), unname(ref), tolerance = 1e-12))
cat(sprintf("B1  S2_base (reps 1-200) == M2 P_d1.0_c6 block 1: %s\n", if (b1) "OK" else "*** FAIL ***"))
ok_all <- ok_all && b1
m2all <- read.csv(file.path(DIR_OUT, "table_M2_all.csv"))
for (id in c("P_d1.0_c6", "S_rho0.3", "S_rho0.9")) {
  got <- colMeans(read_m2(id)[, CHARTS]); exp <- unlist(m2all[m2all$id == id, CHARTS])
  ok <- all(abs(got - exp) < 5e-4); ok_all <- ok_all && ok
  cat(sprintf("B2  %-10s recomputed means == table_M2_all.csv: %s\n", id, if (ok) "OK" else "*** FAIL ***"))
}
if (!ok_all) stop("ANCHOR FAILED: do NOT use these results.")

# --- Assemble --------------------------------------------------------------------
summ <- function(id, M, c) data.frame(
  c[, c("id", "K", "I", "J", "h", "rho", "ob")],
  t(round(colMeans(M[, CHARTS]), 6)),
  t(setNames(round(apply(M[, CHARTS], 2, function(x) sd(x) / sqrt(length(x))), 4), paste0("SE_", CHARTS))),
  check.names = FALSE, row.names = NULL)
NEWROWS <- do.call(rbind, lapply(CFG$id[-1], function(i) {
  M <- do.call(rbind, out[ids == i]); stopifnot(nrow(M) == N_REP)
  summ(i, M, CFG[CFG$id == i, ]) }))
MM <- c(setNames(lapply(CFG$id[-1], function(i) do.call(rbind, out[ids == i])), CFG$id[-1]),
        list(P_d1.0_c6 = read_m2("P_d1.0_c6"), S_rho0.3 = read_m2("S_rho0.3"), S_rho0.9 = read_m2("S_rho0.9")))
old <- rbind(cfg("P_d1.0_c6", "P", "D1", 1, 6), cfg("S_rho0.3", "S", "D1", 1, 6, rho = 0.3),
             cfg("S_rho0.9", "S", "D1", 1, 6, rho = 0.9))
OLDROWS <- do.call(rbind, lapply(old$id, function(i) summ(i, MM[[i]], old[old$id == i, ])))
ord <- c("P_d1.0_c6", "S2_cont30", "S2_cont40", "S2_h0.75", "S2_J2", "S2_J6", "S2_I15", "S_rho0.3", "S_rho0.9")
T6 <- rbind(OLDROWS, NEWROWS); T6 <- T6[match(ord, T6$id), ]
T6$advantage_NEW <- round(T6$e_NEW - T6$a_classical, 6)

PAIR <- do.call(rbind, lapply(ord, function(i) { M <- MM[[i]]
  do.call(rbind, lapply(c("a_classical", "d_V7", "f_MRCD", "g_RMCD_pooled"), function(o) {
    dd <- M[, "e_NEW"] - M[, o]; se <- sd(dd) / sqrt(length(dd)); z <- if (se > 0) mean(dd) / se else NA
    data.frame(id = i, versus = o, NEW = round(mean(M[, "e_NEW"]), 6), other = round(mean(M[, o]), 6),
               diff = round(mean(dd), 6), SE_paired = round(se, 4), z = round(z, 1),
               verdict = if (is.na(z) || abs(z) <= 3) "no difference" else if (z > 0) "NEW better" else "NEW worse")
  })) }))

write.csv(NEWROWS, file.path(DIR_OUT, "table_M2b_all.csv"),            row.names = FALSE)
write.csv(T6,      file.path(DIR_OUT, "table_M2b_T6_sensitivity.csv"), row.names = FALSE)
write.csv(PAIR,    file.path(DIR_OUT, "table_M2b_paired.csv"),         row.names = FALSE)
cat("\n=== Table 6, one seed set (seeds of Table 4) ===\n")
print(T6[, c("id", "a_classical", "SE_a_classical", "e_NEW", "SE_e_NEW", "advantage_NEW")], row.names = FALSE)
cat("\nALL ANCHORS OK: valid run. Files written in 05_revision_R1/02_resultados\n")
