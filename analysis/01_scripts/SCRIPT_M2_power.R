# =====================================================================
# SCRIPT_M2_power.R
# ---------------------------------------------------------------------
# Revision round 1, module M2: detection power with the final method,
# a contemporary robust competitor, shift directions and the size of
# the calibration set (Reviewer 1, comments 8, 10 and 11; Tables 4-6).
#
# CHARTS (all brought to the same ARL0 by the empirical 1 - alpha
# quantile over N_CAL in-control batches, as in the manuscript)
#   a_classical    classical T2 (batch covariances, uniform weights)
#   b_MCD_unif     per-batch MCD, uniform weights, mean of MCD centers
#   c_AFM_cls      classical batch covariances with AFM weights
#   d_V7           published AFM-MCD (robustT2AFM 0.2.0)
#   e_NEW          final AFM-MCD (robustT2AFM 0.3.0: MAD + MCD centers)
#   f_MRCD         MRCD of all Phase 1 observations pooled
#                  (Boudt et al. 2020; rrcov::CovMrcd), alpha = h
#   g_RMCD_pooled  reweighted MCD of all Phase 1 observations pooled
#                  (robustbase::covMcd), alpha = h
#   (a)-(d) are rebuilt from the same fits as the published scripts;
#   (f) and (g) use the same trimming h as the proposed method, so the
#   only difference is that they ignore the batch structure.
#
# BLOCKS (seeds and designs copied from the scripts of each table)
#   P  Tables 4 and 5: shift 0.5/1/1.5/2 in direction D1, Phase 1 with
#      0 or 6 contaminated batches (SCRIPT_FINAL_ABLATIONS_2x2.R,
#      seeds SEED_BASE*200/250/300)
#   D  Directions at delta = 1 and equal Mahalanobis distance: D2, D3,
#      D4 (SCRIPT_R1_10_shift_directions.R) and the new -D1 and -D4,
#      0 or 6 contaminated batches (seeds *200/250/300)
#   S  Table 6: base, h = 0.75, J = 2, J = 6, I = 15 (seeds *400/450/500),
#      30% and 40% contamination (seeds *800/850/900), rho = 0.3 and 0.9
#      (seeds *200/250/300), all with delta = 1 and 6 contaminated
#      batches unless the row says otherwise
#   N_CAL (comment 8): in the two P configurations with delta = 1, the
#      TPR is also computed with the quantile over 50 000 in-control
#      batches (the 5000 of the manuscript plus 45 000 more).
#
# ANCHORS (tolerance 5e-04; the published variants must be reproduced)
#   P  a, b, c, d = table_ablation_2x2_final.csv (8 rows)
#   D  a, b, c, d = table_R1_10_shift_directions.csv (D2, D3, D4)
#   S  classical and V7 = table_3_4_sensitivity.csv,
#      table_3_4b_contamination.csv, table_3_4c_rho.csv
#
# PARALLEL AND RESUMABLE
#   27 configurations x 10 blocks of 200 replicates = 270 tasks over 10
#   workers. Each task saves its file in 03_T2_crudos/M2/; launching
#   the script again skips the finished tasks.
#
# OUTPUT (05_revision_R1/02_resultados/)
#   table_M2_T4_power.csv, table_M2_T5_factorial.csv,
#   table_M2_T6_sensitivity.csv, table_M2_directions.csv,
#   table_M2_competitors.csv, table_M2_ncal.csv, table_M2_all.csv
#   Run from C:/temp_paper. Expected time: 2 to 4 hours.
# =====================================================================
library(parallel)
library(robustT2AFM)
library(robustbase)
library(rrcov)
library(MASS)
stopifnot(packageVersion("robustT2AFM") >= "0.3.0")

DIR_OUT  <- "05_revision_R1/02_resultados"
DIR_TASK <- "05_revision_R1/03_T2_crudos/M2"
dir.create(DIR_TASK, recursive = TRUE, showWarnings = FALSE)

N_REP <- 2000; N_CAL <- 5000; N_CAL2 <- 50000; BLOCK <- 200; N_WORKERS <- 10
K2 <- 100; ALPHA <- 0.001; OR <- 0.20; OS <- 4; SEED_BASE <- 2026
CHARTS <- c("a_classical", "b_MCD_unif", "c_AFM_cls", "d_V7", "e_NEW", "f_MRCD", "g_RMCD_pooled")
make_equi <- function(p, rho) { S <- matrix(rho, p, p); diag(S) <- 1; S }

# --- Shift directions at the Mahalanobis distance of rep(delta, J) -------
shift_vec <- function(dir, delta, J, rho) {
  Si <- solve(make_equi(J, rho)); m2 <- function(v) as.numeric(t(v) %*% Si %*% v)
  raw <- switch(dir, D1 = rep(1, J), D2 = c(1, -1, 0, 0), D3 = c(1, 1, -1, -1), D4 = c(1, 0, 0, 0),
                mD1 = -rep(1, J), mD4 = c(-1, 0, 0, 0))
  if (dir == "D1") return(rep(delta, J))
  if (dir == "mD1") return(-rep(delta, J))
  raw * sqrt(m2(rep(delta, J)) / m2(raw))
}

# --- Configurations ------------------------------------------------------
cfg <- function(id, block, dir, delta, ob, K = 30, I = 20, J = 4, h = 0.67, rho = 0.6, s = c(200, 250, 300), ncal = FALSE)
  data.frame(id = id, block = block, dir = dir, delta = delta, ob = ob, K = K, I = I, J = J, h = h, rho = rho,
             s1 = s[1], s2 = s[2], s3 = s[3], ncal = ncal, stringsAsFactors = FALSE)
CFG <- rbind(
  do.call(rbind, lapply(c(0.5, 1, 1.5, 2), function(d) rbind(
    cfg(sprintf("P_d%.1f_c0", d), "P", "D1", d, 0, ncal = d == 1),
    cfg(sprintf("P_d%.1f_c6", d), "P", "D1", d, 6, ncal = d == 1)))),
  do.call(rbind, lapply(c("D2", "D3", "D4", "mD1", "mD4"), function(dr) rbind(
    cfg(sprintf("D_%s_c0", dr), "D", dr, 1, 0), cfg(sprintf("D_%s_c6", dr), "D", dr, 1, 6)))),
  cfg("S_base",  "S", "D1", 1, 6, s = c(400, 450, 500)),
  cfg("S_h0.75", "S", "D1", 1, 6, h = 0.75, s = c(400, 450, 500)),
  cfg("S_J2",    "S", "D1", 1, 6, J = 2,    s = c(400, 450, 500)),
  cfg("S_J6",    "S", "D1", 1, 6, J = 6,    s = c(400, 450, 500)),
  cfg("S_I15",   "S", "D1", 1, 6, I = 15,   s = c(400, 450, 500)),
  cfg("S_cont30", "S", "D1", 1, 9,  s = c(800, 850, 900)),
  cfg("S_cont40", "S", "D1", 1, 12, s = c(800, 850, 900)),
  cfg("S_rho0.3", "S", "D1", 1, 6, rho = 0.3),
  cfg("S_rho0.9", "S", "D1", 1, 6, rho = 0.9))
rownames(CFG) <- NULL

# --- One replicate ---------------------------------------------------------
one_rep <- function(c, rep) {
  J <- c$J; I <- c$I; vv <- paste0("Var", 1:J); Sg <- make_equi(J, c$rho)
  v <- shift_vec(c$dir, c$delta, J, c$rho)
  set.seed(SEED_BASE * c$s1 + rep)
  sim <- simulate_batch_process(K1 = c$K, K2 = 1, I = I, J = J, rho = c$rho, Sigma = Sg,
                                outlier_batches_F1 = c$ob, outlier_rate = OR,
                                outlier_shift = OS, prop_ooc_F2 = 0)
  f1 <- subset(sim, Phase == "Phase 1"); b <- unique(f1$Batch)
  cal   <- suppressMessages(calibrate_afm_mcd(f1, vv, mcd_alpha = c$h))   # 0.3.0: same per-batch fits as 0.2.0
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
  f <- file.path(DIR_TASK, sprintf("M2_%s_%04d.rds", tk$id, tk$from))
  if (file.exists(f)) return(readRDS(f))
  c <- CFG[CFG$id == tk$id, ]
  out <- t(vapply(tk$from:tk$to, function(r) one_rep(c, r), numeric(2 * length(CHARTS))))
  saveRDS(out, f); out
}

# --- Serial check ---------------------------------------------------------
cat("=== Serial check: one replicate of P_d1.0_c6 ===\n")
t0 <- Sys.time(); print(round(one_rep(CFG[CFG$id == "P_d1.0_c6", ], 1), 4))
cat(sprintf("  OK (%.1f s per replicate). Launching...\n\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# --- Campaign ----------------------------------------------------------------
tasks <- do.call(rbind, lapply(CFG$id, function(i) data.frame(id = i, from = seq(1, N_REP, BLOCK), to = seq(BLOCK, N_REP, BLOCK))))
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

# --- Assemble ----------------------------------------------------------------
ids <- vapply(tasks, function(t) t$id, character(1))
ALL <- do.call(rbind, lapply(CFG$id, function(i) {
  M <- do.call(rbind, out[ids == i]); stopifnot(nrow(M) == N_REP)
  c <- CFG[CFG$id == i, ]
  data.frame(c[, c("id", "block", "dir", "delta", "ob", "K", "I", "J", "h", "rho")],
             t(round(colMeans(M[, CHARTS]), 4)),
             t(setNames(round(apply(M[, CHARTS], 2, function(x) sd(x) / sqrt(length(x))), 4), paste0("SE_", CHARTS))),
             t(setNames(round(colMeans(M[, paste0(CHARTS, "_ncal50k")]), 4), paste0(CHARTS, "_ncal50k"))),
             check.names = FALSE)
}))

# --- Anchors -------------------------------------------------------------------
tol <- 5e-4; ok_all <- TRUE
chk <- function(lbl, got, exp) { ok <- all(abs(got - exp) < tol); ok_all <<- ok_all && ok
  cat(sprintf("  %-26s got %-32s expected %-32s %s\n", lbl, paste(sprintf("%.4f", got), collapse = " "),
              paste(sprintf("%.4f", exp), collapse = " "), if (ok) "OK" else "*** FAIL ***")) }
cat("=== Anchors: published variants must be reproduced ===\n")
abl <- read.csv("02_resultados/table_ablation_2x2_final.csv")
for (k in seq_len(nrow(abl))) {
  r <- ALL[ALL$block == "P" & ALL$delta == abl$shift[k] & ALL$ob == abl$contam[k], ]
  chk(sprintf("P shift %.1f contam %d", abl$shift[k], abl$contam[k]),
      unlist(r[, c("a_classical", "b_MCD_unif", "c_AFM_cls", "d_V7")]),
      unlist(abl[k, c("TPR_a_classical", "TPR_b_MCD_only", "TPR_c_AFM_only", "TPR_d_AFM_MCD")]))
}
dirs <- read.csv("05_revision_R1/02_resultados/table_R1_10_shift_directions.csv")
map  <- c(D2_contrast_a = "D2", D3_contrast_b = "D3", D4_single_var = "D4")
for (k in which(dirs$direction %in% names(map))) {
  r <- ALL[ALL$block == "D" & ALL$dir == map[dirs$direction[k]] & ALL$ob == dirs$contam[k], ]
  chk(sprintf("D %s contam %d", map[dirs$direction[k]], dirs$contam[k]),
      unlist(r[, c("a_classical", "b_MCD_unif", "c_AFM_cls", "d_V7")]),
      unlist(dirs[k, c("TPR_a_classical", "TPR_b_MCD_only", "TPR_c_AFM_only", "TPR_d_AFM_MCD")]))
}
sens <- read.csv("02_resultados/table_3_4_sensitivity.csv")
sid  <- c(base = "S_base", `h=0.75` = "S_h0.75", `J=2` = "S_J2", `J=6` = "S_J6", `I=15` = "S_I15")
for (k in seq_len(nrow(sens))) { r <- ALL[ALL$id == sid[sens$config[k]], ]
  chk(paste("S", sens$config[k]), c(r$a_classical, r$d_V7), c(sens$TPR_Hot[k], sens$TPR_Rob[k])) }
cont <- read.csv("02_resultados/table_3_4b_contamination.csv")
for (k in seq_len(nrow(cont))) { r <- ALL[ALL$id == sprintf("S_cont%d", cont$contam_pct[k]), ]
  chk(paste("S", cont$config[k]), c(r$a_classical, r$d_V7), c(cont$TPR_Hot[k], cont$TPR_Rob[k])) }
rho <- read.csv("02_resultados/table_3_4c_rho.csv"); rho <- rho[rho$rho != 0.6, ]
for (k in seq_len(nrow(rho))) { r <- ALL[ALL$id == sprintf("S_rho%.1f", rho$rho[k]), ]
  chk(paste("S", rho$config[k]), c(r$a_classical, r$d_V7), c(rho$TPR_Hot[k], rho$TPR_Rob[k])) }

# --- Tables --------------------------------------------------------------------
P <- ALL[ALL$block == "P", ]
T4 <- P[, c("delta", "ob", "a_classical", "d_V7", "e_NEW", "f_MRCD", "g_RMCD_pooled",
            "SE_a_classical", "SE_d_V7", "SE_e_NEW", "SE_f_MRCD", "SE_g_RMCD_pooled")]
T5 <- P[, c("delta", "ob", "a_classical", "b_MCD_unif", "c_AFM_cls", "d_V7", "e_NEW",
            "SE_a_classical", "SE_b_MCD_unif", "SE_c_AFM_cls", "SE_d_V7", "SE_e_NEW")]
T6 <- ALL[ALL$block == "S", c("id", "K", "I", "J", "h", "rho", "ob", "a_classical", "d_V7", "e_NEW", "f_MRCD", "g_RMCD_pooled")]
T6$advantage_NEW <- round(T6$e_NEW - T6$a_classical, 4)
DIR <- ALL[(ALL$block == "P" & ALL$delta == 1) | ALL$block == "D",
           c("dir", "ob", "a_classical", "d_V7", "e_NEW", "f_MRCD", "g_RMCD_pooled")]
COMP <- ALL[, c("id", "ob", "e_NEW", "f_MRCD", "g_RMCD_pooled", "a_classical")]
COMP$NEW_minus_MRCD <- round(COMP$e_NEW - COMP$f_MRCD, 4)
COMP$NEW_minus_RMCD <- round(COMP$e_NEW - COMP$g_RMCD_pooled, 4)
NC <- do.call(rbind, lapply(which(CFG$ncal), function(i) { r <- ALL[ALL$id == CFG$id[i], ]
  data.frame(id = r$id, chart = CHARTS, TPR_ncal5000 = unlist(r[, CHARTS]),
             TPR_ncal50000 = unlist(r[, paste0(CHARTS, "_ncal50k")]),
             difference = round(unlist(r[, paste0(CHARTS, "_ncal50k")]) - unlist(r[, CHARTS]), 4), row.names = NULL) }))

write.csv(ALL,  file.path(DIR_OUT, "table_M2_all.csv"),          row.names = FALSE)
write.csv(T4,   file.path(DIR_OUT, "table_M2_T4_power.csv"),     row.names = FALSE)
write.csv(T5,   file.path(DIR_OUT, "table_M2_T5_factorial.csv"), row.names = FALSE)
write.csv(T6,   file.path(DIR_OUT, "table_M2_T6_sensitivity.csv"), row.names = FALSE)
write.csv(DIR,  file.path(DIR_OUT, "table_M2_directions.csv"),   row.names = FALSE)
write.csv(COMP, file.path(DIR_OUT, "table_M2_competitors.csv"),  row.names = FALSE)
write.csv(NC,   file.path(DIR_OUT, "table_M2_ncal.csv"),         row.names = FALSE)

cat("\n=== Table 4 (power) ===\n");        print(T4[, 1:7], row.names = FALSE)
cat("\n=== Directions (delta = 1) ===\n"); print(DIR, row.names = FALSE)
cat("\n=== Table 6 (sensitivity) ===\n");  print(T6, row.names = FALSE)
cat("\n=== N_CAL 5000 vs 50000 ===\n");    print(NC, row.names = FALSE)
cat(if (ok_all) "\nALL ANCHORS OK: valid run.\n" else "\nANCHOR FAILURES: do NOT use these results.\n")
