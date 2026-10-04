# =====================================================================
# SCRIPT_M1_limit.R
# ---------------------------------------------------------------------
# Revision round 1, module M1: the control limit under in-control
# conditions (Reviewer 1, comments 1, 2, 3 and 7; Tables 1, 2, 3, 13).
#
# WHAT IT COMPUTES
#   For every cell, 2000 replicates x 100 in-control Phase 2 batches
#   (200 000 batches per cell), the T2 of five charts, all against their
#   own analytic limit:
#     V7        published AFM-MCD: raw AFM weights, mean of MCD centers
#     SM_std    V7 weights, standardised spatial median of the centers
#               (as in Table 13 of the manuscript)
#     MCDc_V7w  V7 weights, MCD of the centers (isolates the center)
#     NEW       robustT2AFM 0.3.0 defaults: MAD scaling + MCD of centers
#     CLASSICAL Hotelling T2 with the Montgomery Phase II limit
#   The four AFM charts share the Eq. (8) limit, which depends only on
#   K, I, J and m*.
#
# CELLS (seed schemes copied from the scripts that produced each table)
#   T1  K = 30, Phase 1 contaminated, set.seed(2026*100 + r)
#       (ZTEST_TABLE1_FAITHFUL_v2.R)
#   T3  K = 30, seed = base + r, bases 2026/20000/40000/60000/80000,
#       Phase 1 clean and contaminated (20_arl0_clean_vs_contaminated_v2.R;
#       the base-2026 cells are also those of Table 13,
#       arl0_spatial_median_v2.R, and the clean one is the K = 30 row of
#       Table 2)
#   T2  K = 50, 100, 200, 500, Phase 1 clean, seed = 2026 + r
#       (15_arl0_convergence_K_v3.R)
#
# DESIGN DECISIONS
#   - One calibration per replicate with the installed package 0.3.0.
#     The per-batch MCD fits are identical to those of version 0.2.0
#     (same calls, same random stream), so V7 is rebuilt exactly from
#     them: raw first eigenvalues, inverse weights, mean of the centers.
#     The serial check below verifies this against
#     calibrate_afm_mcd(scaling = "none", center = "mean"), and the
#     anchors verify it against every published false-alarm count.
#   - m* (comment 3): in the K = 30 cells, before the calibration, the
#     per-batch MCD is refitted with the same random stream to record
#     the raw subset size (quan) and the number of observations kept by
#     the reweighting (sum of mcd.wt), for clean and contaminated
#     batches. The stream is then restored, so the calibration is
#     unaffected.
#   - All T2 are saved (03_T2_crudos), so the limit can be recomputed
#     afterwards for any m* or correction factor without simulating.
#   - Parallel: 11 cells + 4 cells split into blocks of 200 replicates,
#     spread over 10 workers with load balancing. Every replicate sets
#     its own seed, so the result does not depend on the worker.
#
# ANCHORS (V7 must reproduce the manuscript)
#   T1   V7 125 false alarms, SE 5.55e-05; classical 146, SE 6.35e-05
#   T2   V7 K = 50/100/200/500: 111 / 133 / 115 / 123
#   T3   V7 clean 125/160/145/145/154, contaminated 140/145/142/131/148
#   T13  SM_std clean 142, contaminated 122 (base 2026)
#
# OUTPUT
#   05_revision_R1/02_resultados/table_M1_T1_calibration.csv
#   05_revision_R1/02_resultados/table_M1_T2_K.csv
#   05_revision_R1/02_resultados/table_M1_T3_seeds.csv
#   05_revision_R1/02_resultados/table_M1_T3_pooled.csv
#   05_revision_R1/02_resultados/table_M1_T13_center.csv
#   05_revision_R1/02_resultados/table_M1_mstar.csv
#   05_revision_R1/03_T2_crudos/M1_<cell>.rds
#   Run from C:/temp_paper.
# =====================================================================
library(parallel)
library(robustT2AFM)
library(robustbase)
stopifnot(packageVersion("robustT2AFM") >= "0.3.0")

DIR_OUT <- "05_revision_R1/02_resultados"
DIR_RAW <- "05_revision_R1/03_T2_crudos"
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_RAW, recursive = TRUE, showWarnings = FALSE)

N_REP <- 2000; N_F2 <- 100; BLOCK <- 200; N_WORKERS <- 10
I <- 20; J <- 4; H_MCD <- 0.67; ALPHA <- 0.001; RHO <- 0.6
OB <- 6; OR <- 0.20; OS <- 4; SEED_BASE <- 2026
SEED_BASES <- c(2026, 20000, 40000, 60000, 80000)
VARS <- paste0("Var", 1:J)
METHODS <- c("V7", "SM_std", "MCDc_V7w", "NEW", "CLASSICAL")
make_equi <- function(p, rho) { S <- matrix(rho, p, p); diag(S) <- 1; S }
Sigma_EQ <- make_equi(J, RHO)
ucl_analitico <- function(K, I, J, h, alpha) {
  m <- round(I * h); d2 <- K * m - K - J + 1
  (J * (K + 1) * (m - 1) / d2) * qf(1 - alpha, J, d2)
}

# --- Cells -------------------------------------------------------------
cells <- rbind(
  data.frame(cell = "T1_K30_contam", table = "T1", K = 30, ob = OB, scheme = "set100", base = SEED_BASE),
  do.call(rbind, lapply(SEED_BASES, function(b) rbind(
    data.frame(cell = sprintf("T3_b%d_clean",  b), table = "T3", K = 30, ob = 0,  scheme = "arg", base = b),
    data.frame(cell = sprintf("T3_b%d_contam", b), table = "T3", K = 30, ob = OB, scheme = "arg", base = b)))),
  do.call(rbind, lapply(c(50, 100, 200, 500), function(k)
    data.frame(cell = sprintf("T2_K%d_clean", k), table = "T2", K = k, ob = 0, scheme = "arg", base = SEED_BASE))))
rownames(cells) <- NULL

# --- One replicate -----------------------------------------------------
run_rep <- function(cl, r) {
  if (cl$scheme == "set100") {
    set.seed(SEED_BASE * 100 + r)
    sim <- simulate_batch_process(K1 = cl$K, K2 = N_F2, I = I, J = J, rho = RHO, Sigma = Sigma_EQ,
                                  outlier_batches_F1 = cl$ob, outlier_rate = OR,
                                  outlier_shift = OS, prop_ooc_F2 = 0)
  } else {
    sim <- simulate_batch_process(K1 = cl$K, K2 = N_F2, I = I, J = J, rho = RHO,
                                  outlier_batches_F1 = cl$ob, outlier_rate = OR,
                                  outlier_shift = OS, prop_contam_F1 = 0,
                                  prop_ooc_F2 = 0, seed = cl$base + r)
  }
  ph1 <- subset(sim, Phase == "Phase 1"); ph2 <- subset(sim, Phase == "Phase 2")

  # m*: refit the per-batch MCD with the same stream, then restore it
  mstar <- c(quan = NA, kept_clean = NA, kept_cont = NA)
  if (cl$K == 30) {
    st <- get(".Random.seed", envir = globalenv())
    b  <- unique(ph1$Batch)
    ct <- vapply(b, function(bb) any(as.character(ph1$ContaminationType[ph1$Batch == bb]) != "Clean"), logical(1))
    fit <- lapply(b, function(bb) covMcd(ph1[ph1$Batch == bb, VARS, drop = FALSE], alpha = H_MCD))
    kept <- vapply(fit, function(f) sum(f$mcd.wt), numeric(1))
    mstar <- c(quan = mean(vapply(fit, function(f) f$quan, numeric(1))),
               kept_clean = mean(kept[!ct]),
               kept_cont  = if (any(ct)) mean(kept[ct]) else NA)
    assign(".Random.seed", st, envir = globalenv())
  }

  cal <- suppressMessages(calibrate_afm_mcd(ph1, VARS, mcd_alpha = H_MCD))
  S   <- cal$mcd_covariances
  C   <- do.call(rbind, cal$mcd_centers)[, VARS, drop = FALSE]
  lam <- sapply(S, function(x) max(eigen(x, symmetric = TRUE, only.values = TRUE)$values))
  w7  <- (1 / lam) / sum(1 / lam)
  Sw7 <- Reduce(`+`, Map(function(x, wi) wi * x, S, w7))
  mu7 <- colMeans(C)
  sd1 <- apply(ph1[, VARS], 2, sd)
  sm  <- local({ X <- sweep(C, 2, sd1, "/"); m <- colMeans(X)
    for (it in 1:2000) { d <- sqrt(rowSums((X - matrix(m, nrow(X), ncol(X), byrow = TRUE))^2))
      d[d < 1e-12] <- 1e-12; w <- 1 / d; mn <- colSums(X * w) / sum(w)
      if (max(abs(mn - m)) < 1e-10) { m <- mn; break }; m <- mn }
    m * sd1 })

  Xb <- rowsum(as.matrix(ph2[, VARS]), ph2$Batch, reorder = FALSE) / I
  t2 <- function(mu, Sw) { D <- sweep(Xb, 2, mu); I * rowSums((D %*% solve(Sw)) * D) }
  calH <- hotelling_classical_calibrate(ph1, VARS)
  mH   <- hotelling_classical_monitor(ph2, calH, VARS)
  T2 <- cbind(V7 = t2(mu7, Sw7), SM_std = t2(sm, Sw7), MCDc_V7w = t2(cal$mu_r, Sw7),
              NEW = t2(cal$mu_r, cal$Sw), CLASSICAL = mH$T2[match(rownames(Xb), mH$Batch)])
  list(T2 = T2, mstar = mstar,
       ucl_cls = hotelling_classical_ucl(K = calH$n_batches, I = I, J = J, alpha = ALPHA, phase = "II")$UCL)
}
run_block <- function(tk) {
  cl <- cells[cells$cell == tk$cell, ]
  lapply(tk$from:tk$to, function(r) run_rep(cl, r))
}

# --- Serial check: replicate 1 of T1 -----------------------------------
cat("=== Serial check (replicate 1 of T1) ===\n")
cl1 <- cells[1, ]
set.seed(SEED_BASE * 100 + 1)
sim <- simulate_batch_process(K1 = 30, K2 = N_F2, I = I, J = J, rho = RHO, Sigma = Sigma_EQ,
                              outlier_batches_F1 = OB, outlier_rate = OR, outlier_shift = OS, prop_ooc_F2 = 0)
ph1 <- subset(sim, Phase == "Phase 1"); ph2 <- subset(sim, Phase == "Phase 2")
st  <- .Random.seed
cal_new <- calibrate_afm_mcd(ph1, VARS, mcd_alpha = H_MCD)
.Random.seed <- st
cal_old <- calibrate_afm_mcd(ph1, VARS, mcd_alpha = H_MCD, scaling = "none", center = "mean")
S <- cal_new$mcd_covariances
lam <- sapply(S, function(x) max(eigen(x, symmetric = TRUE, only.values = TRUE)$values))
w7 <- (1 / lam) / sum(1 / lam); Sw7 <- Reduce(`+`, Map(function(x, wi) wi * x, S, w7))
mu7 <- colMeans(do.call(rbind, cal_new$mcd_centers))
d_rebuild <- max(abs(Sw7 - cal_old$Sw), abs(mu7 - cal_old$mu_r))
r1 <- run_rep(cl1, 1)
m_old <- monitor_afm_mcd(ph2, cal_old, VARS); m_new <- monitor_afm_mcd(ph2, cal_new, VARS)
d_t2 <- max(abs(r1$T2[, "V7"] - m_old$T2[match(rownames(r1$T2), m_old$Batch)]),
            abs(r1$T2[, "NEW"] - m_new$T2[match(rownames(r1$T2), m_new$Batch)]))
cat(sprintf("  V7 rebuilt vs calibrate(scaling='none', center='mean'): max diff %.1e\n", d_rebuild))
cat(sprintf("  manual T2 vs monitor_afm_mcd (V7 and NEW): max diff %.1e\n", d_t2))
stopifnot(d_rebuild < 1e-10, d_t2 < 1e-8)
cat("  OK. Launching the campaign...\n\n")

# --- Parallel campaign -------------------------------------------------
tasks <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
  data.frame(cell = cells$cell[i], K = cells$K[i],
             from = seq(1, N_REP, by = BLOCK), to = seq(BLOCK, N_REP, by = BLOCK))))
tasks <- tasks[order(-tasks$K), ]                       # largest K first
tasks <- split(tasks, seq_len(nrow(tasks)))

t0 <- Sys.time()
clu <- makeCluster(N_WORKERS)
clusterEvalQ(clu, { library(robustT2AFM); library(robustbase); NULL })
clusterExport(clu, c("cells", "run_rep", "N_F2", "I", "J", "H_MCD", "ALPHA", "RHO",
                     "OB", "OR", "OS", "SEED_BASE", "VARS", "Sigma_EQ"))
out <- parLapplyLB(clu, tasks, run_block)
stopCluster(clu)
cat(sprintf("Campaign finished: %.1f min, %d workers, %d tasks\n\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), N_WORKERS, length(tasks)))

# --- Assemble per cell ---------------------------------------------------
cell_of <- vapply(tasks, function(t) t$cell, character(1))
from_of <- vapply(tasks, function(t) t$from, numeric(1))
res <- list()
for (cc in cells$cell) {
  idx  <- which(cell_of == cc); idx <- idx[order(from_of[idx])]
  reps <- do.call(c, out[idx]); stopifnot(length(reps) == N_REP)
  K    <- cells$K[cells$cell == cc]
  arr  <- array(NA_real_, c(N_REP, N_F2, length(METHODS)), dimnames = list(NULL, NULL, METHODS))
  for (r in seq_len(N_REP)) arr[r, , ] <- reps[[r]]$T2
  ucl  <- c(rep(ucl_analitico(K, I, J, H_MCD, ALPHA), 4), reps[[1]]$ucl_cls); names(ucl) <- METHODS
  saveRDS(list(T2 = arr, ucl = ucl, cell = cells[cells$cell == cc, ],
               mstar = do.call(rbind, lapply(reps, `[[`, "mstar"))),
          file.path(DIR_RAW, paste0("M1_", cc, ".rds")))
  res[[cc]] <- list(arr = arr, ucl = ucl, mstar = do.call(rbind, lapply(reps, `[[`, "mstar")))
}

resumen <- function(arr, ucl, m) {
  fa_rep <- rowSums(arr[, , m] > ucl[m]); fa <- sum(fa_rep); n <- length(arr[, , m])
  lo <- n / (qchisq(0.975, 2 * fa + 2) / 2); hi <- if (fa > 0) n / (qchisq(0.025, 2 * fa) / 2) else Inf
  data.frame(method = m, UCL = round(ucl[m], 4), false_alarms = fa, FAR = fa / n,
             SE_FAR = sd(fa_rep / N_F2) / sqrt(N_REP), ARL0 = if (fa > 0) round(n / fa) else NA,
             CI95_low = round(lo), CI95_high = round(hi))
}
tab_cell <- function(cc) do.call(rbind, lapply(METHODS, function(m) resumen(res[[cc]]$arr, res[[cc]]$ucl, m)))

T1 <- cbind(cell = "T1_K30_contam", tab_cell("T1_K30_contam"))
T3 <- do.call(rbind, lapply(cells$cell[cells$table == "T3"], function(cc)
  cbind(seed_base = cells$base[cells$cell == cc], phase1 = ifelse(cells$ob[cells$cell == cc] > 0, "contaminated", "clean"), tab_cell(cc))))
T2 <- do.call(rbind, lapply(c("T3_b2026_clean", cells$cell[cells$table == "T2"]), function(cc)
  cbind(K = cells$K[cells$cell == cc], tab_cell(cc))))
pooled <- do.call(rbind, lapply(c("clean", "contaminated"), function(ph) {
  cc <- cells$cell[cells$table == "T3" & (cells$ob > 0) == (ph == "contaminated")]
  do.call(rbind, lapply(METHODS, function(m) {
    fa <- sum(vapply(cc, function(x) sum(res[[x]]$arr[, , m] > res[[x]]$ucl[m]), numeric(1)))
    n  <- length(cc) * N_REP * N_F2
    data.frame(phase1 = ph, method = m, false_alarms = fa, batches = n, ARL0 = round(n / fa),
               CI95_low = round(n / (qchisq(0.975, 2 * fa + 2) / 2)), CI95_high = round(n / (qchisq(0.025, 2 * fa) / 2)))
  }))
}))
T13 <- T3[T3$seed_base == 2026, ]
mst <- do.call(rbind, lapply(cells$cell[cells$K == 30], function(cc) {
  m <- res[[cc]]$mstar
  data.frame(cell = cc, quan_raw_subset = round(mean(m[, "quan"]), 3),
             kept_after_reweighting_clean = round(mean(m[, "kept_clean"]), 3),
             kept_after_reweighting_contaminated = round(mean(m[, "kept_cont"], na.rm = TRUE), 3))
}))

# --- Anchors -------------------------------------------------------------
fa <- function(tab, m, ...) { f <- list(...); s <- tab[tab$method == m, ]
  for (nm in names(f)) s <- s[s[[nm]] == f[[nm]], ]; s$false_alarms }
chk <- function(lbl, got, exp) { ok <- all(got == exp)
  cat(sprintf("  %-38s got %-28s expected %-28s %s\n", lbl, paste(got, collapse = "/"),
              paste(exp, collapse = "/"), if (ok) "OK" else "*** FAIL ***")); ok }
cat("=== Anchors: V7 must reproduce the manuscript ===\n")
a <- c(
  chk("T1 V7 false alarms", fa(T1, "V7"), 125),
  chk("T1 classical false alarms", fa(T1, "CLASSICAL"), 146),
  abs(T1$SE_FAR[T1$method == "V7"] - 5.55e-05) < 2e-06,
  abs(T1$SE_FAR[T1$method == "CLASSICAL"] - 6.35e-05) < 2e-06,
  chk("T2 V7 K=30/50/100/200/500", fa(T2, "V7"), c(125, 111, 133, 115, 123)),
  chk("T3 V7 clean (5 bases)", fa(T3, "V7", phase1 = "clean"), c(125, 160, 145, 145, 154)),
  chk("T3 V7 contaminated (5 bases)", fa(T3, "V7", phase1 = "contaminated"), c(140, 145, 142, 131, 148)),
  chk("T13 spatial median clean/contam", fa(T13, "SM_std"), c(142, 122)))
cat(sprintf("  T1 SE V7 %.3e (5.55e-05) | classical %.3e (6.35e-05) -> %s\n",
            T1$SE_FAR[T1$method == "V7"], T1$SE_FAR[T1$method == "CLASSICAL"],
            if (a[3] && a[4]) "OK" else "*** FAIL ***"))

# --- Output --------------------------------------------------------------
write.csv(T1,     file.path(DIR_OUT, "table_M1_T1_calibration.csv"), row.names = FALSE)
write.csv(T2,     file.path(DIR_OUT, "table_M1_T2_K.csv"),           row.names = FALSE)
write.csv(T3,     file.path(DIR_OUT, "table_M1_T3_seeds.csv"),       row.names = FALSE)
write.csv(pooled, file.path(DIR_OUT, "table_M1_T3_pooled.csv"),      row.names = FALSE)
write.csv(T13,    file.path(DIR_OUT, "table_M1_T13_center.csv"),     row.names = FALSE)
write.csv(mst,    file.path(DIR_OUT, "table_M1_mstar.csv"),          row.names = FALSE)

cat("\n=== Pooled ARL0 over the 5 seed bases (1 000 000 batches per condition) ===\n")
print(pooled, row.names = FALSE)
cat("\n=== m*: observations actually used by the per-batch MCD ===\n")
print(mst, row.names = FALSE)
cat(if (all(a)) "\nALL ANCHORS OK: valid run.\n" else "\nANCHOR FAILURES: do NOT use these results.\n")
