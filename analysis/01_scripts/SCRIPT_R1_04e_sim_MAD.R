# =====================================================================
# SCRIPT_R1_04e_sim_MAD.R
# ---------------------------------------------------------------------
# Revision round 1, Reviewer 1, comment 4. Simulation with the MAD.
#
# WHY THIS SCRIPT
#   SCRIPT_R1_04d showed on the TEP data that scaling each variable by
#   its MAD (computed on all Phase 1 observations) works as well as the
#   median of the within-batch MCD standard deviations. The MAD is the
#   simplest robust scale to explain and cite, so this script checks
#   that it also holds in the simulation of Table 5.
#
# CHARTS (same data, seeds and ARL0 equalisation as SCRIPT_R1_04b)
#   a  classical T2
#   d  AFM-MCD published (no scaling)
#   f  AFM-MCD, divisor = median of the within-batch MCD std. deviations
#   g  AFM-MCD, divisor = MAD of each variable over all Phase 1 data
#      (robust standardisation of the variables before the weights)
#   z  AFM-MCD, divisor = classical standard deviation of each variable
#      over all Phase 1 data (traditional z-score standardisation, the
#      literal "standardized variables" of the reviewer; non-robust)
#
# PARALLELISATION
#   The 2000 replicates of each of the 6 configurations are split into
#   blocks of 200, so 60 tasks are spread over the workers. Every
#   replicate sets its own seeds (SEED_BASE*200/250/300 + rep), so the
#   result does not depend on which worker runs it.
#
# ANCHORS (tolerance 5e-04)
#   E1  a, d, f reproduce table_R1_04b_global_standardization.csv.
#   E2  g and z must not change with the unit of Var1.
#
# OUTPUT
#   05_revision_R1/02_resultados/table_R1_04e_sim_MAD.csv
#   Run from C:/temp_paper.
# =====================================================================
library(parallel)
library(robustT2AFM)
library(MASS)
dir.create("05_revision_R1/02_resultados", recursive = TRUE, showWarnings = FALSE)

N_REP <- 2000; N_CAL <- 5000; BLOCK <- 200
K1 <- 30; K2 <- 100; I <- 20; J <- 4
RHO <- 0.6; H_MCD <- 0.67; ALPHA <- 0.001
OR <- 0.20; OS <- 4; SEED_BASE <- 2026
DELTA <- 1.0
SCALE_GRID <- c(1, 10, 100)
VARS <- paste0("Var", 1:J)
make_equi <- function(p, rho){ S <- matrix(rho, p, p); diag(S) <- 1; S }
Sigma_EQ <- make_equi(J, RHO)
celdas <- c("a_classical", "d_AFM_MCD_cov", "f_AFM_MCD_globalstd", "g_AFM_MCD_MAD", "z_AFM_MCD_zscore")

# --- Anchors: SCRIPT_R1_04b results -----------------------------------
ancla <- data.frame(
  sc     = c(1, 10, 100, 1, 10, 100),
  contam = c(0, 0, 0, 6, 6, 6),
  a = c(0.9010, 0.9010, 0.9010, 0.2700, 0.2700, 0.2700),
  d = c(0.9269, 0.9159, 0.9152, 0.8629, 0.8434, 0.8424),
  f = c(0.9269, 0.9269, 0.9269, 0.8626, 0.8626, 0.8626))

grid  <- expand.grid(sc = SCALE_GRID, ob = c(0, 6), KEEP.OUT.ATTRS = FALSE)
tasks <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
  data.frame(sc = grid$sc[i], ob = grid$ob[i],
             from = seq(1, N_REP, by = BLOCK), to = seq(BLOCK, N_REP, by = BLOCK))))
tasks <- split(tasks, seq_len(nrow(tasks)))

# =====================================================================
# One block of replicates for one configuration (runs on a worker)
# =====================================================================
run_block <- function(tk) {
  sc <- tk$sc; ob <- tk$ob
  esc <- c(sc, rep(1, J - 1))
  t2_manual <- function(Xb, mu, Si, n){ d <- colMeans(Xb) - mu; as.numeric(n * t(d) %*% Si %*% d) }
  rescale <- function(X){ X[, VARS] <- sweep(as.matrix(X[, VARS]), 2, esc, "*"); X }
  lam1 <- function(S) max(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
  comb <- function(Sl, w) Reduce(`+`, Map(function(S, wi) wi * S, Sl, w))
  w_div <- function(Sl, s){ Di <- diag(1 / s); l <- sapply(Sl, function(S) lam1(Di %*% S %*% Di)); (1 / l) / sum(1 / l) }

  refs <- function(f1){
    batches <- unique(f1$Batch)
    calR   <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
    S_mcd  <- calR$mcd_covariances
    mu_mcd <- calR$mu_r
    S_cls  <- lapply(batches, function(b) cov(f1[f1$Batch == b, VARS]))
    mu_cls <- colMeans(do.call(rbind, lapply(batches, function(b) colMeans(f1[f1$Batch == b, VARS]))))
    w_std  <- w_div(S_mcd, apply(sapply(S_mcd, function(S) sqrt(diag(S))), 1, median))
    w_mad  <- w_div(S_mcd, apply(as.matrix(f1[, VARS]), 2, mad))
    w_z    <- w_div(S_mcd, apply(as.matrix(f1[, VARS]), 2, sd))
    list(
      R = list(a = list(mu = mu_cls, Si = solve(comb(S_cls, rep(1 / length(batches), length(batches))))),
               d = list(mu = mu_mcd, Si = solve(calR$Sw)),
               f = list(mu = mu_mcd, Si = solve(comb(S_mcd, w_std))),
               g = list(mu = mu_mcd, Si = solve(comb(S_mcd, w_mad))),
               z = list(mu = mu_mcd, Si = solve(comb(S_mcd, w_z)))),
      w_cov = calR$weights,
      w_std = setNames(w_std, names(calR$weights)),
      w_mad = setNames(w_mad, names(calR$weights)),
      w_z   = setNames(w_z,   names(calR$weights)))
  }

  reps <- tk$from:tk$to
  out <- matrix(NA_real_, length(reps), 9,
                dimnames = list(NULL, c(celdas, "ratio_cov", "ratio_std", "ratio_mad", "ratio_z")))
  for (i in seq_along(reps)) {
    rep <- reps[i]
    set.seed(SEED_BASE * 200 + rep)
    sim <- simulate_batch_process(K1 = K1, K2 = 1, I = I, J = J, rho = RHO, Sigma = Sigma_EQ,
                                  outlier_batches_F1 = ob, outlier_rate = OR,
                                  outlier_shift = OS, prop_ooc_F2 = 0)
    f1 <- subset(sim, Phase == "Phase 1")
    cont <- if (ob > 0) unique(as.character(f1$Batch[as.character(f1$ContaminationType) != "Clean"])) else character(0)
    f1 <- rescale(f1)
    Rf <- refs(f1); R <- Rf$R

    if (ob > 0) {
      ic <- names(Rf$w_cov) %in% cont
      out[i, "ratio_cov"] <- mean(Rf$w_cov[!ic]) / mean(Rf$w_cov[ic])
      out[i, "ratio_std"] <- mean(Rf$w_std[!ic]) / mean(Rf$w_std[ic])
      out[i, "ratio_mad"] <- mean(Rf$w_mad[!ic]) / mean(Rf$w_mad[ic])
      out[i, "ratio_z"]   <- mean(Rf$w_z[!ic])   / mean(Rf$w_z[ic])
    }

    set.seed(SEED_BASE * 250 + rep)
    f1c <- lapply(seq_len(N_CAL), function(k){ X <- MASS::mvrnorm(I, rep(0, J), Sigma_EQ); colnames(X) <- VARS
                                             sweep(X, 2, esc, "*") })
    q <- sapply(names(R), function(nm) quantile(sapply(f1c, t2_manual, R[[nm]]$mu, R[[nm]]$Si, I), 1 - ALPHA))

    set.seed(SEED_BASE * 300 + rep)
    f2 <- lapply(seq_len(K2), function(k){ X <- MASS::mvrnorm(I, rep(DELTA, J), Sigma_EQ); colnames(X) <- VARS
                                          sweep(X, 2, esc, "*") })
    out[i, 1:5] <- sapply(seq_along(R), function(j) mean(sapply(f2, t2_manual, R[[j]]$mu, R[[j]]$Si, I) > q[j]))
  }
  data.frame(sc = sc, ob = ob, out, check.names = FALSE)
}

# =====================================================================
# Cluster
# =====================================================================
t0 <- Sys.time()
n_workers <- 10
cl <- makeCluster(n_workers)
clusterEvalQ(cl, { library(robustT2AFM); library(MASS); NULL })
clusterExport(cl, c("N_REP","N_CAL","K1","K2","I","J","RHO","H_MCD","ALPHA","OR","OS",
                    "SEED_BASE","DELTA","VARS","Sigma_EQ","celdas"))
bloques <- parLapplyLB(cl, tasks, run_block)
stopCluster(cl)
todo <- do.call(rbind, bloques)

# =====================================================================
# Assemble, anchors, print
# =====================================================================
se <- function(x) sd(x) / sqrt(length(x))
res <- do.call(rbind, lapply(split(todo, list(todo$sc, todo$ob), drop = TRUE), function(x)
  data.frame(scale_Var1 = x$sc[1], contam = x$ob[1],
             t(setNames(round(colMeans(x[, celdas]), 4), paste0("TPR_", celdas))),
             t(setNames(round(apply(x[, celdas], 2, se), 4), paste0("SE_", celdas))),
             ratio_w_cov = round(mean(x$ratio_cov), 3), ratio_w_std = round(mean(x$ratio_std), 3),
             ratio_w_mad = round(mean(x$ratio_mad), 3), ratio_w_z = round(mean(x$ratio_z), 3),
             n_rep = nrow(x), check.names = FALSE)))
res <- res[order(res$contam, res$scale_Var1), ]; rownames(res) <- NULL

ok_all <- all(res$n_rep == N_REP)
cat("\n== E1: a, d, f must reproduce SCRIPT_R1_04b ==\n")
for (i in seq_len(nrow(ancla))) {
  r <- res[res$scale_Var1 == ancla$sc[i] & res$contam == ancla$contam[i], ]
  oks <- c(abs(r$TPR_a_classical         - ancla$a[i]) < 5e-4,
           abs(r$TPR_d_AFM_MCD_cov       - ancla$d[i]) < 5e-4,
           abs(r$TPR_f_AFM_MCD_globalstd - ancla$f[i]) < 5e-4)
  ok_all <- ok_all && all(oks)
  cat(sprintf("  c=%3d contam=%d | a %s | d %s | f %s\n", ancla$sc[i], ancla$contam[i],
              ifelse(oks[1],"OK","FAIL"), ifelse(oks[2],"OK","FAIL"), ifelse(oks[3],"OK","FAIL")))
}
cat("\n== E2: g (MAD) and z (z-score) must not change with the unit of Var1 ==\n")
for (ob in c(0, 6)) {
  rg <- diff(range(res$TPR_g_AFM_MCD_MAD[res$contam == ob]))
  rz <- diff(range(res$TPR_z_AFM_MCD_zscore[res$contam == ob]))
  ok_all <- ok_all && rg < 5e-4 && rz < 5e-4
  cat(sprintf("  contam=%d | range over c: (g) %.4f %s | (z) %.4f %s\n", ob, rg, ifelse(rg < 5e-4, "OK", "FAIL"),
              rz, ifelse(rz < 5e-4, "OK", "FAIL")))
}

cat("\n== RESULT: TPR at equal ARL0 by unit of Var1 ==\n")
print(res[, c("scale_Var1","contam","TPR_a_classical","TPR_d_AFM_MCD_cov","TPR_f_AFM_MCD_globalstd",
              "TPR_g_AFM_MCD_MAD","TPR_z_AFM_MCD_zscore","ratio_w_cov","ratio_w_std","ratio_w_mad","ratio_w_z")], row.names = FALSE)

write.csv(res, "05_revision_R1/02_resultados/table_R1_04e_sim_MAD.csv", row.names = FALSE)
cat(sprintf("\n===== COMPLETE | %.1f min | %d workers | %d tasks =====\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), n_workers, length(tasks)))
cat(if (ok_all) "ALL anchors OK: valid run.\n" else "ANCHOR FAILURES: do NOT use these results.\n")
cat("Saved to 05_revision_R1/02_resultados/table_R1_04e_sim_MAD.csv\n")
