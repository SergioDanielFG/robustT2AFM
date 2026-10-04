# =====================================================================
# SCRIPT_R1_04_scale_dependence.R
# ---------------------------------------------------------------------
# Revision round 1, Reviewer 1, comment 4 (and Reviewer 3, note 1):
# scale dependence of the AFM weighting. Simulation part.
#
# WHAT IT COMPUTES
#   The power study of Table 5 (K = 30, I = 20, J = 4, rho = 0.6, delta =
#   1.0, 0 and 6 contaminated batches) when the first variable is
#   expressed in different units: multiplied by c = 1, 10 and 100. The
#   same simulated data are used for every c; only the unit of Var1
#   changes, so any difference between the c values is due to the
#   scale alone. Three charts are compared at a common ARL0:
#     (a) classical T2                       affine equivariant
#     (d) AFM-MCD, weights from lambda1 of the MCD COVARIANCE (published)
#     (e) AFM-MCD, weights from lambda1 of the MCD CORRELATION matrix
#         (the robust correlation-based scheme suggested by the reviewer;
#          Sw is still the weighted sum of the MCD covariances, so only
#          the weights change)
#   For (d) and (e) it also reports the healthy / contaminated mean
#   weight ratio, which shows whether the weights still identify the
#   contaminated batches when Var1 dominates the first eigenvalue.
#
# DESIGN DECISIONS
#   - Data, seeds (SEED_BASE*200 / *250 / *300 + rep), ARL0 equalisation
#     and the Phase 2 shift are those of SCRIPT_FINAL_ABLATIONS_2x2.R.
#     The shift is delta standard deviations in every variable, so in
#     the rescaled data Var1 is shifted by delta*c.
#   - The reference center is the published one (average of the MCD
#     centers), so the c = 1 rows anchor to Table 5. The center decision
#     (comment 7) is independent of the weighting decision studied here.
#   - The Var1 rescaling is applied after generation, to Phase 1, the
#     calibration set and Phase 2 alike.
#
# ANCHORS (tolerance 5e-04)
#   A1  c = 1: (a) 0.9010 / 0.2700 and (d) 0.9269 / 0.8629 (Table 5).
#   A2  (a) is affine equivariant: its TPR must be identical for
#       c = 1, 10, 100 within each contamination level.
#   A3  (e) is scale invariant by construction: same requirement.
#   Only (d) is allowed to change with c. The run is valid if A1-A3 hold.
#
# OUTPUT
#   05_revision_R1/02_resultados/table_R1_04_scale_dependence.csv
# =====================================================================
library(parallel)
dir.create("05_revision_R1/02_resultados", recursive = TRUE, showWarnings = FALSE)

# --- Global parameters (identical to the published campaign) ---------
N_REP <- 2000; N_CAL <- 5000
K1 <- 30; K2 <- 100; I <- 20; J <- 4
RHO <- 0.6; H_MCD <- 0.67; ALPHA <- 0.001
OR <- 0.20; OS <- 4; SEED_BASE <- 2026
DELTA <- 1.0
SCALE_GRID <- c(1, 10, 100)          # unit factor applied to Var1
VARS <- paste0("Var", 1:J)
make_equi <- function(p, rho){ S <- matrix(rho, p, p); diag(S) <- 1; S }
Sigma_EQ <- make_equi(J, RHO)
celdas <- c("a_classical", "d_AFM_MCD_cov", "e_AFM_MCD_cor")

# --- Anchors from Table 5 (delta = 1.0) ------------------------------
ancla <- data.frame(contam = c(0, 6), a = c(0.9010, 0.2700), d = c(0.9269, 0.8629))

grid <- expand.grid(sc = SCALE_GRID, ob = c(0, 6), KEEP.OUT.ATTRS = FALSE)
configs <- split(grid, seq_len(nrow(grid)))

# =====================================================================
# One full configuration (runs on a worker)
# =====================================================================
run_config <- function(cfg) {
  sc <- cfg$sc; ob <- cfg$ob
  esc <- c(sc, rep(1, J - 1))                     # unit factor per variable

  t2_manual <- function(Xb, mu, Si, n){ d <- colMeans(Xb) - mu; as.numeric(n * t(d) %*% Si %*% d) }
  se <- function(x) sd(x) / sqrt(length(x))
  rescale <- function(X){ X[, VARS] <- sweep(as.matrix(X[, VARS]), 2, esc, "*"); X }
  lam1 <- function(S) max(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
  comb <- function(Sl, w) Reduce(`+`, Map(function(S, wi) wi * S, Sl, w))

  refs <- function(f1){
    batches <- unique(f1$Batch)
    calR   <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
    S_mcd  <- calR$mcd_covariances
    mu_mcd <- calR$mu_r
    # (a) classical
    S_cls  <- lapply(batches, function(b) cov(f1[f1$Batch == b, VARS]))
    mu_cls <- colMeans(do.call(rbind, lapply(batches, function(b) colMeans(f1[f1$Batch == b, VARS]))))
    # (e) correlation-based weights on the same MCD covariances
    lam1_cor <- sapply(S_mcd, function(S) lam1(cov2cor(S)))
    w_cor    <- (1 / lam1_cor) / sum(1 / lam1_cor)
    list(
      R = list(a = list(mu = mu_cls, Si = solve(comb(S_cls, rep(1 / length(batches), length(batches))))),
               d = list(mu = mu_mcd, Si = solve(calR$Sw)),
               e = list(mu = mu_mcd, Si = solve(comb(S_mcd, w_cor)))),
      w_cov = calR$weights, w_cor = setNames(w_cor, names(calR$weights)))
  }

  tpr <- matrix(NA, N_REP, 3, dimnames = list(NULL, celdas))
  ratio_cov <- ratio_cor <- rep(NA_real_, N_REP)
  for (rep in seq_len(N_REP)) {
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
      ratio_cov[rep] <- mean(Rf$w_cov[!ic]) / mean(Rf$w_cov[ic])
      ratio_cor[rep] <- mean(Rf$w_cor[!ic]) / mean(Rf$w_cor[ic])
    }

    set.seed(SEED_BASE * 250 + rep)
    f1c <- lapply(seq_len(N_CAL), function(k){ X <- MASS::mvrnorm(I, rep(0, J), Sigma_EQ); colnames(X) <- VARS
                                             sweep(X, 2, esc, "*") })
    q <- sapply(names(R), function(nm) quantile(sapply(f1c, t2_manual, R[[nm]]$mu, R[[nm]]$Si, I), 1 - ALPHA))

    set.seed(SEED_BASE * 300 + rep)
    f2 <- lapply(seq_len(K2), function(k){ X <- MASS::mvrnorm(I, rep(DELTA, J), Sigma_EQ); colnames(X) <- VARS
                                          sweep(X, 2, esc, "*") })
    tpr[rep, ] <- sapply(seq_along(R), function(j) mean(sapply(f2, t2_manual, R[[j]]$mu, R[[j]]$Si, I) > q[j]))
  }

  data.frame(scale_Var1 = sc, contam = ob,
             t(round(colMeans(tpr), 4)), t(round(apply(tpr, 2, se), 4)),
             ratio_w_cov = round(mean(ratio_cov), 3), ratio_w_cor = round(mean(ratio_cor), 3),
             check.names = FALSE)
}

# =====================================================================
# Cluster
# =====================================================================
t0 <- Sys.time()
n_workers <- min(8, length(configs))
cl <- makeCluster(n_workers)
clusterEvalQ(cl, { library(robustT2AFM); library(MASS) })
clusterExport(cl, c("N_REP","N_CAL","K1","K2","I","J","RHO","H_MCD","ALPHA","OR","OS",
                    "SEED_BASE","DELTA","VARS","Sigma_EQ","celdas"))
resultados <- parLapplyLB(cl, configs, run_config)
stopCluster(cl)

# =====================================================================
# Assemble, anchors, print
# =====================================================================
res <- do.call(rbind, resultados)
names(res)[3:5] <- paste0("TPR_", celdas)
names(res)[6:8] <- paste0("SE_",  celdas)
res$diff_d_minus_a <- round(res$TPR_d_AFM_MCD_cov - res$TPR_a_classical, 4)
res$diff_e_minus_a <- round(res$TPR_e_AFM_MCD_cor - res$TPR_a_classical, 4)
res <- res[order(res$contam, res$scale_Var1), ]; rownames(res) <- NULL

ok_all <- TRUE
cat("\n== A1: c = 1 must reproduce Table 5 ==\n")
for (i in seq_len(nrow(ancla))) {
  r <- res[res$scale_Var1 == 1 & res$contam == ancla$contam[i], ]
  oks <- c(abs(r$TPR_a_classical - ancla$a[i]) < 5e-4, abs(r$TPR_d_AFM_MCD_cov - ancla$d[i]) < 5e-4)
  ok_all <- ok_all && all(oks)
  cat(sprintf("  contam=%d | a=%.4f (exp %.4f) %s | d=%.4f (exp %.4f) %s\n", ancla$contam[i],
              r$TPR_a_classical, ancla$a[i], ifelse(oks[1],"OK","FAIL"),
              r$TPR_d_AFM_MCD_cov, ancla$d[i], ifelse(oks[2],"OK","FAIL")))
}
cat("\n== A2 / A3: (a) and (e) must not change with the unit of Var1 ==\n")
for (ob in c(0, 6)) {
  r <- res[res$contam == ob, ]
  ra <- diff(range(r$TPR_a_classical)); re <- diff(range(r$TPR_e_AFM_MCD_cor)); rd <- diff(range(r$TPR_d_AFM_MCD_cov))
  ok_all <- ok_all && ra < 5e-4 && re < 5e-4
  cat(sprintf("  contam=%d | range over c: (a) %.4f %s | (e) %.4f %s | (d) %.4f (free)\n",
              ob, ra, ifelse(ra < 5e-4,"OK","FAIL"), re, ifelse(re < 5e-4,"OK","FAIL"), rd))
}

cat("\n== RESULT: TPR at equal ARL0 by unit of Var1 ==\n")
print(res[, c("scale_Var1","contam","TPR_a_classical","TPR_d_AFM_MCD_cov","TPR_e_AFM_MCD_cor",
              "diff_d_minus_a","diff_e_minus_a","ratio_w_cov","ratio_w_cor")], row.names = FALSE)

write.csv(res, "05_revision_R1/02_resultados/table_R1_04_scale_dependence.csv", row.names = FALSE)
cat(sprintf("\n===== COMPLETE | %.1f min | %d workers =====\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), n_workers))
cat(if (ok_all) "ALL anchors OK: valid run.\n" else "ANCHOR FAILURES: do NOT use these results.\n")
cat("Saved to 05_revision_R1/02_resultados/table_R1_04_scale_dependence.csv\n")
