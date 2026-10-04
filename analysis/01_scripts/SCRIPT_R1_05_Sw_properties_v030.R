# =====================================================================
# SCRIPT_R1_05_Sw_properties_v030.R
# ---------------------------------------------------------------------
# Revision round 1, Reviewer 1, comment 5: statistical properties of
# the weighted covariance estimator Sw.
#
# VERSION v030 (29-sep-2026): same design, seeds and anchor as
# SCRIPT_R1_05_Sw_properties.R, which ran on 25-sep with robustT2AFM
# 0.2.0 (weights from the raw first eigenvalues). This copy runs with
# 0.3.0, so Sw is the one of the final method (weights computed after
# the common MAD scaling). Estimator (c) uses the same scaling for its
# weights, so that (c) and (d) differ only in MCD vs classical
# covariances. The 0.2.0 results are kept; this writes *_v030.csv.
#
# WHAT IT COMPUTES
#   Monte Carlo study of four Phase 1 covariance estimators against the
#   true Sigma (equicorrelation, rho = 0.6, J = 4):
#     (a) S_p      classical pooled covariance
#     (b) MCD-unif uniform average of the batch MCD covariances
#     (c) AFM-cls  AFM weights on the classical batch covariances
#     (d) Sw       AFM weights on the MCD covariances (proposed)
#   For each estimator and replicate: eigenvalues, determinant,
#   Frobenius error ||S - Sigma||_F and the elementwise error S - Sigma.
#   Averages over N_REP replicates give the bias of each quantity.
#   Conditions: Phase 1 clean and contaminated (6 batches, 20% of
#   observations shifted 4 SD), and K = 30, 100, 300 batches with I = 20
#   fixed, to see whether the bias shrinks as K grows.
#
# DESIGN DECISIONS
#   - Phase 1 is generated exactly as in the power study (same
#     simulate_batch_process() call, seeds SEED_BASE*200 + rep), so the
#     K = 30 calibrations are the ones behind Tables 4 and 5.
#   - True values: eigenvalues 2.8, 0.4, 0.4, 0.4; det = 0.1792.
#   - Each replicate seeds its own stream; parallel and serial runs give
#     the same numbers.
#
# ANCHOR (theoretical)
#   With a clean Phase 1, S_p is exactly unbiased: E[S_p] = Sigma. Its
#   largest absolute elementwise bias must be below 3 Monte Carlo
#   standard errors in every K. The script stops otherwise.
#
# OUTPUT
#   05_revision_R1/02_resultados/table_R1_05_Sw_properties_v030.csv
#   05_revision_R1/02_resultados/table_R1_05_Sw_elementwise_bias_v030.csv
# =====================================================================
library(parallel)
library(robustT2AFM)
stopifnot(packageVersion("robustT2AFM") == "0.3.0")
dir.create("05_revision_R1/02_resultados", recursive = TRUE, showWarnings = FALSE)

N_REP <- 2000
K_GRID <- c(30, 100, 300)
I <- 20; J <- 4
RHO <- 0.6; H_MCD <- 0.67
OB <- 6; OR <- 0.20; OS <- 4
SEED_BASE <- 2026
N_WORKERS <- 8
VARS <- paste0("Var", 1:J)

make_equi <- function(p, rho){ S <- matrix(rho, p, p); diag(S) <- 1; S }
Sigma_EQ <- make_equi(J, RHO)
lam_true <- eigen(Sigma_EQ, symmetric = TRUE, only.values = TRUE)$values
det_true <- det(Sigma_EQ)
est_names <- c("a_Sp_classical", "b_MCD_uniform", "c_AFM_classical", "d_Sw_AFM_MCD")

# --- One replicate: the four estimators and their metrics ------------
one_rep <- function(rep, K, ob) {
  set.seed(SEED_BASE * 200 + rep)
  sim <- simulate_batch_process(K1 = K, K2 = 0, I = I, J = J, rho = RHO,
                                Sigma = Sigma_EQ, outlier_batches_F1 = ob,
                                outlier_rate = OR, outlier_shift = OS,
                                prop_ooc_F2 = 0)
  f1 <- subset(sim, Phase == "Phase 1")
  batches <- unique(f1$Batch)

  calR  <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
  S_mcd <- calR$mcd_covariances
  S_cls <- lapply(batches, function(b) cov(f1[f1$Batch == b, VARS]))
  Dinv <- diag(1 / calR$scale)
  lam1_cls  <- sapply(S_cls, function(S) max(eigen(Dinv %*% S %*% Dinv,
                      symmetric = TRUE, only.values = TRUE)$values))
  w_afm_cls <- (1 / lam1_cls) / sum(1 / lam1_cls)
  unif <- rep(1 / length(batches), length(batches))
  comb <- function(Sl, w) Reduce(`+`, Map(function(S, wi) wi * S, Sl, w))

  ests <- list(comb(S_cls, unif), comb(S_mcd, unif), comb(S_cls, w_afm_cls), calR$Sw)

  out <- lapply(ests, function(S) {
    E <- S - Sigma_EQ
    c(eigen(S, symmetric = TRUE, only.values = TRUE)$values,   # 4 eigenvalues
      det(S),
      sqrt(sum(E^2)),                                          # Frobenius error
      E[upper.tri(E, diag = TRUE)])                            # 10 unique elements
  })
  unlist(out)   # 4 estimators x 16 numbers
}

# --- Cluster -----------------------------------------------------------
cl <- makeCluster(N_WORKERS)
clusterEvalQ(cl, { library(robustT2AFM); library(MASS) })
clusterExport(cl, c("I","J","RHO","H_MCD","OR","OS","SEED_BASE","VARS",
                    "Sigma_EQ","one_rep"))

se <- function(x) sd(x) / sqrt(length(x))
ut <- which(upper.tri(Sigma_EQ, diag = TRUE), arr.ind = TRUE)
elem_lab <- paste0("s", ut[, 1], ut[, 2])

summary_rows <- list(); elem_rows <- list()
t_ini <- Sys.time()
for (K in K_GRID) for (ob in c(0, OB)) {
  cond <- if (ob == 0) "clean" else "contaminated"
  cat(sprintf("--- K = %d | Phase 1 %s ---\n", K, cond))
  t0 <- Sys.time()
  M <- parSapply(cl, seq_len(N_REP), one_rep, K = K, ob = ob)   # 64 x N_REP
  M <- t(M)
  for (e in seq_along(est_names)) {
    cols <- (e - 1) * 16 + 1:16
    B <- M[, cols, drop = FALSE]
    lam <- B[, 1:4]; dt <- B[, 5]; fro <- B[, 6]; el <- B[, 7:16]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      K = K, phase1 = cond, estimator = est_names[e],
      lambda1_mean = mean(lam[, 1]), lambda1_true = lam_true[1],
      lambda1_rel_bias = mean(lam[, 1]) / lam_true[1] - 1,
      lambda2_mean = mean(lam[, 2]), lambda3_mean = mean(lam[, 3]),
      lambda4_mean = mean(lam[, 4]), lambda_rest_true = lam_true[2],
      det_mean = mean(dt), det_true = det_true,
      det_rel_bias = mean(dt) / det_true - 1,
      frobenius_mean = mean(fro), frobenius_se = se(fro),
      max_abs_elem_bias = max(abs(colMeans(el))),
      max_elem_bias_se = max(apply(el, 2, se)),
      n_rep = N_REP)
    eb <- colMeans(el); es <- apply(el, 2, se)
    elem_rows[[length(elem_rows) + 1]] <- data.frame(
      K = K, phase1 = cond, estimator = est_names[e],
      element = elem_lab, bias = eb, se = es, z = eb / es, row.names = NULL)
  }
  cat(sprintf("    %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
stopCluster(cl)

res  <- do.call(rbind, summary_rows)
elem <- do.call(rbind, elem_rows)

# --- Anchor: classical S_p unbiased with a clean Phase 1 --------------
cat("\n== ANCHOR: elementwise bias of S_p, clean Phase 1 (must be < 3 SE) ==\n")
ok_all <- TRUE
for (K in K_GRID) {
  z <- elem[elem$K == K & elem$phase1 == "clean" & elem$estimator == "a_Sp_classical", ]
  zmax <- max(abs(z$z)); ok <- zmax < 3; ok_all <- ok_all && ok
  cat(sprintf("  K = %3d | max |bias / SE| = %.2f  %s\n", K, zmax, if (ok) "OK" else "FAIL"))
}

cat("\n== RESULT ==\n")
print(res[, c("K","phase1","estimator","lambda1_mean","lambda1_rel_bias",
              "lambda2_mean","det_rel_bias","frobenius_mean","max_abs_elem_bias")],
      row.names = FALSE, digits = 4)
cat(sprintf("\nTrue eigenvalues: %s | true det: %.4f\n",
            paste(round(lam_true, 4), collapse = ", "), det_true))

write.csv(res,  "05_revision_R1/02_resultados/table_R1_05_Sw_properties_v030.csv", row.names = FALSE)
write.csv(elem, "05_revision_R1/02_resultados/table_R1_05_Sw_elementwise_bias_v030.csv", row.names = FALSE)
cat(sprintf("\n===== COMPLETE | %.1f min | %d workers =====\n",
            as.numeric(difftime(Sys.time(), t_ini, units = "mins")), N_WORKERS))
cat(if (ok_all) "Anchor OK: valid run.\n" else "ANCHOR FAILED: do not use these results.\n")
cat("Saved to 05_revision_R1/02_resultados/table_R1_05_Sw_properties_v030.csv and _elementwise_bias_v030.csv\n")
