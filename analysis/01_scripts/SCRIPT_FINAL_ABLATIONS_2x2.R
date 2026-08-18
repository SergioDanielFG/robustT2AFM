# =====================================================================
# SCRIPT_FINAL_ABLATIONS_2x2.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The factorial decomposition of the contribution of each component of
#   the method. It crosses the two decisions that define it -- whether
#   each batch covariance is estimated classically or by MCD, and
#   whether batches are combined with uniform weights or with the AFM
#   weighting -- giving four variants:
#     (a) classical + uniform  (b) MCD + uniform
#     (c) classical + AFM      (d) AFM-MCD, the proposed method
#   For each of the eight configurations (4 shifts x 2 contamination
#   levels) it reports the true positive rate of the four variants with
#   their standard errors.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - The four variants are evaluated on the SAME simulated data, with
#     the same seeds and the same ARL0 equalisation procedure as the
#     power study, so variants (a) and (d) reproduce its columns exactly
#     and the comparison is like for like.
#   - Each method is brought to a common false-alarm rate through the
#     empirical 1-alpha quantile of its own statistic over an
#     independent set of N_CAL in-control batches, the same set for all
#     four variants within each replicate.
#   - The AFM weights of variant (c) are recomputed from the CLASSICAL
#     batch covariances, so that variant isolates the weighting and does
#     not smuggle in the MCD.
#   - The OUTER loop is parallelised (the eight configurations), not the
#     inner one. Each replicate still seeds its own stream
#     (SEED_BASE*200 / *250 / *300), so the result is bit for bit
#     identical to the serial version; parallelism only reduces wall
#     clock time and changes no arithmetic.
#   - On Windows parLapply spawns new processes that do NOT inherit the
#     environment, so each worker loads its own packages and receives
#     everything it needs explicitly.
#
# ANCHORS AND THEIR EXPECTED VALUES
#   Variants (a) and (d) must reproduce the published power table in the
#   eight configurations, to better than 5e-04:
#     contam=0  shift 0.5/1.0/1.5/2.0 -> (a) 0.1180 0.9010 1.0000 1.0000
#                                        (d) 0.1546 0.9269 1.0000 1.0000
#     contam=6  shift 0.5/1.0/1.5/2.0 -> (a) 0.0046 0.2700 0.9611 1.0000
#                                        (d) 0.0997 0.8629 0.9997 1.0000
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   The run is valid only if the eight anchors report OK; variants (b)
#   and (c) are taken as final only in that case.
#
# OUTPUT
#   02_resultados/table_ablation_2x2_final.csv
# =====================================================================
library(parallel)
dir.create("02_resultados", showWarnings = FALSE)

# --- Global parameters (identical to the campaign) ---
N_REP <- 2000; N_CAL <- 5000
K1 <- 30; K2 <- 100; I <- 20; J <- 4
RHO <- 0.6; H_MCD <- 0.67; ALPHA <- 0.001
OR <- 0.20; OS <- 4; SEED_BASE <- 2026
VARS <- paste0("Var", 1:J)
make_equi <- function(p, rho){ S <- matrix(rho, p, p); diag(S) <- 1; S }
Sigma_EQ <- make_equi(J, RHO)
celdas <- c("a_classical", "b_MCD_only", "c_AFM_only", "d_AFM_MCD")

# --- ANCHOR: published power table (for verification only) ---
ancla <- data.frame(
  contam = c(0, 6, 0, 6, 0, 6, 0, 6),
  shift  = c(0.5, 0.5, 1.0, 1.0, 1.5, 1.5, 2.0, 2.0),
  Hot    = c(0.1180, 0.0046, 0.9010, 0.2700, 1.0000, 0.9611, 1.0000, 1.0000),
  Rob    = c(0.1546, 0.0997, 0.9269, 0.8629, 1.0000, 0.9997, 1.0000, 1.0000)
)

# --- The eight configurations (parallel work units) ---
grid <- expand.grid(sh = c(0.5, 1.0, 1.5, 2.0), ob = c(0, 6),
                    KEEP.OUT.ATTRS = FALSE)
configs <- split(grid, seq_len(nrow(grid)))

# =====================================================================
# FUNCTION THAT PROCESSES ONE FULL CONFIGURATION (runs on a worker)
# It is EXACTLY the inner loop of the serial script, unchanged.
# =====================================================================
run_config <- function(cfg) {
  sh <- cfg$sh; ob <- cfg$ob

  t2_manual <- function(Xb, mu, Si, n){ d <- colMeans(Xb) - mu; as.numeric(n * t(d) %*% Si %*% d) }
  se <- function(x) sd(x) / sqrt(length(x))

  refs_4celdas <- function(f1){
    batches <- unique(f1$Batch)
    calR   <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
    S_mcd  <- calR$mcd_covariances
    mu_mcd <- calR$mu_r
    S_cls  <- lapply(batches, function(b) cov(f1[f1$Batch == b, VARS]))
    mu_cls <- colMeans(do.call(rbind, lapply(batches, function(b)
      colMeans(f1[f1$Batch == b, VARS]))))
    lam1_cls  <- sapply(S_cls, function(S) max(eigen(S, symmetric = TRUE, only.values = TRUE)$values))
    w_afm_cls <- (1 / lam1_cls) / sum(1 / lam1_cls)
    unif <- rep(1 / length(batches), length(batches))
    comb <- function(Sl, w) Reduce(`+`, Map(function(S, wi) wi * S, Sl, w))
    list(
      a = list(mu = mu_cls, Si = solve(comb(S_cls, unif))),
      b = list(mu = mu_mcd, Si = solve(comb(S_mcd, unif))),
      c = list(mu = mu_cls, Si = solve(comb(S_cls, w_afm_cls))),
      d = list(mu = mu_mcd, Si = solve(calR$Sw))
    )
  }

  tpr <- matrix(NA, N_REP, 4, dimnames = list(NULL, celdas))
  for (rep in seq_len(N_REP)) {
    set.seed(SEED_BASE * 200 + rep)
    sim <- simulate_batch_process(K1 = K1, K2 = 1, I = I, J = J, rho = RHO, Sigma = Sigma_EQ,
                                  outlier_batches_F1 = ob, outlier_rate = OR,
                                  outlier_shift = OS, prop_ooc_F2 = 0)
    f1 <- subset(sim, Phase == "Phase 1")
    R <- refs_4celdas(f1)

    set.seed(SEED_BASE * 250 + rep)
    f1c <- lapply(seq_len(N_CAL), function(k){ X <- MASS::mvrnorm(I, rep(0, J), Sigma_EQ); colnames(X) <- VARS; X })
    q <- sapply(names(R), function(nm)
      quantile(sapply(f1c, t2_manual, R[[nm]]$mu, R[[nm]]$Si, I), 1 - ALPHA))

    set.seed(SEED_BASE * 300 + rep)
    f2 <- lapply(seq_len(K2), function(k){ X <- MASS::mvrnorm(I, rep(sh, J), Sigma_EQ); colnames(X) <- VARS; X })
    tpr[rep, ] <- sapply(seq_along(R), function(j)
      mean(sapply(f2, t2_manual, R[[j]]$mu, R[[j]]$Si, I) > q[j]))
  }

  data.frame(contam = ob, shift = sh,
             t(round(colMeans(tpr), 4)), t(round(apply(tpr, 2, se), 4)),
             check.names = FALSE)
}

# =====================================================================
# START THE CLUSTER (8 workers; 4 cores left free to use the machine)
# =====================================================================
t0 <- Sys.time()
n_workers <- min(8, length(configs))
cl <- makeCluster(n_workers)

# Each worker (a new R process) loads its own packages...
clusterEvalQ(cl, { library(robustT2AFM); library(MASS) })
# ...and receives the globals and functions run_config uses.
clusterExport(cl, c("N_REP","N_CAL","K1","K2","I","J","RHO","H_MCD",
                    "ALPHA","OR","OS","SEED_BASE","VARS","Sigma_EQ","celdas"))

# Distribute the eight configurations across the workers.
resultados <- parLapply(cl, configs, run_config)
stopCluster(cl)

# =====================================================================
# ASSEMBLE, CHECK THE ANCHOR AND PRINT
# =====================================================================
res <- do.call(rbind, resultados)
names(res)[3:6]  <- paste0("TPR_", celdas)
names(res)[7:10] <- paste0("SE_",  celdas)
res <- res[order(res$shift, res$contam), ]
rownames(res) <- NULL

cat("\n== Check against the published power table (variants a and d) ==\n")
todas_ok <- TRUE
for (i in seq_len(nrow(res))) {
  ref  <- ancla[ancla$contam == res$contam[i] & ancla$shift == res$shift[i], ]
  ok_a <- abs(res$TPR_a_classical[i] - ref$Hot) < 5e-4
  ok_d <- abs(res$TPR_d_AFM_MCD[i] - ref$Rob) < 5e-4
  todas_ok <- todas_ok && ok_a && ok_d
  cat(sprintf("contam=%d sh=%.1f | a=%.4f b=%.4f c=%.4f d=%.4f | anchor a:%s d:%s\n",
              res$contam[i], res$shift[i], res$TPR_a_classical[i], res$TPR_b_MCD_only[i],
              res$TPR_c_AFM_only[i], res$TPR_d_AFM_MCD[i],
              ifelse(ok_a, "OK", sprintf("FAIL(exp %.4f)", ref$Hot)),
              ifelse(ok_d, "OK", sprintf("FAIL(exp %.4f)", ref$Rob))))
}

write.csv(res, "02_resultados/table_ablation_2x2_final.csv", row.names = FALSE)
cat(sprintf("\n===== COMPLETE | %.1f min | %d workers =====\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), n_workers))
if (todas_ok) {
  cat("ALL anchors OK: valid run. Variants (b) and (c) are final.\n")
} else {
  cat("ANCHOR FAILURES: do NOT use these results. Keep the whole console output.\n")
}
cat("Saved to 02_resultados/table_ablation_2x2_final.csv\n")
