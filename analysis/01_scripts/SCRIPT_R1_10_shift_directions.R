# =====================================================================
# SCRIPT_R1_10_shift_directions.R
# ---------------------------------------------------------------------
# Revision round 1, Reviewer 1, comment 10: Phase 2 mean shifts in more
# than one direction, at the same Mahalanobis distance.
#
# WHAT IT COMPUTES
#   The 2 x 2 decomposition of SCRIPT_FINAL_ABLATIONS_2x2.R, (a)
#   classical, (b) MCD only, (c) AFM only, (d) AFM-MCD, for four
#   directions of the Phase 2 shift, with 0 and 6 contaminated Phase 1
#   batches and delta = 1.0:
#     D1  (1, 1, 1, 1)   the direction of the manuscript
#     D2  (1,-1, 0, 0)   a contrast between two variables
#     D3  (1, 1,-1,-1)   a contrast between two pairs
#     D4  (1, 0, 0, 0)   a single variable
#   D2 to D4 are rescaled so that their Mahalanobis distance equals that
#   of rep(delta, J): D^2 = J / (1 + (J-1) rho) = 1.4286.
#
# DESIGN DECISIONS
#   - Everything except the Phase 2 mean is identical to the published
#     ablation: parameters, seeds (SEED_BASE*200 / *250 / *300 + rep),
#     ARL0 equalisation by empirical quantile, four variants.
#   - D1 uses the vector rep(delta, J) itself, so it reproduces Table 5.
#   - Phase 1 contamination is unchanged (+4 SD in all variables).
#   - Each replicate seeds its own stream.
#
# ANCHORS (Table 5, delta = 1.0; tolerance 5e-04)
#   D1, contam 0 : (a) 0.9010 (b) 0.8877 (c) 0.9215 (d) 0.9269
#   D1, contam 6 : (a) 0.2700 (b) 0.7304 (c) 0.5346 (d) 0.8629
#   The run is valid only if the eight anchors report OK.
#
# OUTPUT
#   05_revision_R1/02_resultados/table_R1_10_shift_directions.csv
# =====================================================================
library(parallel)
dir.create("05_revision_R1/02_resultados", recursive = TRUE, showWarnings = FALSE)

# --- Global parameters (identical to the published campaign) ---------
N_REP <- 2000; N_CAL <- 5000
K1 <- 30; K2 <- 100; I <- 20; J <- 4
RHO <- 0.6; H_MCD <- 0.67; ALPHA <- 0.001
OR <- 0.20; OS <- 4; SEED_BASE <- 2026
DELTAS <- c(1.0)
VARS <- paste0("Var", 1:J)
make_equi <- function(p, rho){ S <- matrix(rho, p, p); diag(S) <- 1; S }
Sigma_EQ <- make_equi(J, RHO)
Sigma_inv <- solve(Sigma_EQ)
celdas <- c("a_classical", "b_MCD_only", "c_AFM_only", "d_AFM_MCD")

# --- Shift directions, rescaled to a common Mahalanobis distance -----
maha2 <- function(v) as.numeric(t(v) %*% Sigma_inv %*% v)
raw_dirs <- list(
  D1_equal      = c(1,  1,  1,  1),
  D2_contrast_a = c(1, -1,  0,  0),
  D3_contrast_b = c(1,  1, -1, -1),
  D4_single_var = c(1,  0,  0,  0)
)
shift_vector <- function(dir_name, delta) {
  target <- maha2(rep(delta, J))
  if (dir_name == "D1_equal") return(rep(delta, J))
  v <- raw_dirs[[dir_name]]
  v * sqrt(target / maha2(v))
}

# --- Checks: D1 is the published shift, all directions share D^2 -----
cat("=== Shift vectors (same Mahalanobis distance within each delta) ===\n")
for (delta in DELTAS) {
  target <- maha2(rep(delta, J))
  cat(sprintf("  delta = %.1f  (D^2 = %.4f)\n", delta, target))
  for (nm in names(raw_dirs)) {
    v <- shift_vector(nm, delta)
    stopifnot(abs(maha2(v) - target) < 1e-10)
    cat(sprintf("    %-14s v = (%s)\n", nm,
                paste(sprintf("%+.4f", v), collapse = ", ")))
  }
  stopifnot(identical(shift_vector("D1_equal", delta), rep(delta, J)))
}
cat("\n")

# --- Anchors: Table 5 rows for the deltas in DELTAS -------------------
ancla <- data.frame(
  delta  = c(0.5, 0.5, 1.0, 1.0),
  contam = c(0,   6,   0,   6),
  a = c(0.1180, 0.0046, 0.9010, 0.2700),
  b = c(0.1118, 0.0535, 0.8877, 0.7304),
  c = c(0.1409, 0.0125, 0.9215, 0.5346),
  d = c(0.1546, 0.0997, 0.9269, 0.8629)
)
ancla <- ancla[ancla$delta %in% DELTAS, ]

# --- Configurations: 4 directions x length(DELTAS) x 2 contam --------
grid <- expand.grid(dir = names(raw_dirs), delta = DELTAS, ob = c(0, 6),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
configs <- split(grid, seq_len(nrow(grid)))

# =====================================================================
# FUNCTION THAT PROCESSES ONE FULL CONFIGURATION (runs on a worker)
# Same as run_config() in SCRIPT_FINAL_ABLATIONS_2x2.R; the Phase 2
# mean is v instead of rep(sh, J).
# =====================================================================
run_config <- function(cfg) {
  v <- shift_vector(cfg$dir, cfg$delta); ob <- cfg$ob

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
    f2 <- lapply(seq_len(K2), function(k){ X <- MASS::mvrnorm(I, v, Sigma_EQ); colnames(X) <- VARS; X })
    tpr[rep, ] <- sapply(seq_along(R), function(j)
      mean(sapply(f2, t2_manual, R[[j]]$mu, R[[j]]$Si, I) > q[j]))
  }

  data.frame(direction = cfg$dir, delta = cfg$delta, contam = ob,
             maha2 = round(maha2(v), 4),
             t(round(colMeans(tpr), 4)), t(round(apply(tpr, 2, se), 4)),
             check.names = FALSE, stringsAsFactors = FALSE)
}

# =====================================================================
# CLUSTER (8 workers, load-balanced over the configurations)
# =====================================================================
t0 <- Sys.time()
n_workers <- 8
cl <- makeCluster(n_workers)
clusterEvalQ(cl, { library(robustT2AFM); library(MASS) })
clusterExport(cl, c("N_REP","N_CAL","K1","K2","I","J","RHO","H_MCD",
                    "ALPHA","OR","OS","SEED_BASE","VARS","Sigma_EQ",
                    "Sigma_inv","celdas","raw_dirs","maha2","shift_vector"))
resultados <- parLapplyLB(cl, configs, run_config)
stopCluster(cl)

# =====================================================================
# ASSEMBLE, CHECK THE ANCHORS AND PRINT
# =====================================================================
res <- do.call(rbind, resultados)
names(res)[5:8]  <- paste0("TPR_", celdas)
names(res)[9:12] <- paste0("SE_",  celdas)
res$Advantage_d_minus_a <- round(res$TPR_d_AFM_MCD - res$TPR_a_classical, 4)
res <- res[order(res$delta, res$contam, match(res$direction, names(raw_dirs))), ]
rownames(res) <- NULL

cat("\n== ANCHOR: direction D1 must reproduce Table 5 ==\n")
todas_ok <- TRUE
for (i in seq_len(nrow(ancla))) {
  a <- ancla[i, ]
  r <- res[res$direction == "D1_equal" & res$delta == a$delta & res$contam == a$contam, ]
  oks <- c(abs(r$TPR_a_classical - a$a) < 5e-4, abs(r$TPR_b_MCD_only - a$b) < 5e-4,
           abs(r$TPR_c_AFM_only  - a$c) < 5e-4, abs(r$TPR_d_AFM_MCD  - a$d) < 5e-4)
  todas_ok <- todas_ok && all(oks)
  cat(sprintf("  delta=%.1f contam=%d | a=%.4f b=%.4f c=%.4f d=%.4f | anchors: %s\n",
              a$delta, a$contam, r$TPR_a_classical, r$TPR_b_MCD_only,
              r$TPR_c_AFM_only, r$TPR_d_AFM_MCD,
              paste(ifelse(oks, "OK", "FAIL"), collapse = " ")))
}

cat("\n== RESULT: TPR by shift direction (same Mahalanobis distance within each delta) ==\n")
print(res[, c("direction","delta","contam","maha2","TPR_a_classical","TPR_b_MCD_only",
              "TPR_c_AFM_only","TPR_d_AFM_MCD","Advantage_d_minus_a")],
      row.names = FALSE)

write.csv(res, "05_revision_R1/02_resultados/table_R1_10_shift_directions.csv", row.names = FALSE)
cat(sprintf("\n===== COMPLETE | %.1f min | %d workers =====\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), n_workers))
if (todas_ok) {
  cat(sprintf("ALL %d anchors OK: valid run. Directions D2-D4 are final.\n", 4 * nrow(ancla)))
} else {
  cat("ANCHOR FAILURES: do NOT use these results. Keep the whole console output.\n")
}
cat("Saved to 05_revision_R1/02_resultados/table_R1_10_shift_directions.csv\n")
