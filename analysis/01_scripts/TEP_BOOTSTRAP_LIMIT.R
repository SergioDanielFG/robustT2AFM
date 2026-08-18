# =====================================================================
# TEP_BOOTSTRAP_LIMIT.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   Two things on the real Tennessee Eastman Phase 1 and Phase 2.
#     (1) REPRODUCIBILITY. covMcd() is random: it explores subsets at
#         random and the base pipeline fixes no seed before calibrating,
#         so the published figures could move between sessions. The
#         calibration is repeated with N_SEEDS different seeds and the
#         spread of the UCL, the covariance inflation, the weight ratio,
#         the detection and the false alarms is measured.
#     (2) THE RESAMPLING LIMIT. The three bootstrap limits are computed
#         on the REAL Phase 1 and it is checked whether the conclusion
#         of the application changes when they are used instead of the
#         analytic limit. It also reports the overlap between the
#         classical and the MCD trimming, and the safe window between
#         the highest healthy T2 and the lowest faulty T2.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - Phase 1 and Phase 2 are rebuilt exactly as in
#     PHASE1_PIPELINE_STEP30.R, so the anchor block measures the same
#     configuration as the manuscript.
#   - The classical method uses no MCD and is therefore deterministic:
#     it is calibrated once and reused across seeds, so the whole spread
#     observed is attributable to the robust estimator.
#   - The MCD residuals are trimmed by ROBUST Mahalanobis distance, not
#     by the classical one; the overlap between the two trimmings is
#     reported so the difference is visible.
#   - B is set high (200 000) because the limit is computed on a single
#     data set, and the bootstrap is repeated over N_BSEEDS seeds so
#     that its own randomness is separated from the result.
#   - The weighted variant reuses the SAME bootstrap statistic as the
#     full-batch variant and only changes the quantile.
#
# ANCHORS AND THEIR EXPECTED VALUES
#   Before looking at anything new, block 1 must reproduce:
#     robust UCL            = 19.69   (tolerance 0.10)
#     classical UCL         = 19.46   (tolerance 0.10)
#     robust detection      = 20 of 20 (exact)
#     classical detection   =  2 of 20 (tolerance 1)
#     robust false alarms   =  0 of 10 (exact)
#     covariance inflation  = 16.14   (tolerance 1.50)
#     weight ratio          =  6.4    (tolerance 1.00)
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   If the anchor fails, STOP: nothing that follows is used. The
#   resampling limit is declared viable only if it gives the same Phase 2
#   result as the analytic one (20 detected, 0 false alarms); otherwise
#   the choice of limit matters and has to be argued explicitly in the
#   manuscript rather than left implicit.
#
# OUTPUT
#   02_resultados/tep_reproducibility.csv
#   02_resultados/tep_limits_comparison.csv
# =====================================================================
library(robustT2AFM)
library(robustbase)
DIR_DATOS <- "04_datos"
DIR_OUT   <- "02_resultados"

VARS   <- c("xmv_9", "xmv_8", "xmeas_19", "xmeas_17")
I_LOTE <- 20; PASO <- 30; J <- 4
ALPHA  <- 0.001; H_MCD <- 0.67; FALLO <- 7
DESDE_SANO <- 1; DESDE_FALLO <- 161
M_STAR  <- round(I_LOTE * H_MCD)
SEED_REF <- 2026                        # reference seed
N_SEEDS  <- 10                          # seeds for the stability test
B_BOOT   <- 200000                      # resamples (a single data set: keep it high)
N_BSEEDS <- 10                          # seeds for the bootstrap limit

muestras_esp <- function(d) seq(d, by = PASO, length.out = I_LOTE)

# ---------------------------------------------------------------------
# BLOCK 0 - Rebuilding the data (identical to the current pipeline)
# ---------------------------------------------------------------------
if (!exists("ff_test")) { cat("Loading FaultFree_Testing...\n"); ff_test <- read.csv(file.path(DIR_DATOS, "TEP_FaultFree_Testing.csv")) }
if (!exists("fy"))      { cat("Loading Faulty_Testing...\n");    fy      <- read.csv(file.path(DIR_DATOS, "TEP_Faulty_Testing.csv")) }

extraer_lote <- function(datos, run_id, fault_id, desde, batch_name) {
  sel <- datos$simulationRun == run_id & datos$faultNumber == fault_id &
         datos$sample %in% muestras_esp(desde)
  sub <- datos[sel, VARS]
  if (nrow(sub) != I_LOTE) return(NULL)
  data.frame(Batch = batch_name, sub, row.names = NULL)
}

f1 <- do.call(rbind, c(
  lapply(1:24, function(r) extraer_lote(ff_test, r, 0,     DESDE_SANO,  sprintf("F1_sano_%02d", r))),
  lapply(1:6,  function(r) extraer_lote(fy,      r, FALLO, DESDE_FALLO, sprintf("F1_fallo_%02d", r)))
))
f2 <- do.call(rbind, c(
  lapply(25:34, function(r) extraer_lote(ff_test, r, 0,     DESDE_SANO,  sprintf("F2_sano_%02d", r))),
  lapply(7:26,  function(r) extraer_lote(fy,      r, FALLO, DESDE_FALLO, sprintf("F2_fallo_%03d", r)))
))
f1$Batch <- as.character(f1$Batch); f2$Batch <- as.character(f2$Batch)
stopifnot(length(unique(f1$Batch)) == 30, length(unique(f2$Batch)) == 30)

b_f2 <- unique(f2$Batch)
es_fallo <- grepl("fallo", b_f2)
stopifnot(sum(es_fallo) == 20, sum(!es_fallo) == 10)

# The classical method uses no MCD: it is deterministic, calibrated once
cal_cls <- hotelling_classical_calibrate(f1, VARS)
ucl_cls <- hotelling_classical_ucl(K = 30, I = I_LOTE, J = J,
                                   alpha = ALPHA, phase = "II")$UCL
t2_cls  <- { m <- hotelling_classical_monitor(f2, cal_cls, VARS)
             m$T2[match(b_f2, m$Batch)] }

# ---------------------------------------------------------------------
# Metrics of a given robust calibration
# ---------------------------------------------------------------------
metricas <- function(cal) {
  ucl <- ucl_F_adjusted(cal, I = I_LOTE, alpha = ALPHA)$UCL
  m   <- monitor_afm_mcd(f2, cal, VARS)
  t2  <- m$T2[match(b_f2, m$Batch)]
  w   <- cal$weights
  cont <- grepl("fallo", names(w))
  ord  <- order(w)[1:6]
  list(ucl = ucl, t2 = t2,
       det = sum(t2[es_fallo] > ucl), fa = sum(t2[!es_fallo] > ucl),
       inflacion = det(cal_cls$Sp) / det(cal$Sw),
       ratio = mean(w[!cont]) / mean(w[cont]),
       seis_mas_bajos_son_contaminados = all(cont[ord]))
}

# ---------------------------------------------------------------------
# BLOCK 1 - ANCHOR with the reference seed
# ---------------------------------------------------------------------
set.seed(SEED_REF)
cal_ref <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
mr <- metricas(cal_ref)

cat("\n================ BLOCK 1 - ANCHOR OF THE APPLICATION ================\n")
chk <- function(nombre, obs, esp, tol) {
  cat(sprintf("  %-26s obs=%-9.2f exp=%-9.2f %s\n", nombre, obs, esp,
              if (abs(obs - esp) <= tol) "OK" else "REVIEW"))
  abs(obs - esp) <= tol
}
ok <- c(
  chk("Robust UCL",           mr$ucl,       19.69, 0.10),
  chk("Classical UCL",        ucl_cls,      19.46, 0.10),
  chk("Robust detection/20",  mr$det,       20,    0),
  chk("Classical detection/20", sum(t2_cls[es_fallo] > ucl_cls), 2, 1),
  chk("Robust false alarms/10", mr$fa,      0,     0),
  chk("Covariance inflation", mr$inflacion, 16.14, 1.50),
  chk("Weight ratio",         mr$ratio,     6.4,   1.00)
)
cat(sprintf("\n  ANCHOR: %s\n", if (all(ok)) "OK - it is safe to continue"
                               else "REVIEW - STOP HERE, do not use what follows"))

# ---------------------------------------------------------------------
# BLOCK 2 - REPRODUCIBILITY: how much everything moves with the seed
# ---------------------------------------------------------------------
cat("\n=========== BLOCK 2 - STABILITY WITH RESPECT TO THE SEED ===========\n")
cat(sprintf("  Recalibrating %d times with different seeds...\n", N_SEEDS))

rep_tab <- do.call(rbind, lapply(seq_len(N_SEEDS), function(s) {
  set.seed(1000 + s)
  cal <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
  m <- metricas(cal)
  data.frame(seed = 1000 + s, UCL_rob = m$ucl,
             det_rob = m$det, fa_rob = m$fa,
             inflation = m$inflacion, weight_ratio = m$ratio,
             six_lowest_ok = m$seis_mas_bajos_son_contaminados)
}))
print(format(rep_tab, digits = 4), row.names = FALSE)
write.csv(rep_tab, file.path(DIR_OUT, "tep_reproducibility.csv"), row.names = FALSE)

cat("\n  Spread summary:\n")
for (v in c("UCL_rob", "inflation", "weight_ratio")) {
  x <- rep_tab[[v]]
  cat(sprintf("    %-12s mean=%8.3f  sd=%7.4f  range=[%8.3f, %8.3f]  CV=%.3f%%\n",
              v, mean(x), sd(x), min(x), max(x), 100 * sd(x) / mean(x)))
}
cat(sprintf("    robust detection  : always %s\n",
            if (length(unique(rep_tab$det_rob)) == 1)
              paste0(unique(rep_tab$det_rob), "/20") else
              paste0("VARIES: ", paste(sort(unique(rep_tab$det_rob)), collapse = ", "))))
cat(sprintf("    robust false alarms: always %s\n",
            if (length(unique(rep_tab$fa_rob)) == 1)
              paste0(unique(rep_tab$fa_rob), "/10") else
              paste0("VARIES: ", paste(sort(unique(rep_tab$fa_rob)), collapse = ", "))))
cat(sprintf("    the 6 lowest weights are the contaminated ones in %d of %d runs\n",
            sum(rep_tab$six_lowest_ok), N_SEEDS))

# ---------------------------------------------------------------------
# BLOCK 3 - BOOTSTRAP LIMIT on the real Phase 1
# ---------------------------------------------------------------------
cat("\n============ BLOCK 3 - BOOTSTRAP LIMIT (real Phase 1) ============\n")

lotes <- lapply(split(f1[, VARS], f1$Batch), as.matrix)
nm    <- names(lotes)
stopifnot(all(vapply(lotes, nrow, integer(1)) == I_LOTE))
stopifnot(all(nm %in% names(cal_ref$mcd_centers)))

Si  <- solve(cal_ref$Sw)
w   <- as.numeric(cal_ref$weights[match(nm, names(cal_ref$weights))])
stopifnot(!anyNA(w))

# Residuals: MCD subset (trimmed by ROBUST distance) and full batch
res_mcd <- lapply(nm, function(k) {
  X <- lotes[[k]]
  S <- cal_ref$mcd_covariances[[k]]
  stopifnot(det(S) > 0)
  d  <- mahalanobis(X, cal_ref$mcd_centers[[k]], S)
  Xm <- X[order(d)[seq_len(M_STAR)], , drop = FALSE]
  sweep(Xm, 2, colMeans(Xm))
})
res_full <- lapply(nm, function(k) { X <- lotes[[k]]; sweep(X, 2, colMeans(X)) })
names(res_mcd) <- nm; names(res_full) <- nm

solap <- vapply(nm, function(k) {
  X <- lotes[[k]]
  dc <- mahalanobis(X, colMeans(X), cov(X))
  dr <- mahalanobis(X, cal_ref$mcd_centers[[k]], cal_ref$mcd_covariances[[k]])
  length(intersect(order(dc)[seq_len(M_STAR)], order(dr)[seq_len(M_STAR)]))
}, integer(1))
cat(sprintf("  Overlap between classical and MCD trimming: median %d of %d (range %d-%d)\n",
            median(solap), M_STAR, min(solap), max(solap)))

# Bootstrap T2, vectorised by batch
t2_boot <- function(res_list, B) {
  K <- length(res_list)
  pick <- sample.int(K, B, replace = TRUE)
  t2 <- numeric(B); wb <- numeric(B); pos <- 1L
  for (k in seq_len(K)) {
    nb <- sum(pick == k)
    if (nb == 0L) next
    R  <- res_list[[k]]; nr <- nrow(R)
    ii <- matrix(sample.int(nr, nb * I_LOTE, replace = TRUE), nrow = nb)
    Mm <- matrix(0, nb, J)
    for (j in seq_len(I_LOTE)) Mm <- Mm + R[ii[, j], , drop = FALSE]
    Mm <- Mm / I_LOTE
    t2[pos:(pos + nb - 1L)] <- I_LOTE * rowSums((Mm %*% Si) * Mm)
    wb[pos:(pos + nb - 1L)] <- w[k]
    pos <- pos + nb
  }
  list(t2 = t2, w = wb)
}
wtd_quantile <- function(x, ww, p) {
  o <- order(x); x <- x[o]; ww <- ww[o]
  x[which(cumsum(ww) / sum(ww) >= p)[1]]
}

cat(sprintf("  Computing limits with B=%s and %d seeds...\n",
            format(B_BOOT, big.mark = ","), N_BSEEDS))
lim <- do.call(rbind, lapply(seq_len(N_BSEEDS), function(s) {
  set.seed(5000 + s)
  bm <- t2_boot(res_mcd,  B_BOOT)
  bf <- t2_boot(res_full, B_BOOT)
  data.frame(seed = 5000 + s,
             bMCD  = quantile(bm$t2, 1 - ALPHA, names = FALSE),
             bFull = quantile(bf$t2, 1 - ALPHA, names = FALSE),
             bPond = wtd_quantile(bf$t2, bf$w, 1 - ALPHA))
}))
print(format(lim, digits = 5), row.names = FALSE)

cat("\n  Limits (mean over seeds):\n")
for (v in c("bMCD", "bFull", "bPond")) {
  x <- lim[[v]]
  cat(sprintf("    %-6s mean=%7.3f  sd=%6.4f  range=[%7.3f, %7.3f]\n",
              v, mean(x), sd(x), min(x), max(x)))
}
cat(sprintf("    Reference: analytic F = %.3f | chi2_%d(0.999) = %.3f\n",
            mr$ucl, J, qchisq(1 - ALPHA, J)))

# ---------------------------------------------------------------------
# BLOCK 4 - COMPARISON: does the conclusion change?
# ---------------------------------------------------------------------
cat("\n========== BLOCK 4 - PHASE 2 DETECTION WITH EACH LIMIT ==========\n")

limites <- c("F analitico (Ec. 8)" = mr$ucl,
             "Bootstrap MCD"        = mean(lim$bMCD),
             "Bootstrap completo"   = mean(lim$bFull),
             "Bootstrap ponderado"  = mean(lim$bPond),
             "chi2_4(0.999)"        = qchisq(1 - ALPHA, J))

comp <- do.call(rbind, lapply(names(limites), function(nombre) {
  u <- limites[[nombre]]
  data.frame(Limit = nombre, UCL = u,
             Detected_of_20 = sum(mr$t2[es_fallo]  > u),
             False_alarms_of_10 = sum(mr$t2[!es_fallo] > u),
             row.names = NULL)
}))
print(format(comp, digits = 4), row.names = FALSE)
write.csv(comp, file.path(DIR_OUT, "tep_limits_comparison.csv"), row.names = FALSE)

cat("\n  Robust T2 of the 10 healthy batches (sorted):\n    ")
cat(paste(round(sort(mr$t2[!es_fallo]), 2), collapse = "  "), "\n")
cat("  Lowest robust T2 among the 20 faulty batches: ",
    round(min(mr$t2[es_fallo]), 2), "\n")
cat("  -> any limit between those two values gives 20/0\n")
cat(sprintf("  -> safe window: (%.2f , %.2f)\n",
            max(mr$t2[!es_fallo]), min(mr$t2[es_fallo])))

# ---------------------------------------------------------------------
# BLOCK 5 - VERDICT
# ---------------------------------------------------------------------
cat("\n==================== BLOCK 5 - VERDICT ====================\n")
fila_b <- comp[comp$Limit == "Bootstrap MCD", ]
if (fila_b$Detected_of_20 == 20 && fila_b$False_alarms_of_10 == 0) {
  cat("  The bootstrap MCD limit gives the SAME result as the analytic one:\n")
  cat("  20 detected, 0 false alarms.\n\n")
  cat("  >>> The analytic limit can be kept as the operating limit, the\n")
  cat("      bootstrap reported as a better calibrated alternative, and the\n")
  cat("      conclusion of the application shown NOT to depend on which one\n")
  cat("      is chosen. Nothing has to be redone.\n")
} else {
  cat(sprintf("  The bootstrap MCD limit gives %d/20 detected and %d/10 false alarms,\n",
              fila_b$Detected_of_20, fila_b$False_alarms_of_10))
  cat("  different from the analytic one (20/0).\n\n")
  cat("  >>> THE CHOICE OF LIMIT DOES MATTER. It has to be decided explicitly\n")
  cat("      and argued in the manuscript, not left implicit.\n")
}
