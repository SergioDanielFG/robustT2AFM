# =====================================================================
# SCRIPT_M4e_TEP_table12_mstar14.R
# ---------------------------------------------------------------------
# Revision round 1, P20 (Section 4.7, Table 12). Reviewer 1, comments 6
# and 7; Omar, points 6 and 7.
#
# WHY
#   Table 12 compares center estimators keeping the weights, S_w and the
#   limit of the basic combination fixed. In table_M4_final_E_center.csv
#   those four rows (chart "V7") were counted with the published limit,
#   Eq. (8) with m* = 13 (UCL 19.69285). The paper now uses m* = 14
#   (UCL 19.64480, Section 2.5). This script recounts the same four rows
#   with the m* = 14 limit, so that every row of Table 12 uses the limit
#   of the paper. NOTHING IS SIMULATED.
#
# DESIGN (copied from SCRIPT_M4_TEP_final_mstar14.R, part E)
#   Phase 1: 24 normal runs + 6 IDV 1 runs, step 30, I = 20;
#   set.seed(2026) before the calibration; MCD of the centers with
#   set.seed(2027) and alpha 0.5; spatial median by Weiszfeld.
#   Phase 2: normal runs 25-34 + IDV 7 runs 7-26.
#   Drift = largest |center - reference| / sd over the four variables,
#   the reference being the same estimator on the 24 healthy batches.
#
# ANCHORS (STOP if any fails)
#   E1  with UCL 19.69285 the four rows reproduce table_M4_final_E_center
#       (phase1 IDV 1, phase2 IDV 7, chart V7): detected, false alarms,
#       median T2 of normal batches, lowest faulty T2 and drift
#   E2  the m* = 14 limit of the basic combination equals the limit of
#       the final method, 19.64480 (it does not depend on the weights)
#
# OUTPUT
#   05_revision_R1/02_resultados/table_M4e_table12.csv
#   Run from C:/temp_paper. Expected time: under one minute.
# =====================================================================
library(robustT2AFM)
library(robustbase)
stopifnot(packageVersion("robustT2AFM") >= "0.3.0")
DIR_DATOS <- "04_datos"
DIR_OUT   <- "05_revision_R1/02_resultados"
VARS   <- c("xmv_9", "xmv_8", "xmeas_19", "xmeas_17")
I_LOTE <- 20; PASO <- 30; J <- 4
ALPHA  <- 0.001; H_MCD <- 0.67; ALPHA_CEN <- 0.5
DESDE_SANO <- 1; DESDE_FALLO <- 161
SEED_CAL <- 2026; SEED_CEN <- 2027

if (!exists("ff_test")) { cat("Loading FaultFree_Testing...\n"); ff_test <- read.csv(file.path(DIR_DATOS, "TEP_FaultFree_Testing.csv")) }
if (!exists("fy"))      { cat("Loading Faulty_Testing...\n");    fy      <- read.csv(file.path(DIR_DATOS, "TEP_Faulty_Testing.csv")) }
TE <- read.csv(file.path(DIR_OUT, "table_M4_final_E_center.csv"), stringsAsFactors = FALSE)

# --- Helpers (copied from M4) ------------------------------------------
lote <- function(dat, run, fault, desde, nom) {
  sel <- dat$simulationRun == run & dat$faultNumber == fault &
         dat$sample %in% seq(desde, by = PASO, length.out = I_LOTE)
  s <- dat[sel, VARS]
  if (nrow(s) != I_LOTE) return(NULL)
  data.frame(Batch = nom, s, row.names = NULL)
}
fase1 <- function(f) do.call(rbind, c(
  lapply(1:24, function(r) lote(ff_test, r, 0, DESDE_SANO,  sprintf("F1_sano_%02d",  r))),
  lapply(1:6,  function(r) lote(fy,      r, f, DESDE_FALLO, sprintf("F1_fallo_%02d", r)))))
fase2 <- function(f) do.call(rbind, c(
  lapply(25:34, function(r) lote(ff_test, r, 0, DESDE_SANO,  sprintf("F2_sano_%02d",  r))),
  lapply(7:26,  function(r) lote(fy,      r, f, DESDE_FALLO, sprintf("F2_fallo_%03d", r)))))
spatial_median <- function(X, tol = 1e-10, maxit = 2000) {
  m <- colMeans(X)
  for (it in seq_len(maxit)) {
    d <- sqrt(rowSums((X - matrix(m, nrow(X), ncol(X), byrow = TRUE))^2)); d[d < 1e-12] <- 1e-12
    w <- 1 / d; m_new <- colSums(X * w) / sum(w)
    if (max(abs(m_new - m)) < tol) return(m_new)
    m <- m_new
  }
  warning("Weiszfeld did not converge"); m
}
mcd_center <- function(C) { set.seed(SEED_CEN); covMcd(C, alpha = ALPHA_CEN)$center }
t2_batches <- function(p2, mu, Sw) {
  Si <- solve(Sw); b <- unique(p2$Batch)
  setNames(sapply(b, function(bb) { d <- colMeans(p2[p2$Batch == bb, VARS]) - mu
                                     as.numeric(I_LOTE * t(d) %*% Si %*% d) }), b)
}
resumen <- function(t2, ucl) {
  ef <- grepl("fallo", names(t2))
  c(det = sum(t2[ef] > ucl), fa = sum(t2[!ef] > ucl),
    T2_med_sano = round(median(t2[!ef]), 1), T2_min_fallo = round(min(t2[ef]), 1),
    T2_max_sano = max(t2[!ef]), faulty_in_window = sum(t2[ef] > 19.6448 & t2[ef] <= 19.69285),
    normal_in_window = sum(t2[!ef] > 19.6448 & t2[!ef] <= 19.69285))
}
cal_v7 <- function(f1) {
  set.seed(SEED_CAL)
  k <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD, scaling = "none", center = "mean"))
  list(k = k, C = do.call(rbind, k$mcd_centers)[, VARS, drop = FALSE])
}
drift <- function(mu, ref, sd) round(max(abs(mu - ref) / sd), 4)

# --- Calibrations --------------------------------------------------------
f1  <- fase1(1); fs <- f1[grepl("sano", f1$Batch), ]
K   <- cal_v7(f1); Ks <- cal_v7(fs)
sdr <- apply(fs[, VARS], 2, sd)
sm_std <- function(C) spatial_median(sweep(C, 2, sdr, "/")) * sdr
cc <- list(mean   = list(mu = K$k$mu_r,              ref = Ks$k$mu_r),
           sm     = list(mu = spatial_median(K$C),   ref = spatial_median(Ks$C)),
           sm_std = list(mu = sm_std(K$C),           ref = sm_std(Ks$C)),
           mcd    = list(mu = mcd_center(K$C),       ref = mcd_center(Ks$C)))
ucl13 <- ucl_F_adjusted(K$k, I = I_LOTE, alpha = ALPHA, m_star = "nominal")$UCL
ucl14 <- ucl_F_adjusted(K$k, I = I_LOTE, alpha = ALPHA)$UCL
p2 <- fase2(7)

okE2 <- abs(ucl14 - 19.64480) < 1e-4 && abs(ucl13 - 19.69285) < 1e-4
cat(sprintf("E2 limits of the basic combination: m* = 13 -> %.5f | m* = 14 -> %.5f -> %s\n",
            ucl13, ucl14, if (okE2) "OK" else "*** FAIL ***"))
if (!okE2) stop("Anchor E2 failed. STOP.")

filas <- list(); okE1 <- TRUE
for (est in names(cc)) {
  t2 <- t2_batches(p2, cc[[est]]$mu, K$k$Sw)
  dr <- drift(cc[[est]]$mu, cc[[est]]$ref, sdr)
  r13 <- resumen(t2, ucl13); r14 <- resumen(t2, ucl14)
  ref <- TE[TE$chart == "V7" & TE$center == est & TE$phase1 == "IDV 1" & TE$phase2 == "IDV 7", ]
  ok <- nrow(ref) == 1 && r13["det"] == ref$det && r13["fa"] == ref$fa &&
        r13["T2_med_sano"] == ref$T2_med_sano && r13["T2_min_fallo"] == ref$T2_min_fallo &&
        abs(dr - ref$drift_SD) < 1e-4
  okE1 <- okE1 && ok
  cat(sprintf("E1 %-7s drift %.4f | m*=13: det %2d FA %d | m*=14: det %2d FA %d | med normal %.1f | min faulty %.1f | max normal %.2f | batches between the two limits: faulty %d, normal %d -> %s\n",
              est, dr, r13["det"], r13["fa"], r14["det"], r14["fa"], r14["T2_med_sano"], r14["T2_min_fallo"],
              r14["T2_max_sano"], r14["faulty_in_window"], r14["normal_in_window"], if (ok) "OK" else "*** DIFFERS ***"))
  filas[[length(filas) + 1]] <- data.frame(center = est, drift_SD = dr, UCL = ucl14,
    det = r14["det"], fa = r14["fa"], T2_med_sano = r14["T2_med_sano"], T2_min_fallo = r14["T2_min_fallo"],
    det_mstar13 = r13["det"], fa_mstar13 = r13["fa"], row.names = NULL)
}
if (!okE1) stop("Anchor E1 failed: the rows do not reproduce table_M4_final_E_center. STOP.")
tab <- do.call(rbind, filas)
write.csv(tab, file.path(DIR_OUT, "table_M4e_table12.csv"), row.names = FALSE)
cat("\nAnchors E1 and E2 OK\n"); print(tab, row.names = FALSE)
cat("Saved: 05_revision_R1/02_resultados/table_M4e_table12.csv\n")
