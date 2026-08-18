# =====================================================================
# 18_phase1_idv1_cross_matrix.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The full 2x2 matrix of calibrating Phase 1 with one Tennessee
#   Eastman fault and monitoring Phase 2 with each fault (IDV 7 and
#   IDV 1). For every cell it reports the robust UCL, the healthy /
#   contaminated weight ratio, the drift of the reference center mu_r,
#   the covariance-determinant inflation, the faulty batches detected
#   out of twenty and the false alarms out of ten, for both methods.
#
#   The matrix separates two competing explanations for the false
#   alarms observed when Phase 1 is contaminated with IDV 1:
#     (a) the CROSS itself: calibrating with one fault and monitoring
#         another misplaces the reference relative to what is watched;
#     (b) the CONTAMINATION itself: mu_r is a SIMPLE average of the
#         batch centers, with no weights, so it carries part of the
#         shift of the contaminated batches.
#   If the diagonal (7->7 and 1->1) is clean and only the crossed cells
#   fail, the answer is (a). If the IDV 1 row fails in both columns,
#   the answer is (b).
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - Batches are built by spaced subsampling, one observation every
#     PASO = 30, so that the observations of a batch can be treated as
#     independent in time.
#   - Phase 1 always holds 24 healthy batches plus 6 contaminated ones,
#     and Phase 2 always holds 10 healthy plus 20 faulty; only the
#     fault identity changes between cells.
#   - Faulty batches start at DESDE_FALLO, once the fault is already
#     present in the run.
#   - The Phase 2 runs do not reuse the Phase 1 runs.
#   - The center drift is measured against the calibration obtained
#     with the healthy batches only, on the scale of the healthy
#     observations' standard deviation.
#
# ANCHORS AND THEIR EXPECTED VALUES
#   Cell IDV 7 -> IDV 7 must reproduce the published application:
#     robust UCL          = 19.69285   (tolerance 1e-04)
#     weight ratio        =  6.4394    (tolerance 0.01)
#     determinant inflation = 16.1351  (tolerance 0.01)
#     detected            = 20 of 20   (exact)
#     false alarms        =  0 of 10   (exact)
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   The script STOPS if the anchor cell does not reproduce those five
#   values, because none of the remaining cells would be comparable.
#   The prediction registered before running is explanation (b): the
#   IDV 1 row fails in both columns. If (a) comes out instead, the
#   prediction was wrong and the reason has to be understood.
#
# OUTPUT
#   02_resultados/phase1_cross_matrix.csv
# =====================================================================

library(robustT2AFM)

VARS    <- c("xmv_9", "xmv_8", "xmeas_19", "xmeas_17")
I_LOTE  <- 20
PASO    <- 30
J       <- 4
ALPHA   <- 0.001
H_MCD   <- 0.67
DESDE_SANO  <- 1
DESDE_FALLO <- 161
N_SANOS_F1  <- 24
N_CONTA_F1  <- 6
N_SANOS_F2  <- 10
N_FALLO_F2  <- 20
DIR_DATOS <- "04_datos"
DIR_OUT   <- "02_resultados"

muestras_esp <- function(d) seq(d, by = PASO, length.out = I_LOTE)

cat("Loading TEP data...\n")
ff <- read.csv(file.path(DIR_DATOS, "TEP_FaultFree_Testing.csv"))
fy <- read.csv(file.path(DIR_DATOS, "TEP_Faulty_Testing.csv"))
cat(sprintf("A batch spans %d samples.\n\n", 1 + (I_LOTE - 1) * PASO))

lote <- function(dat, run, fault, desde, nom) {
  sel <- dat$simulationRun == run & dat$faultNumber == fault &
    dat$sample %in% muestras_esp(desde)
  s <- dat[sel, VARS]
  if (nrow(s) != I_LOTE) return(NULL)
  data.frame(Batch = nom, s, row.names = NULL)
}

# ---------------------------------------------------------------------
#  Calibration: Phase 1 with 24 healthy + 6 contaminated by the fault
# ---------------------------------------------------------------------
calibrar <- function(fallo_f1) {
  f1 <- do.call(rbind, c(
    lapply(1:N_SANOS_F1, function(r)
      lote(ff, r, 0, DESDE_SANO, sprintf("sano_%02d", r))),
    lapply(1:N_CONTA_F1, function(r)
      lote(fy, r, fallo_f1, DESDE_FALLO, sprintf("cont_%02d", r)))
  ))
  stopifnot(nrow(f1) == (N_SANOS_F1 + N_CONTA_F1) * I_LOTE)

  cal_r <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
  cal_c <- hotelling_classical_calibrate(f1, VARS)
  K     <- length(cal_r$weights)

  nm   <- names(cal_r$weights)
  cont <- grepl("cont", nm)
  w    <- as.numeric(cal_r$weights)
  seis <- names(sort(cal_r$weights))[seq_len(N_CONTA_F1)]

  # how far the center drifts from the one obtained with healthy batches only
  cal_s <- suppressMessages(
    calibrate_afm_mcd(f1[grepl("sano", f1$Batch), ], VARS, mcd_alpha = H_MCD)
  )
  sd_ref <- apply(f1[grepl("sano", f1$Batch), VARS], 2, sd)
  corri  <- max(abs(cal_r$mu_r - cal_s$mu_r) / sd_ref)

  # shift of the contaminated batches relative to the healthy center
  despl <- max(abs(colMeans(f1[grepl("cont", f1$Batch), VARS]) -
                   cal_s$mu_r) / sd_ref)

  list(cal_r = cal_r, cal_c = cal_c, K = K,
       ucl_r = ucl_F_adjusted(cal_r, I = I_LOTE, alpha = ALPHA)$UCL,
       ucl_c = hotelling_classical_ucl(K = K, I = I_LOTE, J = J,
                                       alpha = ALPHA, phase = "II")$UCL,
       ratio_w  = mean(w[!cont]) / mean(w[cont]),
       n_abajo  = sum(seis %in% nm[cont]),
       infl     = det(cal_c$Sp) / det(cal_r$Sw),
       despl    = despl,
       corrimiento_centro = corri)
}

# ---------------------------------------------------------------------
#  Monitoring: 10 healthy batches + 20 with the given fault
# ---------------------------------------------------------------------
monitorear <- function(cal, fallo_f2) {
  p2 <- rbind(
    do.call(rbind, lapply(25:(24 + N_SANOS_F2), function(r)
      lote(ff, r, 0, DESDE_SANO, sprintf("F2sano_%02d", r)))),
    do.call(rbind, lapply(7:(6 + N_FALLO_F2), function(r)
      lote(fy, r, fallo_f2, DESDE_FALLO, sprintf("F2fallo_%03d", r))))
  )
  stopifnot(nrow(p2) == (N_SANOS_F2 + N_FALLO_F2) * I_LOTE)

  mr <- monitor_afm_mcd(p2, cal$cal_r, VARS)
  mc <- hotelling_classical_monitor(p2, cal$cal_c, VARS)
  bf <- unique(p2$Batch); ef <- grepl("fallo", bf)
  tr <- mr$T2[match(bf, mr$Batch)]
  tc <- mc$T2[match(bf, mc$Batch)]

  list(det_rob = sum(tr[ef] > cal$ucl_r),
       det_cls = sum(tc[ef] > cal$ucl_c),
       fa_rob  = sum(tr[!ef] > cal$ucl_r),
       fa_cls  = sum(tc[!ef] > cal$ucl_c),
       t2_sano_med  = median(tr[!ef]),
       t2_fallo_med = median(tr[ef]))
}

# =====================================================================
#  The 2x2 matrix
# =====================================================================
cals <- list(`7` = calibrar(7), `1` = calibrar(1))

cat("=== Anchor: the IDV 7 -> IDV 7 cell ===\n")
a <- cals[["7"]]
m <- monitorear(a, 7)
cat(sprintf("  UCL %.5f (exp. 19.69285) | weight ratio %.4f (exp. 6.4394)\n",
            a$ucl_r, a$ratio_w))
cat(sprintf("  inflation %.4f (exp. 16.1351) | %d of 20 detected, %d false alarms\n",
            a$infl, m$det_rob, m$fa_rob))
ok <- abs(a$ucl_r - 19.69285) < 1e-4 && abs(a$ratio_w - 6.4394) < 0.01 &&
      abs(a$infl - 16.1351) < 0.01 && m$det_rob == 20 && m$fa_rob == 0
cat(sprintf("  ANCHOR: %s\n\n", if (ok) "OK" else "*** FAILS, STOP ***"))
stopifnot(ok)

cat("=== Properties of each calibration ===\n")
for (f in c("7", "1")) {
  c_ <- cals[[f]]
  cat(sprintf("  Phase 1 with IDV %s:\n", f))
  cat(sprintf("    robust UCL                     : %.5f\n", c_$ucl_r))
  cat(sprintf("    healthy/contaminated weight ratio: %.3f\n", c_$ratio_w))
  cat(sprintf("    contaminated among the 6 lowest: %d of 6\n", c_$n_abajo))
  cat(sprintf("    shift of the contaminated ones : %.2f SD\n", c_$despl))
  cat(sprintf("    DRIFT OF THE CENTER mu_r       : %.2f SD\n", c_$corrimiento_centro))
  cat(sprintf("    determinant inflation          : %.4f\n\n", c_$infl))
}

cat("=== CROSS MATRIX ===\n")
filas <- list()
for (f1 in c("7", "1")) for (f2 in c("7", "1")) {
  r <- monitorear(cals[[f1]], as.numeric(f2))
  cat(sprintf("  Phase1 IDV %s -> Phase2 IDV %s : robust %2d of 20, %d false alarms of 10 | classical %2d of 20, %d FA\n",
              f1, f2, r$det_rob, r$fa_rob, r$det_cls, r$fa_cls))
  filas[[length(filas) + 1]] <- data.frame(
    phase1 = paste0("IDV ", f1), phase2 = paste0("IDV ", f2),
    UCL_rob = cals[[f1]]$ucl_r,
    weight_ratio = cals[[f1]]$ratio_w,
    center_drift_SD = cals[[f1]]$corrimiento_centro,
    det_rob = r$det_rob, fa_rob = r$fa_rob,
    det_cls = r$det_cls, fa_cls = r$fa_cls,
    T2_median_healthy = r$t2_sano_med,
    T2_median_faulty  = r$t2_fallo_med
  )
}
res <- do.call(rbind, filas)

# =====================================================================
#  VERDICT
# =====================================================================
cat("\n=========================================================\n")
cat(" VERDICT: where do the false alarms come from?\n")
cat("=========================================================\n")

fa <- function(f1, f2) res$fa_rob[res$phase1 == paste0("IDV ", f1) &
                                 res$phase2 == paste0("IDV ", f2)]
cat(sprintf("  Phase1 IDV 7: %d false alarms with Phase2 IDV 7, %d with Phase2 IDV 1\n",
            fa(7, 7), fa(7, 1)))
cat(sprintf("  Phase1 IDV 1: %d false alarms with Phase2 IDV 7, %d with Phase2 IDV 1\n\n",
            fa(1, 7), fa(1, 1)))

diag_limpia  <- fa(7, 7) == 0 && fa(1, 1) == 0
fila1_falla  <- fa(1, 7) > 0 && fa(1, 1) > 0

if (fila1_falla) {
  cat("  EXPLANATION (b): THE CONTAMINATION ITSELF.\n")
  cat("  The IDV 1 row produces false alarms with BOTH Phase 2 sets, so\n")
  cat("  it is not a problem of crossing different faults. The center\n")
  cat("  mu_r is a SIMPLE average, with no weights, and it carries the\n")
  cat("  shift of the six contaminated batches.\n")
  cat("  This is the limitation the manuscript already states: the\n")
  cat("  protection of the average depends on the shift being moderate.\n")
} else if (diag_limpia) {
  cat("  EXPLANATION (a): THE CROSS.\n")
  cat("  The diagonal is clean and only the crossed cells fail.\n")
  cat("  Calibrating with one fault and watching another misplaces the\n")
  cat("  reference. The coherent scenario (same fault in both phases)\n")
  cat("  works.\n")
} else {
  cat("  MIXED: review cell by cell before concluding.\n")
}

cat("\n  The column that decides it is the drift of the center:\n")
for (f in c("7", "1")) {
  c_ <- cals[[f]]
  cat(sprintf("    IDV %s: the contaminated batches shift %.2f SD and the center\n",
              f, c_$despl))
  cat(sprintf("            drifts %.2f SD (about one fifth, because they are\n",
              c_$corrimiento_centro))
  cat("            6 batches out of 30 in a simple average)\n")
}

cat("\n=== FULL TABLE ===\n")
print(res, row.names = FALSE)

if (!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)
ruta <- file.path(DIR_OUT, "phase1_cross_matrix.csv")
write.csv(res, ruta, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", ruta))
cat("\n=== END. Copy the whole output. ===\n\n")
