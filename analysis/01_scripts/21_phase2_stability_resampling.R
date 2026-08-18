# =====================================================================
# 21_phase2_stability_resampling.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   How much the Phase 2 result depends on which runs make up Phase 2.
#   Phase 1 is held frozen and the monitoring is repeated over twenty
#   different Phase 2 compositions. For each method it reports the mean,
#   standard deviation, standard error, minimum and maximum of the
#   faulty batches detected out of twenty and of the false alarms out of
#   ten, plus the margin with which the proposed method separates the
#   faulty batches from its limit and how many faulty batches the
#   classical method leaves within 20% of its own limit.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - Phase 1 is FROZEN: 24 healthy runs and 6 runs with IDV 7, exactly
#     as in PHASE1_PIPELINE_STEP30.R. Both UCLs are therefore constants
#     and the only thing varying between compositions is the Phase 2
#     runs.
#   - Phase 1 runs are never reused in Phase 2: the healthy draw starts
#     at run 25 and the faulty draw at run 7.
#   - Composition 1 is the published one, so it serves as an anchor;
#     compositions 2 to 20 are drawn with a fixed seed and are therefore
#     reproducible.
#   - Batches are built by spaced subsampling, one observation every
#     PASO = 30, so the observations of a batch can be treated as
#     independent in time.
#
# ANCHORS AND THEIR EXPECTED VALUES
#   A1  robust UCL    = 19.69285  (tolerance 1e-04)
#       classical UCL = 19.46440  (tolerance 1e-04)
#   A2  Composition 1 must reproduce 20 of 20 detected with the robust
#       method, 2 of 20 with the classical one, and 0 false alarms in
#       both.
#   If A1 or A2 fails the script STOPS: Phase 1 would not be the one of
#   the manuscript and nothing downstream would be comparable.
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   The proposed method is declared STABLE if it detects 20 of 20 in
#   all twenty compositions (SD = 0) with no false alarms. If it fails
#   in any composition the exact figure is reported instead; the design
#   is not adjusted after seeing the result.
#
# OUTPUT
#   02_resultados/phase2_stability_by_composition.csv
#   02_resultados/phase2_stability_summary.csv
# =====================================================================

library(robustT2AFM)

# ------------------------- Configuration -----------------------------
VARS        <- c("xmv_9", "xmv_8", "xmeas_19", "xmeas_17")
I_LOTE      <- 20
PASO        <- 30
J           <- 4
ALPHA       <- 0.001
H_MCD       <- 0.67
FALLO       <- 7
DESDE_SANO  <- 1
DESDE_FALLO <- 161
N_SANOS_F1  <- 24
N_CONTA_F1  <- 6
N_SANOS_F2  <- 10
N_FALLO_F2  <- 20
N_COMPOS    <- 20
SEMILLA     <- 2026
DIR_DATOS   <- "04_datos"
DIR_OUT     <- "02_resultados"

muestras_esp <- function(d) seq(d, by = PASO, length.out = I_LOTE)

# ----------------------------- Data ----------------------------------
if (!exists("ff_test")) {
  cat("Loading FaultFree_Testing...\n")
  ff_test <- read.csv(file.path(DIR_DATOS, "TEP_FaultFree_Testing.csv"))
}
if (!exists("fy")) {
  cat("Loading Faulty_Testing...\n")
  fy <- read.csv(file.path(DIR_DATOS, "TEP_Faulty_Testing.csv"))
}

extraer_lote <- function(datos, run_id, fault_id, desde, batch_name) {
  sel <- datos$simulationRun == run_id & datos$faultNumber == fault_id &
    datos$sample %in% muestras_esp(desde)
  sub <- datos[sel, VARS]
  if (nrow(sub) != I_LOTE) return(NULL)
  data.frame(Batch = batch_name, sub, row.names = NULL)
}

# ------------------- Pools of runs available for Phase 2 --------------
runs_sanos_disp <- sort(unique(ff_test$simulationRun))
runs_fallo_disp <- sort(unique(fy$simulationRun[fy$faultNumber == FALLO]))

pool_sanos <- setdiff(runs_sanos_disp, 1:N_SANOS_F1)   # excludes Phase 1
pool_fallo <- setdiff(runs_fallo_disp, 1:N_CONTA_F1)   # excludes Phase 1

cat("\n=== Pools available for Phase 2 (Phase 1 runs excluded) ===\n")
cat(sprintf("  healthy available : %d  (from run %d to %d)\n",
            length(pool_sanos), min(pool_sanos), max(pool_sanos)))
cat(sprintf("  faulty available  : %d  (from run %d to %d)\n",
            length(pool_fallo), min(pool_fallo), max(pool_fallo)))
if (length(pool_sanos) < N_SANOS_F2 || length(pool_fallo) < N_FALLO_F2) {
  stop("Not enough runs left to compose Phase 2.")
}

# ===================== FROZEN PHASE 1 + ANCHOR A1 =====================
f1 <- do.call(rbind, c(
  lapply(1:N_SANOS_F1, function(r)
    extraer_lote(ff_test, r, 0, DESDE_SANO, sprintf("F1_sano_%02d", r))),
  lapply(1:N_CONTA_F1, function(r)
    extraer_lote(fy, r, FALLO, DESDE_FALLO, sprintf("F1_fallo_%02d", r)))
))
stopifnot(length(unique(f1$Batch)) == N_SANOS_F1 + N_CONTA_F1)

cal_rob <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
cal_cls <- hotelling_classical_calibrate(f1, VARS)
K       <- length(unique(f1$Batch))
ucl_rob <- ucl_F_adjusted(cal_rob, I = I_LOTE, alpha = ALPHA)$UCL
ucl_cls <- hotelling_classical_ucl(K = K, I = I_LOTE, J = J,
                                   alpha = ALPHA, phase = "II")$UCL

cat("\n=== Anchor A1: the two UCLs of the frozen Phase 1 ===\n")
cat(sprintf("  robust UCL    = %.5f  (expected 19.69285)  %s\n",
            ucl_rob, if (abs(ucl_rob - 19.69285) < 1e-4) "OK" else "*** FAILS ***"))
cat(sprintf("  classical UCL = %.5f  (expected 19.46440)  %s\n",
            ucl_cls, if (abs(ucl_cls - 19.46440) < 1e-4) "OK" else "*** FAILS ***"))
if (abs(ucl_rob - 19.69285) > 1e-4 || abs(ucl_cls - 19.46440) > 1e-4) {
  stop("Anchor A1 failed. Phase 1 is not the one of the manuscript. STOP.")
}

# ------------------ Evaluate one Phase 2 composition ------------------
evaluar <- function(runs_sanos, runs_fallo) {
  f2 <- do.call(rbind, c(
    lapply(runs_sanos, function(r)
      extraer_lote(ff_test, r, 0, DESDE_SANO, sprintf("F2_sano_%03d", r))),
    lapply(runs_fallo, function(r)
      extraer_lote(fy, r, FALLO, DESDE_FALLO, sprintf("F2_fallo_%03d", r)))
  ))
  if (length(unique(f2$Batch)) != (N_SANOS_F2 + N_FALLO_F2)) return(NULL)

  b  <- unique(f2$Batch)
  ef <- grepl("fallo", b)
  mr <- monitor_afm_mcd(f2, cal_rob, VARS)
  mc <- hotelling_classical_monitor(f2, cal_cls, VARS)
  tr <- mr$T2[match(b, mr$Batch)]
  tc <- mc$T2[match(b, mc$Batch)]

  # how many faulty batches the classical method leaves within +-20% of its limit
  rozando <- sum(tc[ef] >= 0.8 * ucl_cls & tc[ef] <= 1.2 * ucl_cls)

  list(det_rob = sum(tr[ef]  > ucl_rob),
       det_cls = sum(tc[ef]  > ucl_cls),
       fa_rob  = sum(tr[!ef] > ucl_rob),
       fa_cls  = sum(tc[!ef] > ucl_cls),
       t2_rob_fallo_min = min(tr[ef]),
       t2_cls_fallo_med = median(tc[ef]),
       cls_rozando_ucl  = rozando)
}

# ================= Composition 1: the published one (ANCHOR A2) =======
comp <- list(list(sanos = 25:34, fallo = 7:26, etiqueta = "publicada"))

set.seed(SEMILLA)
for (i in 2:N_COMPOS) {
  comp[[i]] <- list(
    sanos    = sort(sample(pool_sanos, N_SANOS_F2)),
    fallo    = sort(sample(pool_fallo, N_FALLO_F2)),
    etiqueta = sprintf("sorteo_%02d", i - 1)
  )
}

# ============================ Execution ===============================
res <- data.frame()
cat("\n=== Evaluating", N_COMPOS, "Phase 2 compositions ===\n")
for (i in seq_len(N_COMPOS)) {
  r <- evaluar(comp[[i]]$sanos, comp[[i]]$fallo)
  if (is.null(r)) {
    cat(sprintf("  %2d  %-12s  SKIPPED (incomplete runs)\n",
                i, comp[[i]]$etiqueta))
    next
  }
  cat(sprintf("  %2d  %-12s  robust %2d/20 (FA %d/10)  |  classical %2d/20 (FA %d/10)\n",
              i, comp[[i]]$etiqueta, r$det_rob, r$fa_rob, r$det_cls, r$fa_cls))
  res <- rbind(res, data.frame(
    composition  = i,
    label        = comp[[i]]$etiqueta,
    healthy_runs = paste(comp[[i]]$sanos, collapse = "|"),
    faulty_runs  = paste(comp[[i]]$fallo, collapse = "|"),
    det_rob = r$det_rob, det_cls = r$det_cls,
    fa_rob  = r$fa_rob,  fa_cls  = r$fa_cls,
    t2_rob_faulty_min    = r$t2_rob_fallo_min,
    t2_cls_faulty_median = r$t2_cls_fallo_med,
    cls_near_ucl         = r$cls_rozando_ucl,
    stringsAsFactors = FALSE
  ))
}

# ========================== ANCHOR A2 =================================
cat("\n=== Anchor A2: the published composition ===\n")
a2 <- res[res$label == "publicada", ]
ok2 <- nrow(a2) == 1 && a2$det_rob == 20 && a2$det_cls == 2 &&
  a2$fa_rob == 0 && a2$fa_cls == 0
cat(sprintf("  robust %d/20 (exp 20) | classical %d/20 (exp 2) | FA %d and %d (exp 0 and 0)  %s\n",
            a2$det_rob, a2$det_cls, a2$fa_rob, a2$fa_cls,
            if (ok2) "OK" else "*** FAILS ***"))
if (!ok2) stop("Anchor A2 failed. STOP: it does not reproduce the published table.")

# ============================ SUMMARY =================================
n <- nrow(res)
resumir <- function(x) c(media = mean(x), sd = sd(x), se = sd(x) / sqrt(length(x)),
                         min = min(x), max = max(x))

rr <- resumir(res$det_rob)
rc <- resumir(res$det_cls)

cat("\n=========== STABILITY OVER", n, "COMPOSITIONS ===========\n")
cat(sprintf("  Faulty batches detected, out of 20:\n"))
cat(sprintf("    Proposed   mean %.2f | SD %.2f | SE %.2f | range [%d, %d]\n",
            rr["media"], rr["sd"], rr["se"], rr["min"], rr["max"]))
cat(sprintf("    Classical  mean %.2f | SD %.2f | SE %.2f | range [%d, %d]\n",
            rc["media"], rc["sd"], rc["se"], rc["min"], rc["max"]))
cat(sprintf("\n  False alarms, out of 10:\n"))
cat(sprintf("    Proposed   mean %.2f | range [%d, %d]\n",
            mean(res$fa_rob), min(res$fa_rob), max(res$fa_rob)))
cat(sprintf("    Classical  mean %.2f | range [%d, %d]\n",
            mean(res$fa_cls), min(res$fa_cls), max(res$fa_cls)))
cat(sprintf("\n  Margin of the proposed method: the lowest T2 of a faulty batch,\n"))
cat(sprintf("  over the %d compositions, is %.1f, against a UCL of %.2f.\n",
            n, min(res$t2_rob_faulty_min), ucl_rob))
cat(sprintf("  Faulty batches the classical method leaves within 20%% of its\n"))
cat(sprintf("  limit: mean %.1f of 20 per composition.\n",
            mean(res$cls_near_ucl)))

# ============================ VERDICT =================================
estable_rob <- rr["sd"] == 0 && rr["min"] == 20 && max(res$fa_rob) == 0
cat("\n=========================================================\n")
if (estable_rob) {
  cat(" VERDICT: the proposed method is STABLE. It detects 20 of 20 in\n")
  cat(" the", n, "compositions, with SD = 0 and no false alarms.\n")
  cat(" The sentence in the manuscript can be replaced by figures.\n")
} else {
  cat(" VERDICT: the proposed method is NOT stable over the", n, "\n")
  cat(" compositions. Report the real range and adjust the sentence in\n")
  cat(" the manuscript. It remains far more stable than the classical\n")
  cat(" one, but the current claim would have to be qualified.\n")
}
cat(sprintf(" Classical: SD = %.2f over a range of %d to %d out of 20.\n",
            rc["sd"], rc["min"], rc["max"]))
cat("=========================================================\n")

# ============================ SAVING ==================================
if (!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)

write.csv(res, file.path(DIR_OUT, "phase2_stability_by_composition.csv"),
          row.names = FALSE)

resumen <- data.frame(
  method          = c("Propuesto (AFM-MCD)", "Clasico (Hotelling)"),
  n_compositions  = n,
  det_mean  = c(rr["media"], rc["media"]),
  det_sd    = c(rr["sd"],    rc["sd"]),
  det_se    = c(rr["se"],    rc["se"]),
  det_min   = c(rr["min"],   rc["min"]),
  det_max   = c(rr["max"],   rc["max"]),
  fa_mean   = c(mean(res$fa_rob), mean(res$fa_cls)),
  fa_min    = c(min(res$fa_rob),  min(res$fa_cls)),
  fa_max    = c(max(res$fa_rob),  max(res$fa_cls)),
  UCL       = c(ucl_rob, ucl_cls),
  row.names = NULL
)
write.csv(resumen, file.path(DIR_OUT, "phase2_stability_summary.csv"),
          row.names = FALSE)

cat("\nSaved:\n")
cat("  phase2_stability_by_composition.csv\n")
cat("  phase2_stability_summary.csv\n")
cat("\n=== END. Copy the whole output. ===\n\n")
