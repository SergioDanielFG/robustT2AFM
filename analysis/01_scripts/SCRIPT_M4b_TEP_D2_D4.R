# =====================================================================
# SCRIPT_M4b_TEP_D2_D4.R
# ---------------------------------------------------------------------
# Revision round 1, module M4, second part (Reviewer 1, comments 4, 6, 7
# and 11; items D2, D3 and D4 of REGISTRO_CAMBIOS_v9.md).
#
# WHY
#   Three figures of Section 4 were computed with the published method
#   (v7: no scaling, mean of the batch centers, m* = 13) and have no
#   counterpart with the final method in 02_resultados:
#     D2  covariance inflation det(S_classical) / det(S_w), Phase 1 IDV 7
#         (v7: 16.14; Section 4.4)
#     D3  stability of Phase 2 over the 20 compositions of
#         21_phase2_stability_resampling.R (v7: 19.95 of 20, SD 0.22,
#         range 19-20; classical 1.40, SD 1.31, range 0-5; lowest T2 of
#         a faulty batch 26.6 except in one composition; Section 4.5)
#     D4  weight ratio healthy / contaminated, Phase 1 IDV 1
#         (v7: 7.1; Section 4.7)
#   Decision of 02-oct-2026: recompute the three with the final method
#   (robustT2AFM 0.3.0 defaults: MAD scaling for the weights, MCD of the
#   batch centers, Eq. (8) with m* = 14) and report what comes out.
#
# WHAT IT COMPUTES
#   NOTHING IS SIMULATED. The 20 compositions are read from
#   02_resultados/phase2_stability_by_composition.csv (the published
#   ones), so Phase 2 is exactly the one of the manuscript. The
#   simulated limit of NEW is read from table_M4_final_A_limits.csv.
#
# DESIGN
#   - Batches, runs, starting samples and step as in
#     SCRIPT_M4_TEP_final_mstar14.R (copied without changes).
#   - Every calibration is preceded by set.seed(2026), as in M4 final,
#     so the per-batch MCD fits are those of table_M4_final_*.
#   - V7  = calibrate_afm_mcd(scaling = "none", center = "mean"),
#           Eq. (8) with m* = 13 (m_star = "nominal").
#   - NEW = calibrate_afm_mcd() with the 0.3.0 defaults, Eq. (8) with
#           m* = 14; in D3 also counted with its simulated limit.
#   - Values are saved without rounding (the convention of the paper
#     is applied when writing the text).
#
# ANCHORS (STOP if any fails)
#   A1  V7, Phase 1 IDV 7: UCL 19.69285; weight ratio 6.44;
#       inflation 16.14
#   A2  NEW, Phase 1 IDV 7: UCL 19.64480; weight ratio 4.02; 6 of 6
#       contaminated batches lowest; published Phase 2: 20/20, 0/10,
#       lowest faulty T2 41.9
#   A3  V7, Phase 1 IDV 1: weight ratio 7.1 (one decimal)
#   A4  V7 and classical over the 20 compositions: detected and false
#       alarms equal to phase2_stability_by_composition.csv in every
#       composition (means 19.95 and 1.40)
#   A5  NEW, Phase 1 IDV 1, published Phase 2: 20/20, 0/10
#       (table_M4_final_E_center.csv)
# CHECKS (printed, do not stop)
#   C1  V7 lowest faulty T2 per composition vs the CSV (max difference)
#   C2  the 20 compositions regenerated with set.seed(2026), as in
#       script 21, coincide with the CSV
#
# OUTPUT (05_revision_R1/02_resultados/)
#   table_M4b_D2_D4.csv            inflation and weight ratios
#   table_M4b_D3_compositions.csv  one row per composition and chart
#   table_M4b_D3_summary.csv       summary per chart and limit
#   Run from C:/temp_paper. Expected time: under one minute.
# =====================================================================
library(robustT2AFM)
stopifnot(packageVersion("robustT2AFM") >= "0.3.0")
DIR_DATOS <- "04_datos"
DIR_OUT   <- "05_revision_R1/02_resultados"
VARS   <- c("xmv_9", "xmv_8", "xmeas_19", "xmeas_17")
I_LOTE <- 20; PASO <- 30; J <- 4
ALPHA  <- 0.001; H_MCD <- 0.67
DESDE_SANO <- 1; DESDE_FALLO <- 161
SEED_CAL <- 2026

if (!exists("ff_test")) { cat("Loading FaultFree_Testing...\n"); ff_test <- read.csv(file.path(DIR_DATOS, "TEP_FaultFree_Testing.csv")) }
if (!exists("fy"))      { cat("Loading Faulty_Testing...\n");    fy      <- read.csv(file.path(DIR_DATOS, "TEP_Faulty_Testing.csv")) }

comp_v7 <- read.csv("02_resultados/phase2_stability_by_composition.csv", stringsAsFactors = FALSE)
lim_A   <- read.csv(file.path(DIR_OUT, "table_M4_final_A_limits.csv"), stringsAsFactors = FALSE)
UCL_NEW_SIM <- lim_A$UCL[lim_A$phase1 == "IDV 7, step 30" & lim_A$chart == "NEW" & lim_A$limit == "simulated"]
stopifnot(nrow(comp_v7) == 20, length(UCL_NEW_SIM) == 1)

# ---------------------------------------------------------------------
# Batches (as in M4 final)
# ---------------------------------------------------------------------
lote <- function(dat, run, fault, desde, nom) {
  sel <- dat$simulationRun == run & dat$faultNumber == fault &
         dat$sample %in% seq(desde, by = PASO, length.out = I_LOTE)
  s <- dat[sel, VARS]
  if (nrow(s) != I_LOTE) return(NULL)
  data.frame(Batch = nom, s, row.names = NULL)
}
fase1 <- function(fault_f1) do.call(rbind, c(
  lapply(1:24, function(r) lote(ff_test, r, 0,        DESDE_SANO,  sprintf("F1_sano_%02d",  r))),
  lapply(1:6,  function(r) lote(fy,      r, fault_f1, DESDE_FALLO, sprintf("F1_fallo_%02d", r)))))
fase2 <- function(sanos, fallo, fault_f2 = 7) do.call(rbind, c(
  lapply(sanos, function(r) lote(ff_test, r, 0,        DESDE_SANO,  sprintf("F2_sano_%03d",  r))),
  lapply(fallo, function(r) lote(fy,      r, fault_f2, DESDE_FALLO, sprintf("F2_fallo_%03d", r)))))
runs <- function(s) as.integer(strsplit(s, "|", fixed = TRUE)[[1]])

# ---------------------------------------------------------------------
# Calibration of V7, NEW and classical on one Phase 1
# ---------------------------------------------------------------------
calibrar <- function(f1) {
  set.seed(SEED_CAL)
  v7  <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD, scaling = "none", center = "mean"))
  set.seed(SEED_CAL)
  new <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
  cls <- hotelling_classical_calibrate(f1, VARS)
  list(v7 = v7, new = new, cls = cls,
       ucl_v7  = ucl_F_adjusted(v7,  I = I_LOTE, alpha = ALPHA, m_star = "nominal")$UCL,
       ucl_new = ucl_F_adjusted(new, I = I_LOTE, alpha = ALPHA)$UCL,
       ucl_cls = hotelling_classical_ucl(K = cls$n_batches, I = I_LOTE, J = J, alpha = ALPHA, phase = "II")$UCL)
}
ratio_w <- function(w) { ic <- grepl("fallo", names(w)); mean(w[!ic]) / mean(w[ic]) }
low6    <- function(w) { ic <- grepl("fallo", names(w)); sum(ic[order(w)][1:6]) }
t2_of <- function(p2, K, chart) {
  b <- unique(p2$Batch)
  m <- switch(chart,
              V7        = monitor_afm_mcd(p2, K$v7,  VARS),
              NEW       = monitor_afm_mcd(p2, K$new, VARS),
              classical = hotelling_classical_monitor(p2, K$cls, VARS))
  setNames(m$T2[match(b, m$Batch)], b)
}
cuenta <- function(t2, ucl) {
  ef <- grepl("fallo", names(t2))
  c(det = sum(t2[ef] > ucl), fa = sum(t2[!ef] > ucl),
    t2_min_faulty = min(t2[ef]),
    near_limit = sum(t2[ef] >= 0.8 * ucl & t2[ef] <= 1.2 * ucl))
}
same <- function(x, v, d) abs(round(x, d) - v) < 1e-9
p2_pub <- fase2(25:34, 7:26)

# =====================================================================
# Phase 1 with IDV 7: D2 and anchors A1, A2
# =====================================================================
cat("== Phase 1 with IDV 7 ==\n")
K7 <- calibrar(fase1(7))
infl7_v7  <- det(K7$cls$Sp) / det(K7$v7$Sw)
infl7_new <- det(K7$cls$Sp) / det(K7$new$Sw)
r7_v7  <- ratio_w(K7$v7$weights);  r7_new <- ratio_w(K7$new$weights)
c7_new <- cuenta(t2_of(p2_pub, K7, "NEW"), K7$ucl_new)

okA1 <- abs(K7$ucl_v7 - 19.69285) < 1e-4 && same(r7_v7, 6.44, 2) && same(infl7_v7, 16.14, 2)
cat(sprintf("A1 V7 : UCL %.5f | ratio %.4f | inflation %.4f -> %s\n",
            K7$ucl_v7, r7_v7, infl7_v7, if (okA1) "OK" else "*** FAIL ***"))
if (!okA1) stop("Anchor A1 failed: V7 does not reproduce the manuscript. STOP.")

okA2 <- abs(K7$ucl_new - 19.64480) < 1e-4 && same(r7_new, 4.02, 2) && low6(K7$new$weights) == 6 &&
        c7_new["det"] == 20 && c7_new["fa"] == 0 && same(c7_new["t2_min_faulty"], 41.9, 1)
cat(sprintf("A2 NEW: UCL %.5f | ratio %.4f | %d of 6 lowest | %d/20, FA %d | min T2 %.4f -> %s\n",
            K7$ucl_new, r7_new, low6(K7$new$weights), c7_new["det"], c7_new["fa"],
            c7_new["t2_min_faulty"], if (okA2) "OK" else "*** FAIL ***"))
if (!okA2) stop("Anchor A2 failed: NEW does not reproduce table_M4_final_B. STOP.")

# =====================================================================
# Phase 1 with IDV 1: D4 and anchors A3, A5
# =====================================================================
cat("\n== Phase 1 with IDV 1 ==\n")
K1 <- calibrar(fase1(1))
infl1_v7  <- det(K1$cls$Sp) / det(K1$v7$Sw)
infl1_new <- det(K1$cls$Sp) / det(K1$new$Sw)
r1_v7  <- ratio_w(K1$v7$weights);  r1_new <- ratio_w(K1$new$weights)
c1_new <- cuenta(t2_of(p2_pub, K1, "NEW"), K1$ucl_new)

okA3 <- same(r1_v7, 7.1, 1)
cat(sprintf("A3 V7 : ratio %.4f (expected 7.1) -> %s\n", r1_v7, if (okA3) "OK" else "*** FAIL ***"))
if (!okA3) stop("Anchor A3 failed: V7 with IDV 1 does not reproduce 7.1. STOP.")
okA5 <- c1_new["det"] == 20 && c1_new["fa"] == 0
cat(sprintf("A5 NEW: published Phase 2 %d/20, FA %d -> %s\n",
            c1_new["det"], c1_new["fa"], if (okA5) "OK" else "*** FAIL ***"))
if (!okA5) stop("Anchor A5 failed: NEW with IDV 1 does not reproduce table_M4_final_E. STOP.")

# =====================================================================
# D3: the 20 published compositions, Phase 1 IDV 7 frozen
# =====================================================================
cat("\n== D3: 20 Phase 2 compositions ==\n")
LIMS <- data.frame(chart = c("V7", "NEW", "NEW", "classical"),
                   limit = c("Eq8_mstar13_published", "Eq8_mstar14", "simulated", "F"),
                   UCL   = c(K7$ucl_v7, K7$ucl_new, UCL_NEW_SIM, K7$ucl_cls),
                   stringsAsFactors = FALSE)
filas <- list()
for (i in seq_len(nrow(comp_v7))) {
  p2 <- fase2(runs(comp_v7$healthy_runs[i]), runs(comp_v7$faulty_runs[i]))
  stopifnot(length(unique(p2$Batch)) == 30)
  t2 <- list(V7 = t2_of(p2, K7, "V7"), NEW = t2_of(p2, K7, "NEW"), classical = t2_of(p2, K7, "classical"))
  for (j in seq_len(nrow(LIMS))) {
    cc <- cuenta(t2[[LIMS$chart[j]]], LIMS$UCL[j])
    filas[[length(filas) + 1]] <- data.frame(
      composition = comp_v7$composition[i], label = comp_v7$label[i],
      chart = LIMS$chart[j], limit = LIMS$limit[j], UCL = LIMS$UCL[j], t(cc), row.names = NULL)
  }
}
tabD3 <- do.call(rbind, filas)

v7d <- tabD3[tabD3$chart == "V7", ]; cld <- tabD3[tabD3$chart == "classical", ]
okA4 <- all(v7d$det == comp_v7$det_rob) && all(v7d$fa == comp_v7$fa_rob) &&
        all(cld$det == comp_v7$det_cls) && all(cld$fa == comp_v7$fa_cls)
cat(sprintf("A4 V7 and classical, 20 compositions = CSV: %s (means %.2f and %.2f)\n",
            if (okA4) "OK" else "*** FAIL ***", mean(v7d$det), mean(cld$det)))
if (!okA4) stop("Anchor A4 failed: V7 or classical do not reproduce the 20 compositions. STOP.")
cat(sprintf("C1 V7 lowest faulty T2 vs CSV, max difference: %.2e\n",
            max(abs(v7d$t2_min_faulty - comp_v7$t2_rob_faulty_min))))

pool_s <- setdiff(sort(unique(ff_test$simulationRun)), 1:24)
pool_f <- setdiff(sort(unique(fy$simulationRun[fy$faultNumber == 7])), 1:6)
set.seed(2026); okC2 <- TRUE
for (i in 2:20) {
  s <- sort(sample(pool_s, 10)); f <- sort(sample(pool_f, 20))
  okC2 <- okC2 && all(s == runs(comp_v7$healthy_runs[i])) && all(f == runs(comp_v7$faulty_runs[i]))
}
cat(sprintf("C2 compositions regenerated with set.seed(2026) = CSV: %s\n", if (okC2) "OK" else "REVIEW"))

resumir <- function(d) data.frame(
  chart = d$chart[1], limit = d$limit[1], UCL = d$UCL[1], n_compositions = nrow(d),
  det_mean = mean(d$det), det_sd = sd(d$det), det_min = min(d$det), det_max = max(d$det),
  faulty_detected_total = sum(d$det), n_comp_20of20 = sum(d$det == 20), n_comp_zero = sum(d$det == 0),
  fa_total = sum(d$fa), fa_max = max(d$fa),
  t2min_lowest = min(d$t2_min_faulty), t2min_second = sort(d$t2_min_faulty)[2],
  near_limit_mean = mean(d$near_limit))
tabS <- do.call(rbind, lapply(seq_len(nrow(LIMS)), function(j)
  resumir(tabD3[tabD3$chart == LIMS$chart[j] & tabD3$limit == LIMS$limit[j], ])))
cat("\n== D3 summary ==\n"); print(tabS, row.names = FALSE, digits = 6)

# =====================================================================
# D2 and D4 table, output
# =====================================================================
tabDD <- data.frame(
  item     = c("D2", "D2", "D2", "D2", "D4", "D4", "D4", "D4"),
  phase1   = c("IDV 7", "IDV 7", "IDV 1", "IDV 1", "IDV 7", "IDV 7", "IDV 1", "IDV 1"),
  chart    = rep(c("V7", "NEW"), 4),
  quantity = c(rep("inflation det(S_classical)/det(S_w)", 4), rep("weight ratio healthy/contaminated", 4)),
  value    = c(infl7_v7, infl7_new, infl1_v7, infl1_new, r7_v7, r7_new, r1_v7, r1_new),
  contaminated_in_6_lowest = c(NA, NA, NA, NA,
                               low6(K7$v7$weights), low6(K7$new$weights),
                               low6(K1$v7$weights), low6(K1$new$weights)))
cat("\n== D2 and D4 ==\n"); print(tabDD, row.names = FALSE, digits = 8)

write.csv(tabDD, file.path(DIR_OUT, "table_M4b_D2_D4.csv"),           row.names = FALSE)
write.csv(tabD3, file.path(DIR_OUT, "table_M4b_D3_compositions.csv"), row.names = FALSE)
write.csv(tabS,  file.path(DIR_OUT, "table_M4b_D3_summary.csv"),      row.names = FALSE)
cat("\nAnchors A1-A5 OK | C2", if (okC2) "OK" else "REVIEW",
    "\nSaved: table_M4b_D2_D4.csv, table_M4b_D3_compositions.csv, table_M4b_D3_summary.csv\n")
