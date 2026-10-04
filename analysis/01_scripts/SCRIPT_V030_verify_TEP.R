# =====================================================================
# SCRIPT_V030_verify_TEP.R
# ---------------------------------------------------------------------
# Checks that the installed package robustT2AFM 0.3.0 reproduces, on the
# Tennessee Eastman data, what SCRIPT_M4_TEP_final_method.R computed with
# the method programmed inside the script. Nothing new is estimated.
#
# CHECKS
#   V1  Legacy options (scaling = "none", center = "mean") reproduce the
#       published application exactly: UCL 19.69285, weight ratio 6.44,
#       IDV 7 detected 20/20 with 0/10 false alarms.
#   V2  Default options (MAD + MCD of the centers), Phase 1 with IDV 7:
#       weight ratio 4.02 and 6 of 6 contaminated batches lowest (same as
#       M4, since the per-batch MCD fits are the same); IDV 7 20/20, 0/10.
#   V3  Default options, Phase 1 with IDV 1 (whole shifted batches):
#       IDV 7 detected 20/20 with 0/10 false alarms (M4: 20/20, 0/10).
#   V4  The minimum T2 of the faulty batches is close to M4 (41.9 and
#       32.4). It is not identical: M4 seeded the MCD of the centers with
#       2027, the package uses its internal seed.
#   V5  Invariance on real data: xmeas_17 multiplied by 100 leaves the
#       weights and the Phase 2 T2 unchanged.
#   V6  The same Phase 1 calibrated twice gives the same mu_r.
#
# OUTPUT: printed only. Run from C:/temp_paper.
# =====================================================================
library(robustT2AFM)
stopifnot(packageVersion("robustT2AFM") >= "0.3.0")
DIR_DATOS <- "04_datos"
VARS <- c("xmv_9", "xmv_8", "xmeas_19", "xmeas_17")
I_LOTE <- 20; PASO <- 30; ALPHA <- 0.001
DESDE_SANO <- 1; DESDE_FALLO <- 161; SEED_CAL <- 2026

if (!exists("ff_test")) ff_test <- read.csv(file.path(DIR_DATOS, "TEP_FaultFree_Testing.csv"))
if (!exists("fy"))      fy      <- read.csv(file.path(DIR_DATOS, "TEP_Faulty_Testing.csv"))

lote <- function(dat, run, fault, desde, nom) {
  sel <- dat$simulationRun == run & dat$faultNumber == fault &
         dat$sample %in% seq(desde, by = PASO, length.out = I_LOTE)
  s <- dat[sel, VARS]; if (nrow(s) != I_LOTE) return(NULL)
  data.frame(Batch = nom, s, row.names = NULL)
}
fase1 <- function(f) do.call(rbind, c(
  lapply(1:24, function(r) lote(ff_test, r, 0, DESDE_SANO,  sprintf("F1_sano_%02d",  r))),
  lapply(1:6,  function(r) lote(fy,      r, f, DESDE_FALLO, sprintf("F1_fallo_%02d", r)))))
fase2 <- function(f) do.call(rbind, c(
  lapply(25:34, function(r) lote(ff_test, r, 0, DESDE_SANO,  sprintf("F2_sano_%02d",  r))),
  lapply(7:26,  function(r) lote(fy,      r, f, DESDE_FALLO, sprintf("F2_fallo_%03d", r)))))

evalua <- function(f1, p2, ...) {
  set.seed(SEED_CAL)
  cal <- calibrate_afm_mcd(f1, VARS, ...)
  ucl <- ucl_F_adjusted(cal, I = I_LOTE, alpha = ALPHA)$UCL
  mon <- monitor_afm_mcd(p2, cal, VARS)
  t2  <- setNames(mon$T2, mon$Batch); ef <- grepl("fallo", names(t2))
  w   <- cal$weights; ic <- grepl("fallo", names(w))
  list(cal = cal, ucl = ucl, t2 = t2,
       det = sum(t2[ef] > ucl), fa = sum(t2[!ef] > ucl), t2min = min(t2[ef]),
       ratio = mean(w[!ic]) / mean(w[ic]), low6 = sum(ic[order(w)][1:6]))
}
ok <- function(x) if (isTRUE(x)) "OK" else "*** FAIL ***"

f1_7 <- fase1(7); f1_1 <- fase1(1); p2_7 <- fase2(7)

leg <- evalua(f1_7, p2_7, scaling = "none", center = "mean")
cat(sprintf("V1 legacy   : UCL %.5f | ratio %.2f | %d/20 | FA %d  -> %s\n",
            leg$ucl, leg$ratio, leg$det, leg$fa,
            ok(abs(leg$ucl - 19.69285) < 1e-4 && round(leg$ratio, 2) == 6.44 && leg$det == 20 && leg$fa == 0)))

n7 <- evalua(f1_7, p2_7)
cat(sprintf("V2 new, F1 IDV7: ratio %.2f | %d of 6 lowest | %d/20 | FA %d -> %s\n",
            n7$ratio, n7$low6, n7$det, n7$fa,
            ok(round(n7$ratio, 2) == 4.02 && n7$low6 == 6 && n7$det == 20 && n7$fa == 0)))

n1 <- evalua(f1_1, p2_7)
cat(sprintf("V3 new, F1 IDV1: %d/20 | FA %d -> %s\n", n1$det, n1$fa, ok(n1$det == 20 && n1$fa == 0)))

cat(sprintf("V4 min T2 of faulty batches: F1 IDV7 %.1f (M4 41.9) | F1 IDV1 %.1f (M4 32.4) -> %s\n",
            n7$t2min, n1$t2min, ok(abs(n7$t2min - 41.9) < 3 && abs(n1$t2min - 32.4) < 3)))

f1s <- f1_7; f1s$xmeas_17 <- 100 * f1s$xmeas_17
p2s <- p2_7; p2s$xmeas_17 <- 100 * p2s$xmeas_17
sc  <- evalua(f1s, p2s)
cat(sprintf("V5 invariance (xmeas_17 x100): max weight diff %.1e | max T2 diff %.1e -> %s\n",
            max(abs(sc$cal$weights - n7$cal$weights)), max(abs(sc$t2 - n7$t2)),
            ok(max(abs(sc$cal$weights - n7$cal$weights)) < 1e-10 && max(abs(sc$t2 - n7$t2)) < 1e-6)))

again <- evalua(f1_7, p2_7)
cat(sprintf("V6 same Phase 1 twice: identical mu_r -> %s\n", ok(identical(again$cal$mu_r, n7$cal$mu_r))))

cat("\nReference center (package 0.3.0, F1 IDV7):\n"); print(round(n7$cal$mu_r, 4))
cat("Options recorded:", n7$cal$scaling, "/", n7$cal$center, "/ center_alpha", n7$cal$center_alpha, "\n")
# ---------------------------------------------------------------------
# V7 (Reviewer 1, comment 4): final method under a change of units,
#   saved to CSV. Each TEP variable, one at a time, is multiplied by 100
#   in Phase 1 and Phase 2. The calibration uses the package defaults
#   (MAD scaling + MCD of the batch centers) and the same seed as V2,
#   so the only difference between rows is the unit of one variable.
# ANCHORS
#   - The "original" row reproduces V2: ratio 4.02, 6 of 6 lowest,
#     20/20 detected, 0/10 false alarms.
#   - Every rescaled row has the same weights (< 1e-10) and the same
#     Phase 2 T2 (< 1e-6) as the original units.
# OUTPUT
#   05_revision_R1/02_resultados/table_R1_04f_TEP_final_units.csv
# ---------------------------------------------------------------------
fila <- function(unidades, r) data.frame(
  unidades             = unidades,
  ratio_sano_cont      = round(r$ratio, 2),
  fallo_en_6_menores   = r$low6,
  deteccion_de_20      = r$det,
  falsas_alarmas_de_10 = r$fa,
  T2_min_fallo         = round(r$t2min, 1),
  UCL                  = round(r$ucl, 4),
  max_dif_pesos        = max(abs(r$cal$weights - n7$cal$weights)),
  max_dif_T2           = max(abs(r$t2 - n7$t2)))

tab_unid <- fila("original", n7)
for (v in VARS) {
  f1s <- f1_7; f1s[[v]] <- 100 * f1s[[v]]
  p2s <- p2_7; p2s[[v]] <- 100 * p2s[[v]]
  tab_unid <- rbind(tab_unid, fila(paste(v, "x100"), evalua(f1s, p2s)))
}
print(tab_unid, row.names = FALSE)

stopifnot(tab_unid$ratio_sano_cont[1] == 4.02,
          tab_unid$fallo_en_6_menores[1] == 6,
          tab_unid$deteccion_de_20[1] == 20,
          tab_unid$falsas_alarmas_de_10[1] == 0,
          all(tab_unid$max_dif_pesos < 1e-10),
          all(tab_unid$max_dif_T2 < 1e-6))

write.csv(tab_unid,
          "05_revision_R1/02_resultados/table_R1_04f_TEP_final_units.csv",
          row.names = FALSE)
cat("V7 saved: table_R1_04f_TEP_final_units.csv\n")
