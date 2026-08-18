# =====================================================================
# SUBSAMPLING_STEP_SENSITIVITY.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The sensitivity of the application to the subsampling step. For each
#   step (10, 20, 30) and each of the four stripper variables it reports
#   the lag-1 autocorrelation of the subsampled series and the Ljung-Box
#   test of temporal independence. It then rebuilds the whole Phase 1
#   and Phase 2 at each step and reports, for both methods, the two
#   control limits, the faulty batches detected out of twenty and the
#   false alarms out of ten.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - The question "why step 30 and not less" is answered by the
#     Ljung-Box test, which is the same test reported in the manuscript,
#     and not by an ARL0: with only ten in-control Phase 2 batches an
#     ARL0 of 1000 cannot be estimated, only false alarms can be
#     counted.
#   - Nothing new is simulated: the same runs as the manuscript are
#     used, so the comparison across steps isolates the step.
#   - A step is only feasible if a batch of I_LOTE observations still
#     fits inside the run; the script checks that before evaluating it.
#   - The Ljung-Box test uses LB_LAG = 10 lags at every step.
#
# ANCHOR AND ITS EXPECTED VALUE
#   Step 30 must reproduce the published application: robust
#   UCL = 19.69 (tolerance 0.02), 20 of 20 faulty batches detected and
#   0 false alarms.
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   A step is considered sufficient if the Ljung-Box test does NOT
#   reject (p > 0.05) in ALL FOUR variables. Registered prediction:
#   steps 10 and 20 reject in at least one variable, step 30 rejects in
#   none. If it comes out otherwise, the sentence of the manuscript
#   stating that a shorter step would not suffice has to be revised.
#
# OUTPUT
#   02_resultados/step_sensitivity_ljungbox.csv
#   02_resultados/step_sensitivity_detection.csv
# =====================================================================
library(robustT2AFM)
DIR_DATOS <- "04_datos"
DIR_OUT   <- "02_resultados"
VARS   <- c("xmv_9", "xmv_8", "xmeas_19", "xmeas_17")
I_LOTE <- 20
J      <- 4
ALPHA  <- 0.001
H_MCD  <- 0.67
DESDE_SANO  <- 1
DESDE_FALLO <- 161
PASOS  <- c(10, 20, 30)
LB_LAG <- 10          # lags of the Ljung-Box test
cat("Loading data...\n")
ff <- read.csv(file.path(DIR_DATOS, "TEP_FaultFree_Testing.csv"))
fy <- read.csv(file.path(DIR_DATOS, "TEP_Faulty_Testing.csv"))
max_obs <- max(ff$sample)
cat("Samples per run:", max_obs, "\n\n")
# =====================================================================
#  PART 1 - Temporal independence at each step
# =====================================================================
cat("#########################################################\n")
cat("#  PART 1 - LJUNG-BOX AT EACH STEP\n")
cat("#########################################################\n")
serie_sana <- ff[ff$simulationRun == 1 & ff$faultNumber == 0, ]
serie_sana <- serie_sana[order(serie_sana$sample), ]
filas <- list()
for (paso in PASOS) {
idx <- seq(1, nrow(serie_sana), by = paso)
sub <- serie_sana[idx, VARS]
tramo <- 1 + (I_LOTE - 1) * paso
cat(sprintf("\n--- STEP %d  (a batch spans %d samples; %s) ---\n",
paso, tramo,
if (tramo <= max_obs) "fits in the run" else "DOES NOT FIT"))
cat(sprintf("  observations available after subsampling: %d\n\n", nrow(sub)))
for (v in VARS) {
x  <- sub[[v]]
a1 <- as.numeric(acf(x, lag.max = 1, plot = FALSE)$acf[2])
lb <- Box.test(x, lag = LB_LAG, type = "Ljung-Box")
veredicto <- if (lb$p.value > 0.05) "independent" else "*** AUTOCORRELATED ***"
cat(sprintf("  %-10s ACF(1) = %7.3f   Ljung-Box p = %.4f   %s\n",
v, a1, lb$p.value, veredicto))
filas[[length(filas) + 1]] <- data.frame(
step = paso, variable = v,
acf1 = round(a1, 4), p_ljungbox = round(lb$p.value, 4),
independent = lb$p.value > 0.05
)
}
}
tabla_lb <- do.call(rbind, filas)
cat("\n\n===== SUMMARY BY STEP =====\n")
for (paso in PASOS) {
s <- tabla_lb[tabla_lb$step == paso, ]
n_ok <- sum(s$independent)
cat(sprintf("  step %2d : %d of 4 variables independent   %s\n",
paso, n_ok,
if (n_ok == 4) "<-- SUFFICIENT" else "insufficient"))
}
# =====================================================================
#  PART 2 - Detection and false alarms at each step
# =====================================================================
cat("\n\n#########################################################\n")
cat("#  PART 2 - FULL PHASE 2 AT EACH STEP\n")
cat("#########################################################\n")
evaluar_paso <- function(paso) {
muestras <- function(d) seq(d, by = paso, length.out = I_LOTE)
if (DESDE_FALLO + (I_LOTE - 1) * paso > max_obs) return(NULL)
lote <- function(dat, run, fault, desde, nom) {
sel <- dat$simulationRun == run & dat$faultNumber == fault &
dat$sample %in% muestras(desde)
s <- dat[sel, VARS]
if (nrow(s) != I_LOTE) return(NULL)
data.frame(Batch = nom, s, row.names = NULL)
}
f1 <- do.call(rbind, c(
lapply(1:24, function(r) lote(ff, r, 0, DESDE_SANO,  sprintf("F1_sano_%02d", r))),
lapply(1:6,  function(r) lote(fy, r, 7, DESDE_FALLO, sprintf("F1_fallo_%02d", r)))
))
if (is.null(f1) || nrow(f1) != 30 * I_LOTE) return(NULL)
cal_r <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
cal_c <- hotelling_classical_calibrate(f1, VARS)
ucl_r <- ucl_F_adjusted(cal_r, I = I_LOTE, alpha = ALPHA)$UCL
ucl_c <- hotelling_classical_ucl(K = 30, I = I_LOTE, J = J,
alpha = ALPHA, phase = "II")$UCL
p2 <- rbind(
do.call(rbind, lapply(25:34, function(r) lote(ff, r, 0, DESDE_SANO,  sprintf("F2_sano_%02d", r)))),
do.call(rbind, lapply(7:26,  function(r) lote(fy, r, 7, DESDE_FALLO, sprintf("F2_fallo_%03d", r))))
)
mr <- monitor_afm_mcd(p2, cal_r, VARS)
mc <- hotelling_classical_monitor(p2, cal_c, VARS)
bf <- unique(p2$Batch); ef <- grepl("fallo", bf)
tr <- mr$T2[match(bf, mr$Batch)]; tc <- mc$T2[match(bf, mc$Batch)]
data.frame(
step = paso,
UCL_rob = round(ucl_r, 2), UCL_cls = round(ucl_c, 2),
Det_Rob = sum(tr[ef] > ucl_r), Det_Cls = sum(tc[ef] > ucl_c),
FA_Rob  = sum(tr[!ef] > ucl_r), FA_Cls = sum(tc[!ef] > ucl_c),
T2med_Rob = round(median(tr[ef]), 1)
)
}
res <- do.call(rbind, lapply(PASOS, evaluar_paso))
cat("\n")
print(res, row.names = FALSE)
cat("\n  Det: faulty batches detected out of 20 | FA: false alarms out of 10\n")
# =====================================================================
#  OUTPUT
# =====================================================================
write.csv(tabla_lb, file.path(DIR_OUT, "step_sensitivity_ljungbox.csv"), row.names = FALSE)
write.csv(res,      file.path(DIR_OUT, "step_sensitivity_detection.csv"), row.names = FALSE)
cat("\nSaved to:\n")
cat("  02_resultados/step_sensitivity_ljungbox.csv\n")
cat("  02_resultados/step_sensitivity_detection.csv\n")
cat("\n===== ANCHOR =====\n")
a <- res[res$step == 30, ]
if (nrow(a) == 1) {
ok <- abs(a$UCL_rob - 19.69) < 0.02 && a$Det_Rob == 20 && a$FA_Rob == 0
cat(sprintf("  step 30 -> UCL %.2f, %d of 20 detected, %d false alarms  %s\n",
a$UCL_rob, a$Det_Rob, a$FA_Rob,
if (ok) "MATCHES THE MANUSCRIPT" else "*** DOES NOT MATCH, REVIEW ***"))
}
cat("\n=== END. Copy the whole output. ===\n\n")
