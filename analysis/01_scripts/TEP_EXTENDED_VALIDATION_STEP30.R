# =====================================================================
# TEP_EXTENDED_VALIDATION_STEP30.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The validation of the application over five Tennessee Eastman faults
#   of deliberately different nature: a strong step (IDV 1), the
#   moderate step of the main analysis (IDV 7), a slow drift (IDV 13), a
#   random disturbance (IDV 8) and a fault known in the literature to be
#   hard to detect (IDV 3). For each one it reports the mean shift it
#   induces, the faulty batches detected out of twenty and the
#   percentage of false alarms over ten healthy batches, for both
#   methods.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - Phase 1 is FIXED and contaminated with IDV 7, and is calibrated
#     ONCE. The only thing that varies between rows is the fault each
#     method has to detect, so the comparison isolates the fault type.
#   - The ten healthy Phase 2 batches are the same in every row.
#   - The mean shift is measured as the largest change across the four
#     variables, in standard deviations of the pre-fault stretch of the
#     same runs, so the shift and the detection refer to the same data.
#   - Batches are built by spaced subsampling, one observation every
#     PASO = 30, so the observations of a batch can be treated as
#     independent in time; faulty batches start at DESDE_FALLO, once the
#     fault is already present.
#   - Rows are sorted by decreasing shift, which is what makes the
#     pattern of the classical method readable.
#
# ANCHOR AND ITS EXPECTED VALUE
#   The IDV 7 row must reproduce the published application: 20 of 20
#   detected with the robust method and 2 of 20 with the classical one.
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   The advantage of the proposed method is confirmed as independent of
#   the chosen fault only if it detects every faulty batch in all the
#   faults that shift the mean vector, with no false alarms. A fault
#   that barely moves the mean is outside the scope of the method and
#   its result is reported as such, not as a failure.
#
# OUTPUT
#   02_resultados/table_multiple_faults_step30.csv
# =====================================================================
library(robustT2AFM)
DIR_DATOS <- "04_datos"
DIR_OUT   <- "02_resultados"
VARS   <- c("xmv_9", "xmv_8", "xmeas_19", "xmeas_17")
I_LOTE <- 20; PASO <- 30; J <- 4
ALPHA  <- 0.001; H_MCD <- 0.67
DESDE_SANO  <- 1
DESDE_FALLO <- 161
N_SANO_F2   <- 10
muestras_esp <- function(d) seq(d, by = PASO, length.out = I_LOTE)
if (!exists("ff_test")) ff_test <- read.csv(file.path(DIR_DATOS, "TEP_FaultFree_Testing.csv"))
if (!exists("fy"))      fy      <- read.csv(file.path(DIR_DATOS, "TEP_Faulty_Testing.csv"))
extraer <- function(datos, run, fault, desde, nm) {
sel <- datos$simulationRun==run & datos$faultNumber==fault & datos$sample %in% muestras_esp(desde)
sub <- datos[sel, VARS]; if (nrow(sub)!=I_LOTE) return(NULL)
data.frame(Batch=nm, sub, row.names=NULL)
}
# --- FIXED PHASE 1: 24 healthy (Testing) + 6 IDV 7. Calibrate once. ---
f1 <- do.call(rbind, c(
lapply(1:24, function(r) extraer(ff_test, r, 0, DESDE_SANO, sprintf("F1s%02d",r))),
lapply(1:6,  function(r) extraer(fy, r, 7, DESDE_FALLO, sprintf("F1f%02d",r)))))
cal_rob <- calibrate_afm_mcd(f1, VARS, mcd_alpha=H_MCD)
cal_cls <- hotelling_classical_calibrate(f1, VARS)
K <- length(unique(f1$Batch))
ucl_rob <- ucl_F_adjusted(cal_rob, I=I_LOTE, alpha=ALPHA)$UCL
ucl_cls <- hotelling_classical_ucl(K=K, I=I_LOTE, J=J, alpha=ALPHA, phase="II")$UCL
cat(sprintf("UCL_rob=%.2f UCL_cls=%.2f\n", ucl_rob, ucl_cls))
sanos_f2 <- do.call(rbind, lapply(25:34, function(r) extraer(ff_test, r, 0, DESDE_SANO, sprintf("F2s%02d",r))))
# --- Detection of one fault in Phase 2 (20 batches, runs 7-26) ---
evaluar <- function(fault_id) {
fallo_f2 <- do.call(rbind, lapply(7:26, function(r) extraer(fy, r, fault_id, DESDE_FALLO, sprintf("F2f%03d",r))))
p2 <- rbind(sanos_f2, fallo_f2)
mr <- monitor_afm_mcd(p2, cal_rob, VARS); mc <- hotelling_classical_monitor(p2, cal_cls, VARS)
b <- unique(p2$Batch); ef <- grepl("F2f", b)
tr <- mr$T2[match(b, mr$Batch)]; tc <- mc$T2[match(b, mc$Batch)]
c(det_rob=sum(tr[ef]>ucl_rob), det_cls=sum(tc[ef]>ucl_cls),
far_rob=100*sum(tr[!ef]>ucl_rob)/N_SANO_F2, far_cls=100*sum(tc[!ef]>ucl_cls)/N_SANO_F2)
}
# --- Mean shift of one fault (same method as the IDV 3 check) ---
desplaz <- function(fault_id) {
ms <- muestras_esp(DESDE_FALLO)
normal  <- fy[fy$faultNumber==fault_id & fy$simulationRun %in% 7:26 & fy$sample %in% 1:150, VARS]
fallado <- fy[fy$faultNumber==fault_id & fy$simulationRun %in% 7:26 & fy$sample %in% ms, VARS]
max(abs((colMeans(fallado)-colMeans(normal))/apply(normal,2,sd)))
}
fallos <- c(1, 7, 13, 8, 3)
tipo   <- c("Escalon fuerte","Escalon moderado","Deriva lenta","Aleatorio","Casi indetectable")
cat("\nEvaluating 5 faults with step 30...\n")
tabla <- data.frame()
for (i in seq_along(fallos)) {
d <- evaluar(fallos[i]); dm <- desplaz(fallos[i])
tabla <- rbind(tabla, data.frame(
IDV=fallos[i], Type=tipo[i], Shift_SD=round(dm,2),
Det_Rob=d["det_rob"], Det_Cls=d["det_cls"],
FAR_Rob=d["far_rob"], FAR_Cls=d["far_cls"]))
}
rownames(tabla) <- NULL
# Sort by decreasing shift, as in the published table
tabla <- tabla[order(-tabla$Shift_SD), ]
cat("\n================ MULTIPLE-FAULT TABLE (step 30) ================\n")
print(tabla, row.names=FALSE)
# --- IDV 7 anchor ---
a <- tabla[tabla$IDV==7, ]
cat(sprintf("\n--- ANCHOR IDV 7: rob=%d (exp 20) %s | cls=%d (exp 2) %s ---\n",
a$Det_Rob, ifelse(a$Det_Rob==20,"OK","REVIEW"),
a$Det_Cls, ifelse(a$Det_Cls==2,"OK","REVIEW")))
write.csv(tabla, file.path(DIR_OUT, "table_multiple_faults_step30.csv"), row.names=FALSE)
cat("\nSaved to 02_resultados/table_multiple_faults_step30.csv\n")
