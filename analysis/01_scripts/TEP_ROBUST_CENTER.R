# =====================================================================
# TEP_ROBUST_CENTER.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The three shift figures of the reference-center paragraph of the
#   Tennessee Eastman application, per variable and as a maximum:
#     (a) shift of the ARITHMETIC MEAN of the faulty batches with
#         respect to the center of the healthy batches
#     (b) shift of the MCD CENTER of those same batches
#     (c) shift of the reference center mu_r estimated with the
#         contaminated Phase 1 against the one estimated with the
#         healthy batches only
#   It also reports the recentering factor (a)/(b), which is the actual
#   measurement of what the MCD corrects inside each batch.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - All shifts are expressed in standard deviations of the healthy
#     observations, so the four stripper variables, which live on very
#     different scales, are comparable.
#   - The healthy reference center is the average of the MCD centers of
#     the 24 healthy batches, not the raw mean, so (a) and (b) are
#     measured against the same reference.
#   - Batches are built by spaced subsampling, one observation every
#     PASO = 30, so the observations of a batch can be treated as
#     independent in time.
#   - Quantity (c) equals (b)/5 by construction, because 6 contaminated
#     batches out of 30 enter a simple average with weight 1/5. It is
#     arithmetic, not an independent check; the real measurement is (a)
#     against (b). The script prints that verification explicitly.
#   - Nothing here depends on the previous state of the R session.
#
# ANCHOR AND ITS EXPECTED VALUE
#   Robust UCL = 19.69285 (tolerance 1e-04), checked before computing
#   anything else.
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   The script stops if the UCL anchor fails; the calibration would not
#   be the one of the manuscript and the three figures would not refer
#   to the published application.
#
# OUTPUT
#   02_resultados/tep_robust_center_by_variable.csv
#   02_resultados/tep_robust_center_summary.csv
# =====================================================================
library(robustT2AFM)
DIR_DATOS <- "04_datos"
DIR_OUT   <- "02_resultados"
VARS   <- c("xmv_9", "xmv_8", "xmeas_19", "xmeas_17")
I_LOTE <- 20; PASO <- 30; J <- 4
ALPHA  <- 0.001; H_MCD <- 0.67
DESDE_SANO <- 1; DESDE_FALLO <- 161
muestras_esp <- function(d) seq(d, by = PASO, length.out = I_LOTE)
cat("Loading data...\n")
ff <- read.csv(file.path(DIR_DATOS, "TEP_FaultFree_Testing.csv"))
fy <- read.csv(file.path(DIR_DATOS, "TEP_Faulty_Testing.csv"))
extraer_lote <- function(datos, run_id, fault_id, desde, batch_name) {
sel <- datos$simulationRun == run_id & datos$faultNumber == fault_id &
datos$sample %in% muestras_esp(desde)
sub <- datos[sel, VARS]
if (nrow(sub) != I_LOTE) return(NULL)
data.frame(Batch = batch_name, sub, row.names = NULL)
}
# --- Phase 1: 24 healthy + 6 with IDV 7 ------------------------------
f1 <- do.call(rbind, c(
lapply(1:24, function(r) extraer_lote(ff, r, 0, DESDE_SANO,  sprintf("F1_sano_%02d", r))),
lapply(1:6,  function(r) extraer_lote(fy, r, 7, DESDE_FALLO, sprintf("F1_fallo_%02d", r)))
))
stopifnot(nrow(f1) == 30 * I_LOTE)
cal_cont <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
# --- ANCHOR ----------------------------------------------------------
ucl <- ucl_F_adjusted(cal_cont, I = I_LOTE, alpha = ALPHA)$UCL
cat(sprintf("\nANCHOR  UCL = %.5f  (expected 19.69285)  %s\n\n",
ucl, if (abs(ucl - 19.69285) < 1e-4) "OK" else "*** FAILS ***"))
if (abs(ucl - 19.69285) > 1e-4) stop("Anchor failed. STOP.")
# --- References ------------------------------------------------------
sanos_nm <- unique(f1$Batch[grepl("sano",  f1$Batch)])
fallo_nm <- unique(f1$Batch[grepl("fallo", f1$Batch)])
# scale: sd of the observations of the 24 healthy batches
sd_ref <- apply(f1[f1$Batch %in% sanos_nm, VARS], 2, sd)
# healthy reference center = average of the MCD centers of the 24 healthy batches
c_sano <- colMeans(do.call(rbind, cal_cont$mcd_centers[sanos_nm]))
# --- (a) and (b): simple mean against MCD center of the faulty batches
med_simple <- t(vapply(fallo_nm,
function(k) colMeans(f1[f1$Batch == k, VARS]),
numeric(J)))
cen_mcd <- do.call(rbind, cal_cont$mcd_centers[fallo_nm])
d_simple <- abs(colMeans(med_simple) - c_sano) / sd_ref
d_mcd    <- abs(colMeans(cen_mcd)    - c_sano) / sd_ref
# --- (c): contaminated mu_r against mu_r with healthy batches only ---
cal_sano <- suppressMessages(
calibrate_afm_mcd(f1[f1$Batch %in% sanos_nm, ], VARS, mcd_alpha = H_MCD)
)
d_mur <- abs(cal_cont$mu_r - cal_sano$mu_r) / sd_ref
# =====================================================================
#  OUTPUT
# =====================================================================
tabla <- data.frame(
Variable                = VARS,
reference_sd            = round(as.numeric(sd_ref), 4),
d_simple_mean_SD        = round(as.numeric(d_simple), 4),
d_MCD_center_SD         = round(as.numeric(d_mcd), 4),
d_reference_center_SD   = round(as.numeric(d_mur), 4),
row.names = NULL
)
cat("=========== SHIFTS, BY VARIABLE (in SD) ===========\n")
print(tabla, row.names = FALSE)
a <- max(d_simple); b <- max(d_mcd); cc <- max(d_mur)
cat("\n=========== THE THREE FIGURES OF THE PARAGRAPH ===========\n")
cat(sprintf("  (a) arithmetic mean of the faulty batches : %.2f SD\n", a))
cat(sprintf("  (b) MCD center of those same batches      : %.2f SD\n", b))
cat(sprintf("      recentering factor (a/b)              : %.2f times\n", a / b))
cat(sprintf("  (c) reference center mu_r                 : %.2f SD\n", cc))
cat("\n  NOTE: (c) = (b)/5 by construction, because 6 batches out of 30\n")
cat("  enter the average with weight 1/5. That is arithmetic, not an\n")
cat("  independent check. The real measurement is (a) against (b).\n")
cat(sprintf("  verification: (b)/5 = %.4f   (c) = %.4f\n", b / 5, cc))
resumen <- data.frame(
concept = c("media simple lotes con fallo",
"centro MCD lotes con fallo",
"factor de recentrado",
"centro de referencia mu_r"),
value_SD = c(round(a, 4), round(b, 4), round(a / b, 4), round(cc, 4)),
variable_of_maximum = c(VARS[which.max(d_simple)], VARS[which.max(d_mcd)],
"-", VARS[which.max(d_mur)])
)
write.csv(tabla,   file.path(DIR_OUT, "tep_robust_center_by_variable.csv"), row.names = FALSE)
write.csv(resumen, file.path(DIR_OUT, "tep_robust_center_summary.csv"),      row.names = FALSE)
cat("\nSaved to:\n")
cat("  02_resultados/tep_robust_center_by_variable.csv\n")
cat("  02_resultados/tep_robust_center_summary.csv\n")
cat("\n=========== EXPECTED VALUES ===========\n")
cat("  (a) 1.10   (b) 0.31   (c) 0.06\n")
cat("  If they do not match, do NOT use the published ones until clarified.\n\n")
