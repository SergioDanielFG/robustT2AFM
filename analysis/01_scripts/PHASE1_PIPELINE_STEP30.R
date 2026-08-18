# =====================================================================
# PHASE1_PIPELINE_STEP30.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The base pipeline of the Tennessee Eastman application: it builds
#   the contaminated Phase 1 and the Phase 2 of the manuscript, calibrates
#   both methods, and reports the two control limits, the faulty batches
#   detected out of twenty, the false alarms out of ten, the inflation of
#   the classical covariance determinant relative to the robust one, the
#   healthy / contaminated AFM weight ratio, the sorted weights and the
#   Phase 2 T2 values of each method.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - Batches are built by spaced subsampling, one observation every
#     PASO = 30, which is the smallest of the steps examined that leaves
#     all four variables independent under Ljung-Box. A batch therefore
#     spans 1 + (I-1)*30 = 571 samples of the run.
#   - The healthy runs come from the fault-free TEST set, not the
#     training set: a step-30 batch needs 571 observations and the
#     training runs are shorter than that.
#   - Faulty batches start at DESDE_FALLO, once the fault is already
#     present in the run.
#   - Phase 1 holds 24 healthy batches and 6 with IDV 7, so the
#     historical sample is contaminated by whole anomalous batches;
#     Phase 2 holds 10 healthy and 20 faulty batches and reuses none of
#     the Phase 1 runs.
#   - calibrate_afm_mcd() is called without fixing a seed, because
#     covMcd() explores subsets at random; the stability of the reported
#     figures across seeds is measured separately in
#     TEP_BOOTSTRAP_LIMIT.R.
#
# ANCHORS AND THEIR EXPECTED VALUES
#   robust UCL                  = 19.69   (tolerance 0.10)
#   classical UCL               = 19.46   (tolerance 0.10)
#   robust detection            = 20 of 20 (exact)
#   classical detection         =  2 of 20 (tolerance 1)
#   robust false alarms         =  0 of 10 (exact)
#   classical false alarms      =  0 of 10 (exact)
#   covariance inflation        = 16       (tolerance 3)
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   All seven checks must report OK; only then is Phase 1 considered
#   anchored and its objects usable by the scripts that depend on it.
#
# OUTPUT
#   No CSV file. The script leaves in memory the objects the control
#   chart and AFM weight figures consume: tr, tc, ucl_rob, ucl_cls, ef
#   and cal_rob.
# =====================================================================
library(robustT2AFM)
DIR_DATOS <- "04_datos"

VARS   <- c("xmv_9", "xmv_8", "xmeas_19", "xmeas_17")
I_LOTE <- 20; PASO <- 30; J <- 4
ALPHA  <- 0.001; H_MCD <- 0.67; FALLO <- 7
DESDE_SANO  <- 1
DESDE_FALLO <- 161
muestras_esp <- function(d) seq(d, by = PASO, length.out = I_LOTE)

if (!exists("ff_test")) { cat("Loading FaultFree_Testing...\n"); ff_test <- read.csv(file.path(DIR_DATOS, "TEP_FaultFree_Testing.csv")) }
if (!exists("fy"))      { cat("Loading Faulty_Testing...\n");    fy <- read.csv(file.path(DIR_DATOS, "TEP_Faulty_Testing.csv")) }

extraer_lote <- function(datos, run_id, fault_id, desde, batch_name) {
  sel <- datos$simulationRun==run_id & datos$faultNumber==fault_id &
         datos$sample %in% muestras_esp(desde)
  sub <- datos[sel, VARS]
  if (nrow(sub) != I_LOTE) return(NULL)
  data.frame(Batch = batch_name, sub, row.names = NULL)
}

# --- PHASE 1: 24 healthy (Testing 1-24) + 6 with IDV 7 (Faulty 1-6) ---
f1 <- do.call(rbind, c(
  lapply(1:24, function(r) extraer_lote(ff_test, r, 0, DESDE_SANO,  sprintf("F1_sano_%02d", r))),
  lapply(1:6,  function(r) extraer_lote(fy,      r, FALLO, DESDE_FALLO, sprintf("F1_fallo_%02d", r)))
))
# --- PHASE 2: 10 healthy (Testing 25-34) + 20 with IDV 7 (Faulty 7-26) ---
f2 <- do.call(rbind, c(
  lapply(25:34, function(r) extraer_lote(ff_test, r, 0, DESDE_SANO,  sprintf("F2_sano_%02d", r))),
  lapply(7:26,  function(r) extraer_lote(fy,      r, FALLO, DESDE_FALLO, sprintf("F2_fallo_%03d", r)))
))
cat(sprintf("Phase 1: %d batches | Phase 2: %d batches\n",
            length(unique(f1$Batch)), length(unique(f2$Batch))))
stopifnot(length(unique(f1$Batch))==30, length(unique(f2$Batch))==30)

# --- Calibration ---
cal_rob <- calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD)
cal_cls <- hotelling_classical_calibrate(f1, VARS)
K <- length(unique(f1$Batch))
ucl_rob <- ucl_F_adjusted(cal_rob, I = I_LOTE, alpha = ALPHA)$UCL
ucl_cls <- hotelling_classical_ucl(K = K, I = I_LOTE, J = J, alpha = ALPHA, phase = "II")$UCL

# --- Phase 2 detection ---
b  <- unique(f2$Batch); ef <- grepl("fallo", b)
mr <- monitor_afm_mcd(f2, cal_rob, VARS); mc <- hotelling_classical_monitor(f2, cal_cls, VARS)
tr <- mr$T2[match(b, mr$Batch)]; tc <- mc$T2[match(b, mc$Batch)]
det_rob <- sum(tr[ef]>ucl_rob);  det_cls <- sum(tc[ef]>ucl_cls)
fa_rob  <- sum(tr[!ef]>ucl_rob); fa_cls  <- sum(tc[!ef]>ucl_cls)

# --- Masking metrics ---
inflacion <- det(cal_cls$Sp)/det(cal_rob$Sw)
w <- cal_rob$weights
idx_cont <- grepl("fallo", names(w))
ratio_pesos <- mean(w[!idx_cont])/mean(w[idx_cont])

# =====================================================================
# RESULTS + VERIFICATION (anchor of the application chapter)
# =====================================================================
cat("\n================ BASE RESULT, STEP 30 (IDV 7) ================\n")
cat(sprintf("  Robust UCL: %.2f  | classical UCL: %.2f\n", ucl_rob, ucl_cls))
cat(sprintf("  Robust detection: %d/20  | classical: %d/20\n", det_rob, det_cls))
cat(sprintf("  Robust false alarms: %d/10 | classical: %d/10\n", fa_rob, fa_cls))
cat(sprintf("  Covariance inflation: %.2f\n", inflacion))
cat(sprintf("  Healthy/contaminated weight ratio: %.1f\n", ratio_pesos))

cat("\n--- Sorted AFM weights (the 6 contaminated must be the lowest) ---\n")
print(round(sort(w), 4))

cat("\n--- Phase 2 T2 ---\n")
cat("Rob faulty :", round(tr[ef],1), "\n")
cat("Rob healthy:", round(tr[!ef],1), "| UCL", round(ucl_rob,2), "\n")
cat("Cls faulty :", round(tc[ef],1), "\n")
cat("Cls healthy:", round(tc[!ef],1), "| UCL", round(ucl_cls,2), "\n")

cat("\n===== VERIFICATION =====\n")
chk <- function(n,o,e,tol) cat(sprintf("  %-22s obs=%-8.2f exp=%-8.2f %s\n",
                                n,o,e,ifelse(abs(o-e)<=tol,"OK","REVIEW")))
chk("Robust UCL", ucl_rob, 19.69, 0.1)
chk("Classical UCL", ucl_cls, 19.46, 0.1)
chk("Robust detection", det_rob, 20, 0)
chk("Classical detection", det_cls, 2, 1)
chk("Robust false alarms", fa_rob, 0, 0)
chk("Classical false alarms", fa_cls, 0, 0)
chk("Inflation", inflacion, 16, 3)

cat("\nIf everything reports OK, Phase 1 is anchored and its objects can be\n")
cat("used by the figures and by the multiple-fault validation.\n")
