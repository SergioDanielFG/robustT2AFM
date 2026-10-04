# =====================================================================
# SCRIPT_M4c_TEP_center_weights.R
# ---------------------------------------------------------------------
# Revision round 1, module M4, third part (P17, Section 4.4; P20,
# Section 4.7; Figure 5. Reviewer 1, comments 4, 6 and 7).
#
# WHY
#   Section 4.4 quotes figures that come from TEP_ROBUST_CENTER.R, run
#   with robustT2AFM 0.2.0 and without a fixed seed:
#     - shift of the arithmetic mean of the six faulty batches (1.10 and
#       0.84 SD), of their MCD centers (0.31 and 0.26 SD), and the other
#       two variables below 0.21 SD (tep_robust_center_by_variable.csv);
#     - Figure 5 and the sentence "all below the weight that would
#       correspond to them if all contributed equally" (no CSV holds the
#       weights batch by batch).
#   In 0.3.0 the MCD of each batch is computed on the original data and
#   the MAD scaling enters only the first eigenvalue, so the per-batch
#   centers should not change; this script checks it with the seed of
#   the final method instead of assuming it. It also records which batch
#   centers the MCD of the centers (the new reference center, Eq. (6))
#   keeps, to explain the drift of mu_r in Sections 4.4 and 4.7.
#
# WHAT IT COMPUTES (nothing is simulated)
#   For Phase 1 contaminated with IDV 7 and with IDV 1:
#   A. Shift per variable, in SD of the healthy observations, of
#      (a) the arithmetic mean of the six faulty batches, (b) their MCD
#      centers, (c) mu_r of V7 and (d) mu_r of NEW, each against the
#      same estimator applied to the 24 healthy batches only (as in
#      TEP_ROBUST_CENTER.R and table_M4_final_E_center.csv).
#   B. MFA weight of every batch, V7 and NEW, against the uniform 1/30.
#   C. MCD of the 30 batch centers (center_alpha = 0.5, internal seed of
#      the package): whether each center is in the raw MCD subset and
#      its reweighting weight (0 = excluded from mu_r).
#
# DESIGN
#   - Batches, runs, starting samples and step as in
#     SCRIPT_M4_TEP_final_mstar14.R; set.seed(2026) before every
#     calibration, as in M4 final and M4b.
#   - V7 = calibrate_afm_mcd(scaling = "none", center = "mean");
#     NEW = calibrate_afm_mcd() with the 0.3.0 defaults.
#   - Values are saved without rounding.
#
# ANCHORS (STOP if any fails)
#   B1  UCL V7 19.69285 (m* = 13), NEW 19.64480 (m* = 14)
#   B2  IDV 7, (a) per variable = tep_robust_center_by_variable.csv
#       (1.1033 / 0.0598 / 0.8422 / 0.1914), tolerance 0.005
#   B3  weight ratios = table_M4b_D2_D4.csv: IDV 7 6.4394 / 4.0226,
#       IDV 1 7.0851 / 4.3089 (4 decimals)
#   B4  maximum drift of mu_r = table_M4_final_E_center.csv:
#       V7 0.0614 (IDV 7), 0.8167 (IDV 1); NEW 0.0859 (IDV 7),
#       0.0265 (IDV 1) (4 decimals)
#   B5  the MCD of the centers recomputed here equals mu_r of NEW
#       (max abs difference < 1e-10)
# CHECK (printed, does not stop)
#   C1  IDV 7, (b) per variable vs tep_robust_center_by_variable.csv
#       (0.307 / 0.0926 / 0.2627 / 0.2048)
#
# OUTPUT (05_revision_R1/02_resultados/)
#   table_M4c_center_by_variable.csv   part A
#   table_M4c_weights.csv              part B (input of Figure 5)
#   table_M4c_center_membership.csv    part C
#   Run from C:/temp_paper. Expected time: seconds.
# =====================================================================
library(robustT2AFM)
library(robustbase)
stopifnot(packageVersion("robustT2AFM") >= "0.3.0")
DIR_DATOS <- "04_datos"
DIR_OUT   <- "05_revision_R1/02_resultados"
VARS   <- c("xmv_9", "xmv_8", "xmeas_19", "xmeas_17")
I_LOTE <- 20; PASO <- 30; J <- 4
ALPHA  <- 0.001; H_MCD <- 0.67
DESDE_SANO <- 1; DESDE_FALLO <- 161
SEED_CAL <- 2026; CENTER_SEED <- 20260927L      # the latter as in R/utils_seed.R

if (!exists("ff_test")) { cat("Loading FaultFree_Testing...\n"); ff_test <- read.csv(file.path(DIR_DATOS, "TEP_FaultFree_Testing.csv")) }
if (!exists("fy"))      { cat("Loading Faulty_Testing...\n");    fy      <- read.csv(file.path(DIR_DATOS, "TEP_Faulty_Testing.csv")) }

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
cal <- function(f1, ...) { set.seed(SEED_CAL); suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD, ...)) }
ratio_w <- function(w) { ic <- grepl("fallo", names(w)); mean(w[!ic]) / mean(w[ic]) }
same <- function(x, v, d) abs(round(x, d) - v) < 1e-9

partA <- list(); partB <- list(); partC <- list(); chk <- list()
for (ff in c(7, 1)) {
  f1   <- fase1(ff); stopifnot(nrow(f1) == 30 * I_LOTE)
  sano <- unique(f1$Batch[grepl("sano",  f1$Batch)])
  fall <- unique(f1$Batch[grepl("fallo", f1$Batch)])
  fs   <- f1[f1$Batch %in% sano, ]
  v7  <- cal(f1, scaling = "none", center = "mean"); new  <- cal(f1)
  v7s <- cal(fs, scaling = "none", center = "mean"); news <- cal(fs)
  sd_ref <- apply(fs[, VARS], 2, sd)
  c_sano <- colMeans(do.call(rbind, v7$mcd_centers[sano]))      # as TEP_ROBUST_CENTER.R
  d_simple <- abs(colMeans(t(vapply(fall, function(k) colMeans(f1[f1$Batch == k, VARS]), numeric(J)))) - c_sano) / sd_ref
  d_mcd    <- abs(colMeans(do.call(rbind, v7$mcd_centers[fall])) - c_sano) / sd_ref
  d_v7     <- abs(v7$mu_r  - v7s$mu_r)  / sd_ref
  d_new    <- abs(new$mu_r - news$mu_r) / sd_ref
  partA[[length(partA) + 1]] <- data.frame(phase1 = paste0("IDV ", ff), variable = VARS, sd_ref = sd_ref,
    d_simple_mean = d_simple, d_MCD_centers = d_mcd, d_mu_r_V7 = d_v7, d_mu_r_NEW = d_new, row.names = NULL)
  same_centers <- max(abs(do.call(rbind, v7$mcd_centers) - do.call(rbind, new$mcd_centers)))
  b <- names(new$weights)
  partB[[length(partB) + 1]] <- data.frame(phase1 = paste0("IDV ", ff), batch = b,
    contaminated = grepl("fallo", b), w_V7 = v7$weights[b], w_NEW = new$weights[b], uniform = 1 / length(b),
    below_uniform_V7 = v7$weights[b] < 1 / length(b), below_uniform_NEW = new$weights[b] < 1 / length(b), row.names = NULL)
  C <- do.call(rbind, new$mcd_centers)[, VARS, drop = FALSE]
  set.seed(CENTER_SEED); mc <- robustbase::covMcd(C, alpha = 0.5)
  raw_in <- seq_len(nrow(C)) %in% mc$best
  partC[[length(partC) + 1]] <- data.frame(phase1 = paste0("IDV ", ff), batch = rownames(C),
    contaminated = grepl("fallo", rownames(C)), in_raw_subset = raw_in, reweight_weight = mc$mcd.wt, row.names = NULL)
  chk[[as.character(ff)]] <- list(ucl_v7 = ucl_F_adjusted(v7, I = I_LOTE, alpha = ALPHA, m_star = "nominal")$UCL,
    ucl_new = ucl_F_adjusted(new, I = I_LOTE, alpha = ALPHA)$UCL,
    r_v7 = ratio_w(v7$weights), r_new = ratio_w(new$weights),
    dr_v7 = max(d_v7), dr_new = max(d_new), d_simple = d_simple, d_mcd = d_mcd,
    b5 = max(abs(mc$center - new$mu_r)), same_centers = same_centers)
}
tabA <- do.call(rbind, partA); tabB <- do.call(rbind, partB); tabC <- do.call(rbind, partC)

# --- Anchors ---------------------------------------------------------
k7 <- chk[["7"]]; k1 <- chk[["1"]]
okB1 <- abs(k7$ucl_v7 - 19.69285) < 1e-4 && abs(k7$ucl_new - 19.64480) < 1e-4
okB2 <- all(abs(k7$d_simple - c(1.1033, 0.0598, 0.8422, 0.1914)) < 0.005)
okB3 <- same(k7$r_v7, 6.4394, 4) && same(k7$r_new, 4.0226, 4) && same(k1$r_v7, 7.0851, 4) && same(k1$r_new, 4.3089, 4)
okB4 <- same(k7$dr_v7, 0.0614, 4) && same(k1$dr_v7, 0.8167, 4) && same(k7$dr_new, 0.0859, 4) && same(k1$dr_new, 0.0265, 4)
okB5 <- k7$b5 < 1e-10 && k1$b5 < 1e-10
cat(sprintf("B1 UCL V7 %.5f, NEW %.5f -> %s\n", k7$ucl_v7, k7$ucl_new, if (okB1) "OK" else "*** FAIL ***"))
cat(sprintf("B2 IDV 7 shift of the arithmetic mean: %s -> %s\n", paste(sprintf("%.4f", k7$d_simple), collapse = " / "), if (okB2) "OK" else "*** FAIL ***"))
cat(sprintf("B3 ratios IDV 7 %.4f / %.4f, IDV 1 %.4f / %.4f -> %s\n", k7$r_v7, k7$r_new, k1$r_v7, k1$r_new, if (okB3) "OK" else "*** FAIL ***"))
cat(sprintf("B4 drift of mu_r: V7 %.4f / %.4f, NEW %.4f / %.4f -> %s\n", k7$dr_v7, k1$dr_v7, k7$dr_new, k1$dr_new, if (okB4) "OK" else "*** FAIL ***"))
cat(sprintf("B5 MCD of the centers recomputed = mu_r NEW: %.1e / %.1e -> %s\n", k7$b5, k1$b5, if (okB5) "OK" else "*** FAIL ***"))
if (!(okB1 && okB2 && okB3 && okB4 && okB5)) stop("An anchor failed. Nothing is saved. STOP.")
cat(sprintf("   per-batch MCD centers V7 vs NEW, max difference: %.1e / %.1e\n", k7$same_centers, k1$same_centers))
cat(sprintf("C1 IDV 7 shift of the MCD centers: %s (CSV 0.3070 / 0.0926 / 0.2627 / 0.2048), max diff %.4f\n",
            paste(sprintf("%.4f", k7$d_mcd), collapse = " / "), max(abs(k7$d_mcd - c(0.307, 0.0926, 0.2627, 0.2048)))))

cat("\n== A. Shifts by variable (SD of healthy observations) ==\n"); print(tabA, row.names = FALSE, digits = 4)
cat("\n== B. Weights: contaminated batches below 1/30 ==\n")
print(aggregate(cbind(below_uniform_V7, below_uniform_NEW) ~ phase1 + contaminated, data = tabB, FUN = sum))
cat("\n== C. MCD of the 30 centers: contaminated batches excluded (reweight weight 0) ==\n")
print(aggregate(cbind(excluded = reweight_weight == 0, in_raw = in_raw_subset) ~ phase1 + contaminated, data = tabC, FUN = sum))

write.csv(tabA, file.path(DIR_OUT, "table_M4c_center_by_variable.csv"), row.names = FALSE)
write.csv(tabB, file.path(DIR_OUT, "table_M4c_weights.csv"),            row.names = FALSE)
write.csv(tabC, file.path(DIR_OUT, "table_M4c_center_membership.csv"),  row.names = FALSE)
cat("\nAnchors B1-B5 OK\nSaved: table_M4c_center_by_variable.csv, table_M4c_weights.csv, table_M4c_center_membership.csv\n")
