# =====================================================================
# SCRIPT_M4d_TEP_competitors_verify.R
# ---------------------------------------------------------------------
# Revision round 1, P18 (3/3). Reviewer 1, comment 11 (competitors);
# Omar, point 11 ("be honest if the competitor is comparable or better").
#
# WHY
#   The author decided (03-oct-2026) to report the pooled competitors
#   (MRCD and RMCD) in the Tennessee Eastman application, and to narrow
#   one sentence of Section 3.8 ("when whole batches are contaminated,
#   RMCD ... stops detecting"). Every word of the new text must rest on a
#   number. This script verifies the numbers and the mechanism.
#
# WHAT IT DOES (NOTHING NEW IS SIMULATED FOR THE TEP; the simulated
#               limits are read from table_M4_final_A_limits.csv)
#   PART A  Re-fits MRCD and RMCD exactly as SCRIPT_M4_TEP_final_mstar14.R
#           (same batches, set.seed(2026) before each fit) and recounts,
#           with their simulated limits:
#             A1  Phase 1 IDV 7, Phase 2 IDV 7  (Table 10)
#             A2  Phase 1 IDV 7, Phase 2 IDV 1, 13, 8, 3  (Table 11)
#             A3  Phase 1 IDV 1, Phase 2 IDV 7  (Section 4.7)
#           and checks each count against the CSV of M4.
#   PART B  Mechanism. Why RMCD works on the TEP but collapses in block W
#           of the simulation (Section 3.7). Same quantities as
#           SCRIPT_M2M3_diagnostics.R (table_M3_W_diagnostic.csv):
#             - share of the contaminated observations kept by the raw
#               MCD subset of RMCD and by its reweighting step
#             - share of the observations finally used by RMCD that come
#               from contaminated batches
#           plus two descriptors of how different the contaminated
#           batches are from the healthy ones:
#             - ratio of within-batch standard deviation, contaminated /
#               healthy, per variable (simulation, block W: 1 by design)
#             - Mahalanobis distance of each contaminated batch mean to
#               the mean of the healthy batches, with the pooled
#               within-batch covariance of the healthy batches
#               (simulation, block W: shift x sqrt(1' Sigma^-1 1) =
#               1.79 and 3.59 for 1.5 and 3 sd)
#
# ANCHORS (STOP if any fails)
#   N1  MRCD and RMCD reproduce, in A1, A2 and A3, the detections, false
#       alarms, median T2 of the normal batches and lowest T2 of the
#       faulty batches of table_M4_final_B_idv7.csv,
#       table_M4_final_C_faults.csv and table_M4_final_E_center.csv
#   N2  table_M3_W_diagnostic.csv is the one expected (RMCD reweighting
#       keeps 0.942 and 0.655 of the shifted observations)
#
# OUTPUT (05_revision_R1/02_resultados/)
#   table_M4d_competitors_TEP.csv  one row per setting and chart
#   table_M4d_mechanism.csv        Part B, TEP and simulation side by side
#   Run from C:/temp_paper. Expected time: one to two minutes.
# =====================================================================
library(robustT2AFM)
library(robustbase)
library(rrcov)
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

LIM <- read.csv(file.path(DIR_OUT, "table_M4_final_A_limits.csv"), stringsAsFactors = FALSE)
TB  <- read.csv(file.path(DIR_OUT, "table_M4_final_B_idv7.csv"),   stringsAsFactors = FALSE)
TC  <- read.csv(file.path(DIR_OUT, "table_M4_final_C_faults.csv"), stringsAsFactors = FALSE)
TE  <- read.csv(file.path(DIR_OUT, "table_M4_final_E_center.csv"), stringsAsFactors = FALSE)
WD  <- read.csv(file.path(DIR_OUT, "table_M3_W_diagnostic.csv"),   stringsAsFactors = FALSE)

# --- Batches (copied from SCRIPT_M4_TEP_final_mstar14.R) -----------------
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
t2_batches <- function(p2, mu, S) {
  Si <- solve(S); b <- unique(p2$Batch)
  setNames(sapply(b, function(bb) { d <- colMeans(p2[p2$Batch == bb, VARS]) - mu
                                     as.numeric(I_LOTE * t(d) %*% Si %*% d) }), b)
}
resumen <- function(t2, ucl) {
  ef <- grepl("fallo", names(t2))
  c(det = sum(t2[ef] > ucl), fa = sum(t2[!ef] > ucl),
    T2_med_sano = round(median(t2[!ef]), 1), T2_min_fallo = round(min(t2[ef]), 1))
}
ajustar <- function(f1) {
  X <- as.matrix(f1[, VARS])
  set.seed(SEED_CAL); mr  <- suppressWarnings(rrcov::CovMrcd(X, alpha = H_MCD))
  set.seed(SEED_CAL); rm_ <- robustbase::covMcd(X, alpha = H_MCD)
  list(X = X, mr = mr, rm = rm_,
       MRCD        = list(mu = rrcov::getCenter(mr), S = rrcov::getCov(mr)),
       RMCD_pooled = list(mu = rm_$center,           S = rm_$cov))
}
ucl_of <- function(label, chart) {
  u <- LIM$UCL[LIM$phase1 == label & LIM$chart == chart & LIM$limit == "simulated"]
  stopifnot(length(u) == 1); u
}
ref_row <- function(tab, chart, ...) {
  k <- tab$chart == chart & tab$limit == "simulated"
  extra <- list(...)
  for (nm in names(extra)) k <- k & tab[[nm]] == extra[[nm]]
  r <- tab[k, ]; stopifnot(nrow(r) == 1); r
}

# =====================================================================
# PART A: counts with the simulated limits, against the M4 CSVs
# =====================================================================
cat("\n===== PART A: MRCD and RMCD on the Tennessee Eastman =====\n")
fit7 <- ajustar(fase1(7))
fit1 <- ajustar(fase1(1))
settings <- list(
  list(id = "A1", p1 = 7, p2 = 7,  fit = fit7, lab = "IDV 7, step 30", src = "B"),
  list(id = "A2", p1 = 7, p2 = 1,  fit = fit7, lab = "IDV 7, step 30", src = "C"),
  list(id = "A2", p1 = 7, p2 = 13, fit = fit7, lab = "IDV 7, step 30", src = "C"),
  list(id = "A2", p1 = 7, p2 = 8,  fit = fit7, lab = "IDV 7, step 30", src = "C"),
  list(id = "A2", p1 = 7, p2 = 3,  fit = fit7, lab = "IDV 7, step 30", src = "C"),
  list(id = "A3", p1 = 1, p2 = 7,  fit = fit1, lab = "IDV 1, step 30", src = "E"))
filas <- list(); okN1 <- TRUE
for (s in settings) {
  p2 <- fase2(s$p2)
  stopifnot(length(unique(p2$Batch)) == 30)
  for (ch in c("MRCD", "RMCD_pooled")) {
    u  <- ucl_of(s$lab, ch)
    rs <- resumen(t2_batches(p2, s$fit[[ch]]$mu, s$fit[[ch]]$S), u)
    ref <- switch(s$src,
      B = ref_row(TB, ch),
      C = ref_row(TC, ch, IDV = s$p2),
      E = ref_row(TE, ch, phase1 = paste("IDV", s$p1), phase2 = paste("IDV", s$p2)))
    ok <- rs["det"] == ref$det && rs["fa"] == ref$fa &&
          rs["T2_med_sano"] == ref$T2_med_sano && rs["T2_min_fallo"] == ref$T2_min_fallo
    okN1 <- okN1 && ok
    cat(sprintf("%s Phase 1 IDV %-2d Phase 2 IDV %-2d %-12s UCL %.2f | det %2d | FA %d | med normal %.1f | min faulty %.1f | ratio %.2f -> %s\n",
                s$id, s$p1, s$p2, ch, u, rs["det"], rs["fa"], rs["T2_med_sano"], rs["T2_min_fallo"],
                rs["T2_min_fallo"] / u, if (ok) "OK" else "*** DIFFERS FROM M4 CSV ***"))
    filas[[length(filas) + 1]] <- data.frame(setting = s$id, phase1 = paste("IDV", s$p1), phase2 = paste("IDV", s$p2),
      chart = ch, limit = "simulated", UCL = u, t(rs), min_faulty_over_UCL = unname(rs["T2_min_fallo"] / u),
      row.names = NULL)
  }
}
if (!okN1) stop("Anchor N1 failed: MRCD or RMCD do not reproduce the M4 CSVs. STOP.")
cat("N1 OK: every count reproduces the M4 CSVs.\n")
TABA <- do.call(rbind, filas)

# =====================================================================
# PART B: mechanism
# =====================================================================
cat("\n===== PART B: why RMCD works on the TEP =====\n")
okN2 <- nrow(WD) == 2 && all(abs(WD$RMCD_rew_keeps_shifted - c(0.942, 0.655)) < 1e-9)
cat(sprintf("N2 table_M3_W_diagnostic.csv as expected: %s\n", if (okN2) "OK" else "*** FAIL ***"))
if (!okN2) stop("Anchor N2 failed. STOP.")

mecanismo <- function(f1, fit, label) {
  cont <- grepl("fallo", f1$Batch)
  rm_ <- fit$rm
  stopifnot(!is.null(rm_$best), !is.null(rm_$mcd.wt))
  raw_keep <- sum(cont[rm_$best]) / sum(cont)
  rew_keep <- sum(rm_$mcd.wt[cont] > 0) / sum(cont)
  rew_share <- sum(rm_$mcd.wt[cont] > 0) / sum(rm_$mcd.wt > 0)
  hb <- unique(f1$Batch[!cont]); cb <- unique(f1$Batch[cont])
  sd_b <- function(bb) apply(f1[f1$Batch == bb, VARS], 2, sd)
  sd_ratio <- colMeans(do.call(rbind, lapply(cb, sd_b))) / colMeans(do.call(rbind, lapply(hb, sd_b)))
  Sh <- Reduce(`+`, lapply(hb, function(bb) cov(f1[f1$Batch == bb, VARS]))) / length(hb)
  mh <- colMeans(do.call(rbind, lapply(hb, function(bb) colMeans(f1[f1$Batch == bb, VARS]))))
  dM <- sapply(cb, function(bb) { d <- colMeans(f1[f1$Batch == bb, VARS]) - mh; sqrt(as.numeric(t(d) %*% solve(Sh) %*% d)) })
  data.frame(source = label,
             RMCD_raw_keeps_contaminated = raw_keep,
             RMCD_rew_keeps_contaminated = rew_keep,
             RMCD_rew_share_contaminated = rew_share,
             sd_ratio_xmv_9 = sd_ratio[1], sd_ratio_xmv_8 = sd_ratio[2],
             sd_ratio_xmeas_19 = sd_ratio[3], sd_ratio_xmeas_17 = sd_ratio[4],
             sd_ratio_max = max(sd_ratio),
             mahal_contaminated_min = min(dM), mahal_contaminated_mean = mean(dM),
             row.names = NULL)
}
f1_7 <- fase1(7); f1_1 <- fase1(1)
stopifnot(all(f1_7[, VARS] == fit7$X), all(f1_1[, VARS] == fit1$X))
MEC <- rbind(mecanismo(f1_7, fit7, "TEP, Phase 1 with 6 IDV 7 batches"),
             mecanismo(f1_1, fit1, "TEP, Phase 1 with 6 IDV 1 batches"))
q <- sqrt(sum(solve(matrix(c(1, .6, .6, .6, .6, 1, .6, .6, .6, .6, 1, .6, .6, .6, .6, 1), 4))))
SIM <- data.frame(source = sprintf("Simulation, block W, 6 batches shifted %.1f sd (table_M3_W_diagnostic.csv)", WD$shift),
                  RMCD_raw_keeps_contaminated = WD$RMCD_raw_keeps_shifted,
                  RMCD_rew_keeps_contaminated = WD$RMCD_rew_keeps_shifted,
                  RMCD_rew_share_contaminated = WD$RMCD_rew_share_shifted,
                  sd_ratio_xmv_9 = 1, sd_ratio_xmv_8 = 1, sd_ratio_xmeas_19 = 1, sd_ratio_xmeas_17 = 1,
                  sd_ratio_max = 1,
                  mahal_contaminated_min = WD$shift * q, mahal_contaminated_mean = WD$shift * q)
MECH <- rbind(MEC, SIM)
cat("\nShare of the contaminated observations that RMCD keeps, and how different the contaminated batches are:\n")
print(data.frame(source = MECH$source,
                 raw_keeps = round(MECH$RMCD_raw_keeps_contaminated, 3),
                 reweighted_keeps = round(MECH$RMCD_rew_keeps_contaminated, 3),
                 share_in_final_fit = round(MECH$RMCD_rew_share_contaminated, 3),
                 sd_ratio_max = round(MECH$sd_ratio_max, 2),
                 mahal_min = round(MECH$mahal_contaminated_min, 2),
                 mahal_mean = round(MECH$mahal_contaminated_mean, 2)), row.names = FALSE)
cat("\nWithin-batch SD ratio per variable on the TEP (contaminated / healthy):\n")
print(cbind(source = MEC$source, round(MEC[, grep("^sd_ratio_x", names(MEC))], 2)), row.names = FALSE)

write.csv(TABA, file.path(DIR_OUT, "table_M4d_competitors_TEP.csv"), row.names = FALSE)
write.csv(MECH, file.path(DIR_OUT, "table_M4d_mechanism.csv"),       row.names = FALSE)
cat("\nAnchors N1 and N2 OK\nSaved: table_M4d_competitors_TEP.csv, table_M4d_mechanism.csv\n")
