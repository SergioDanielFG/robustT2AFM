# =====================================================================
# SCRIPT_M4_TEP_final_mstar14.R
# ---------------------------------------------------------------------
# Revision round 1: the Tennessee Eastman application with the final
# package (robustT2AFM 0.3.0: MAD scaling, MCD of the batch centers,
# m* = MCD subset size = 14, simulated limit), the published method (v7)
# and the two robust competitors of M2/M3.
# It replaces SCRIPT_M4_TEP_final_method.R for the paper. M4 and its
# tables (table_M4_*.csv) are kept untouched as the record of v7 with
# m* = 13; this script writes table_M4_final_*.csv.
#
# CHARTS AND LIMITS
#   classical    classical Hotelling T2 (package); F limit, exact under
#                normality, and its simulated limit as a check
#   V7           published method: calibrate_afm_mcd(scaling = "none",
#                center = "mean"), Eq. (8) with m* = 13 (m_star =
#                "nominal"), i.e. exactly the published application
#   NEW          final method: calibrate_afm_mcd() with the defaults of
#                0.3.0; Eq. (8) with m* = 14 and ucl_simulated()
#   MRCD         MRCD of all Phase 1 observations pooled (rrcov::CovMrcd,
#                alpha = 0.67), as in M2/M3
#   RMCD_pooled  reweighted MCD of all Phase 1 observations pooled
#                (robustbase::covMcd, alpha = 0.67), as in M2/M3
#   The pooled competitors have no analytic limit for batch means, so
#   each gets a simulated limit built exactly like ucl_simulated(): B
#   synthetic Phase 1 samples from a normal with its own estimates, the
#   estimator refitted on each, and the 1 - alpha quantile of the T2 of
#   new in-control batch means. Every simulated limit uses the same B,
#   n_new and seed, so all charts are brought to the same nominal ARL0.
#   For affine equivariant estimators (classical, RMCD) the simulated
#   limit does not depend on the estimates used to simulate.
#
# PARTS                                             paper table  comments
#   A  limits of every chart (Phase 1 IDV 7 with step 10,
#      20 and 30; Phase 1 IDV 1)                         new      R1-2, R1-3
#   B  Phase 2 with IDV 7                               Table 10   R1-4, R1-6, R1-11
#   C  five faults, Phase 1 fixed (IDV 7)               Table 11   R1-11, R1-12
#   D  subsampling step 10 / 20 / 30                    Table 9    --
#   E  reference center, Phase 1 with IDV 7 or IDV 1    Table 12   R1-6, R1-7, R1-11
#
# DESIGN DECISIONS
#   - Batch construction, runs, starting samples, step, seeds and the V7
#     rows are those of SCRIPT_M4_TEP_final_method.R, copied without
#     changes, so the V7 anchors of M4 must be reproduced.
#   - V7 and NEW come from the package itself (no hand-built NEW as in
#     M4): what the paper reports is what a user of 0.3.0 obtains. Both
#     are calibrated after set.seed(SEED_CAL), so their per-batch MCD fits
#     are the same.
#   - In part E the drift of each center is measured against the same
#     estimator applied to the 24 healthy batches only, in standard
#     deviations of the healthy observations, as in the published table.
#
# ANCHORS AND CHECKS
#   A  (STOP if it fails) V7, Phase 1 IDV 7: UCL 19.69285, weight ratio
#      6.44, covariance inflation 16.14, IDV 7 detected 20/20 with 0/10
#      false alarms; classical 2/20.
#   A2 (STOP if it fails) NEW analytic limit with m* = 14: 19.64480.
#   A3 classical simulated limit within 0.6 of its exact F limit (the
#      simulation reproduces a known answer).
#   C  V7 row IDV 7: 20 (classical 2); shifts 8.33 / 4.89 / 1.78 / 1.76 / 0.30.
#   D  V7 steps 10 / 20 / 30: detected 20 / 20 / 20, false alarms 5 / 0 / 0.
#   E  V7 mean, Phase 1 IDV 1 -> Phase 2 IDV 7: 14 of 20, 6 false alarms,
#      drift 0.82; spatial median: 20 of 20, 0 false alarms, drift 0.06.
#
# OUTPUT (05_revision_R1/02_resultados/)
#   table_M4_final_A_limits.csv   table_M4_final_B_idv7.csv
#   table_M4_final_C_faults.csv   table_M4_final_D_step.csv
#   table_M4_final_E_center.csv
#   Run from C:/temp_paper. The TEP files are not reloaded if the
#   objects ff_test and fy are already in memory. Expected time: 15 to
#   30 minutes (the simulated limits: 4 Phase 1 settings x 4 charts x
#   500 refits).
# =====================================================================
library(robustT2AFM)
library(robustbase)
library(rrcov)
library(MASS)
stopifnot(packageVersion("robustT2AFM") >= "0.3.0",
          exists("ucl_simulated", where = asNamespace("robustT2AFM")))
DIR_DATOS <- "04_datos"
DIR_OUT   <- "05_revision_R1/02_resultados"
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

VARS   <- c("xmv_9", "xmv_8", "xmeas_19", "xmeas_17")
I_LOTE <- 20; PASO <- 30; J <- 4
ALPHA  <- 0.001; H_MCD <- 0.67; ALPHA_CEN <- 0.5
DESDE_SANO <- 1; DESDE_FALLO <- 161
SEED_CAL <- 2026; SEED_CEN <- 2027
SEED_SIM <- 2028; B_SIM <- 500; N_NEW <- 200

if (!exists("ff_test")) { cat("Loading FaultFree_Testing...\n"); ff_test <- read.csv(file.path(DIR_DATOS, "TEP_FaultFree_Testing.csv")) }
if (!exists("fy"))      { cat("Loading Faulty_Testing...\n");    fy      <- read.csv(file.path(DIR_DATOS, "TEP_Faulty_Testing.csv")) }
max_obs <- max(ff_test$sample)
t_start <- Sys.time()

# ---------------------------------------------------------------------
# Helpers (lote, fase1, fase2, spatial_median, mcd_center, resumen,
# pesos_info, cls_t2: copied from M4)
# ---------------------------------------------------------------------
lote <- function(dat, run, fault, desde, nom, paso = PASO) {
  sel <- dat$simulationRun == run & dat$faultNumber == fault &
         dat$sample %in% seq(desde, by = paso, length.out = I_LOTE)
  s <- dat[sel, VARS]
  if (nrow(s) != I_LOTE) return(NULL)
  data.frame(Batch = nom, s, row.names = NULL)
}
fase1 <- function(fault_f1, paso = PASO) do.call(rbind, c(
  lapply(1:24, function(r) lote(ff_test, r, 0,        DESDE_SANO,  sprintf("F1_sano_%02d",  r), paso)),
  lapply(1:6,  function(r) lote(fy,      r, fault_f1, DESDE_FALLO, sprintf("F1_fallo_%02d", r), paso))))
fase2 <- function(fault_f2, paso = PASO) do.call(rbind, c(
  lapply(25:34, function(r) lote(ff_test, r, 0,        DESDE_SANO,  sprintf("F2_sano_%02d",  r), paso)),
  lapply(7:26,  function(r) lote(fy,      r, fault_f2, DESDE_FALLO, sprintf("F2_fallo_%03d", r), paso))))

spatial_median <- function(X, tol = 1e-10, maxit = 2000) {    # as in spatial_median_center.R
  m <- colMeans(X)
  for (it in seq_len(maxit)) {
    d <- sqrt(rowSums((X - matrix(m, nrow(X), ncol(X), byrow = TRUE))^2)); d[d < 1e-12] <- 1e-12
    w <- 1 / d; m_new <- colSums(X * w) / sum(w)
    if (max(abs(m_new - m)) < tol) return(m_new)
    m <- m_new
  }
  warning("Weiszfeld did not converge"); m
}
mcd_center <- function(C) { set.seed(SEED_CEN); covMcd(C, alpha = ALPHA_CEN)$center }

t2_batches <- function(p2, r) {
  Si <- solve(r$Sw); b <- unique(p2$Batch)
  setNames(sapply(b, function(bb) { d <- colMeans(p2[p2$Batch == bb, VARS]) - r$mu
                                     as.numeric(I_LOTE * t(d) %*% Si %*% d) }), b)
}
resumen <- function(t2, ucl) {
  ef <- grepl("fallo", names(t2))
  c(det = sum(t2[ef] > ucl), fa = sum(t2[!ef] > ucl),
    T2_med_sano = round(median(t2[!ef]), 1), T2_min_fallo = round(min(t2[ef]), 1))
}
pesos_info <- function(w) { ic <- grepl("fallo", names(w))
  c(ratio = round(mean(w[!ic]) / mean(w[ic]), 2), en_6_menores = sum(ic[order(w)][1:6])) }
cls_t2 <- function(p2, cls) { m <- hotelling_classical_monitor(p2, cls, VARS); setNames(m$T2, m$Batch)[unique(p2$Batch)] }

# ---------------------------------------------------------------------
# Calibration of every chart on one Phase 1
# ---------------------------------------------------------------------
calibrar <- function(f1) {
  set.seed(SEED_CAL)
  cal7 <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD, scaling = "none", center = "mean"))
  set.seed(SEED_CAL)
  caln <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
  X <- as.matrix(f1[, VARS])
  set.seed(SEED_CAL); mr  <- suppressWarnings(rrcov::CovMrcd(X, alpha = H_MCD))
  set.seed(SEED_CAL); rm_ <- robustbase::covMcd(X, alpha = H_MCD)
  cls <- hotelling_classical_calibrate(f1, VARS)
  list(cal7 = cal7, caln = caln, cls = cls,
       C       = do.call(rbind, cal7$mcd_centers)[, VARS, drop = FALSE],
       ucl_v7  = ucl_F_adjusted(cal7, I = I_LOTE, alpha = ALPHA, m_star = "nominal")$UCL,
       ucl_new = ucl_F_adjusted(caln, I = I_LOTE, alpha = ALPHA)$UCL,
       ucl_cls = hotelling_classical_ucl(K = cls$n_batches, I = I_LOTE, J = J, alpha = ALPHA, phase = "II")$UCL,
       v7   = list(mu = cal7$mu_r, Sw = cal7$Sw, w = cal7$weights),
       new  = list(mu = caln$mu_r, Sw = caln$Sw, w = caln$weights),
       mrcd = list(mu = rrcov::getCenter(mr), Sw = rrcov::getCov(mr)),
       rmcd = list(mu = rm_$center, Sw = rm_$cov),
       clsr = list(mu = cls$mu_global, Sw = cls$Sp))
}

# ---------------------------------------------------------------------
# Simulated limits (same construction as robustT2AFM::ucl_simulated)
# ---------------------------------------------------------------------
fit_cls <- function(f1) {
  b <- unique(f1$Batch)
  list(mu = colMeans(do.call(rbind, lapply(b, function(x) colMeans(f1[f1$Batch == x, VARS])))),
       Sw = Reduce(`+`, lapply(b, function(x) cov(f1[f1$Batch == x, VARS]))) / length(b))
}
fit_mrcd <- function(f1) { m <- suppressWarnings(rrcov::CovMrcd(as.matrix(f1[, VARS]), alpha = H_MCD))
  list(mu = rrcov::getCenter(m), Sw = rrcov::getCov(m)) }
fit_rmcd <- function(f1) { m <- robustbase::covMcd(as.matrix(f1[, VARS]), alpha = H_MCD)
  list(mu = m$center, Sw = m$cov) }
sim_limit <- function(K, fit, r) {
  set.seed(SEED_SIM); t2 <- numeric(0)
  for (b in seq_len(B_SIM)) {
    Xs <- MASS::mvrnorm(K * I_LOTE, r$mu, r$Sw)
    f1 <- data.frame(Batch = rep(sprintf("S%03d", seq_len(K)), each = I_LOTE), Xs); names(f1)[-1] <- VARS
    ft <- fit(f1)
    xb <- MASS::mvrnorm(N_NEW, r$mu, r$Sw / I_LOTE)
    D  <- sweep(xb, 2, ft$mu)
    t2 <- c(t2, I_LOTE * rowSums((D %*% solve(ft$Sw)) * D))
  }
  unname(quantile(t2, 1 - ALPHA, type = 8))
}
limites <- function(Kc, label) {
  K <- length(Kc$caln$weights)
  t0 <- Sys.time()
  u <- c(Kc$ucl_cls,
         sim_limit(K, fit_cls, Kc$clsr),
         Kc$ucl_v7,
         Kc$ucl_new,
         ucl_simulated(Kc$caln, I = I_LOTE, alpha = ALPHA, B = B_SIM, n_new = N_NEW, seed = SEED_SIM)$UCL,
         sim_limit(K, fit_mrcd, Kc$mrcd),
         sim_limit(K, fit_rmcd, Kc$rmcd))
  cat(sprintf("   limits for %-18s %.1f min\n", label, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  data.frame(phase1 = label,
             chart = c("classical", "classical", "V7", "NEW", "NEW", "MRCD", "RMCD_pooled"),
             limit = c("F", "simulated", "Eq8_mstar13_published", "Eq8_mstar14", "simulated", "simulated", "simulated"),
             UCL = round(u, 5), stringsAsFactors = FALSE)
}
# T2 of every chart and its count against each of its limits
evaluar <- function(p2, Kc, L) {
  t2 <- list(classical = cls_t2(p2, Kc$cls), V7 = t2_batches(p2, Kc$v7), NEW = t2_batches(p2, Kc$new),
             MRCD = t2_batches(p2, Kc$mrcd), RMCD_pooled = t2_batches(p2, Kc$rmcd))
  do.call(rbind, lapply(seq_len(nrow(L)), function(i)
    data.frame(chart = L$chart[i], limit = L$limit[i], UCL = round(L$UCL[i], 2),
               t(resumen(t2[[L$chart[i]]], L$UCL[i])), row.names = NULL)))
}

# =====================================================================
# A. Anchors and limits
# =====================================================================
cat("== Calibrating Phase 1 with IDV 7 ==\n")
f1_7 <- fase1(7); f2_7 <- fase2(7)
K7   <- calibrar(f1_7)
rv7  <- resumen(t2_batches(f2_7, K7$v7), K7$ucl_v7); rcl <- resumen(cls_t2(f2_7, K7$cls), K7$ucl_cls)
pv7  <- pesos_info(K7$v7$w); infl <- det(K7$cls$Sp) / det(K7$v7$Sw)
okA  <- abs(K7$ucl_v7 - 19.69285) < 1e-4 && abs(pv7["ratio"] - 6.44) < 0.011 &&
        abs(infl - 16.14) < 0.05 && rv7["det"] == 20 && rv7["fa"] == 0 && rcl["det"] == 2
cat(sprintf("\n== A. ANCHOR V7: UCL %.5f | ratio %.2f | inflation %.2f | IDV7 %d/20, FA %d | classical %d/20  -> %s\n",
            K7$ucl_v7, pv7["ratio"], infl, rv7["det"], rv7["fa"], rcl["det"], if (okA) "OK" else "*** FAIL ***"))
if (!okA) stop("Anchor A failed: V7 does not reproduce the published application. STOP.")
okA2 <- abs(K7$ucl_new - 19.64480) < 1e-4
cat(sprintf("== A2. NEW analytic limit, m* = 14: %.5f (expected 19.64480) -> %s\n",
            K7$ucl_new, if (okA2) "OK" else "*** FAIL ***"))
if (!okA2) stop("Anchor A2 failed: the installed package does not use m* = 14. STOP.")

cat("\n== Simulated limits (the long part) ==\n")
L7 <- limites(K7, "IDV 7, step 30")
dA3 <- abs(L7$UCL[L7$chart == "classical" & L7$limit == "simulated"] - K7$ucl_cls)
okA3 <- dA3 < 0.6
cat(sprintf("== A3. Classical simulated vs exact F limit: difference %.3f -> %s\n",
            dA3, if (okA3) "OK" else "REVIEW"))

# =====================================================================
# B. Phase 2 with IDV 7 (Table 10)
# =====================================================================
tabB <- evaluar(f2_7, K7, L7)
tabB$ratio <- NA; tabB$en_6_menores <- NA
tabB[tabB$chart == "V7",  c("ratio", "en_6_menores")] <- rep(pv7, each = sum(tabB$chart == "V7"))
pnw <- pesos_info(K7$new$w)
tabB[tabB$chart == "NEW", c("ratio", "en_6_menores")] <- matrix(pnw, nrow = sum(tabB$chart == "NEW"), ncol = 2, byrow = TRUE)
cat("\n== B. Phase 2 with IDV 7 (Table 10) ==\n"); print(tabB, row.names = FALSE)

# =====================================================================
# C. Five faults, Phase 1 fixed with IDV 7 (Table 11)
# =====================================================================
desplaz <- function(fault_id) {                   # as in TEP_EXTENDED_VALIDATION_STEP30.R
  ms <- seq(DESDE_FALLO, by = PASO, length.out = I_LOTE)
  normal  <- fy[fy$faultNumber == fault_id & fy$simulationRun %in% 7:26 & fy$sample %in% 1:150, VARS]
  fallado <- fy[fy$faultNumber == fault_id & fy$simulationRun %in% 7:26 & fy$sample %in% ms, VARS]
  max(abs((colMeans(fallado) - colMeans(normal)) / apply(normal, 2, sd)))
}
fallos <- c(1, 13, 7, 8, 3)
tabC <- do.call(rbind, lapply(fallos, function(f)
  data.frame(IDV = f, Shift_SD = round(desplaz(f), 2), evaluar(fase2(f), K7, L7))))
cat("\n== C. Five faults (Table 11) ==\n"); print(tabC, row.names = FALSE)
v7C <- tabC[tabC$chart == "V7", ]; clC <- tabC[tabC$chart == "classical" & tabC$limit == "F", ]
okC <- all(abs(v7C$Shift_SD - c(8.33, 4.89, 1.78, 1.76, 0.30)) < 0.006) &&
       v7C$det[v7C$IDV == 7] == 20 && clC$det[clC$IDV == 7] == 2
cat(sprintf("   Anchor C (V7 and shifts as published): %s\n", if (okC) "OK" else "REVIEW"))

# =====================================================================
# D. Subsampling step (Table 9, detection part)
# =====================================================================
LD <- list(`30` = L7)
tabD <- do.call(rbind, lapply(c(10, 20, 30), function(p) {
  if (DESDE_FALLO + (I_LOTE - 1) * p > max_obs) return(NULL)
  f1 <- fase1(7, p); p2 <- fase2(7, p)
  if (nrow(f1) != 30 * I_LOTE || nrow(p2) != 30 * I_LOTE) return(NULL)
  if (p == 30) { Kp <- K7; Lp <- L7 } else {
    Kp <- calibrar(f1); Lp <- limites(Kp, sprintf("IDV 7, step %d", p)); LD[[as.character(p)]] <<- Lp }
  data.frame(step = p, evaluar(p2, Kp, Lp))
}))
cat("\n== D. Subsampling step (Table 9) ==\n"); print(tabD, row.names = FALSE)
v7D <- tabD[tabD$chart == "V7", ]
okD <- all(v7D$det == c(20, 20, 20)) && all(v7D$fa == c(5, 0, 0))
cat(sprintf("   Anchor D (V7 detected 20/20/20, false alarms 5/0/0): %s\n", if (okD) "OK" else "REVIEW"))

# =====================================================================
# E. Reference center (Table 12): Phase 1 with IDV 7 or IDV 1
# =====================================================================
cat("\n== Calibrating Phase 1 with IDV 1 ==\n")
K1 <- calibrar(fase1(1))
L1 <- limites(K1, "IDV 1, step 30")
sano_ref <- function(fault_f1) {                  # same estimators on the 24 healthy batches only
  f1 <- fase1(fault_f1); fs <- f1[grepl("sano", f1$Batch), ]
  Ks <- calibrar(fs)
  list(sd = apply(fs[, VARS], 2, sd), C = Ks$C, mean = Ks$v7$mu, new = Ks$new$mu,
       mrcd = Ks$mrcd$mu, rmcd = Ks$rmcd$mu, cls = Ks$clsr$mu)
}
Ks   <- list(`7` = K7, `1` = K1)
Ls   <- list(`7` = L7, `1` = L1)
refs <- list(`7` = sano_ref(7), `1` = sano_ref(1))
drift <- function(mu, ref, sd) round(max(abs(mu - ref) / sd), 4)
filas <- list()
for (f1 in c("7", "1")) {
  Kc <- Ks[[f1]]; ref <- refs[[f1]]; L <- Ls[[f1]]
  sm_std <- function(C) spatial_median(sweep(C, 2, ref$sd, "/")) * ref$sd
  # V7 weights with four centers, V7 published limit (rows of M4, anchors)
  cc <- list(mean   = list(mu = Kc$v7$mu,             ref = ref$mean),
             sm     = list(mu = spatial_median(Kc$C), ref = spatial_median(ref$C)),
             sm_std = list(mu = sm_std(Kc$C),         ref = sm_std(ref$C)),
             mcd    = list(mu = mcd_center(Kc$C),     ref = mcd_center(ref$C)))
  for (f2 in c(7, 1)) {
    p2 <- fase2(f2)
    for (est in names(cc)) {
      r <- resumen(t2_batches(p2, list(mu = cc[[est]]$mu, Sw = Kc$v7$Sw)), Kc$ucl_v7)
      filas[[length(filas) + 1]] <- data.frame(
        chart = "V7", center = est, limit = "Eq8_mstar13_published", UCL = round(Kc$ucl_v7, 2),
        phase1 = paste0("IDV ", f1), phase2 = paste0("IDV ", f2),
        drift_SD = drift(cc[[est]]$mu, cc[[est]]$ref, ref$sd), t(r))
    }
    # Final method and competitors, each with its own center and limits
    otros <- list(NEW = list(r = Kc$new, ref = ref$new), MRCD = list(r = Kc$mrcd, ref = ref$mrcd),
                  RMCD_pooled = list(r = Kc$rmcd, ref = ref$rmcd), classical = list(r = Kc$clsr, ref = ref$cls))
    for (ch in names(otros)) {
      t2 <- if (ch == "classical") cls_t2(p2, Kc$cls) else t2_batches(p2, otros[[ch]]$r)
      for (i in which(L$chart == ch)) {
        filas[[length(filas) + 1]] <- data.frame(
          chart = ch, center = "own", limit = L$limit[i], UCL = round(L$UCL[i], 2),
          phase1 = paste0("IDV ", f1), phase2 = paste0("IDV ", f2),
          drift_SD = drift(otros[[ch]]$r$mu, otros[[ch]]$ref, ref$sd), t(resumen(t2, L$UCL[i])))
      }
    }
  }
}
tabE <- do.call(rbind, filas); rownames(tabE) <- NULL
cat("\n== E. Reference center (Table 12) ==\n"); print(tabE, row.names = FALSE)
e1 <- tabE[tabE$chart == "V7" & tabE$center == "mean" & tabE$phase1 == "IDV 1" & tabE$phase2 == "IDV 7", ]
e2 <- tabE[tabE$chart == "V7" & tabE$center == "sm"   & tabE$phase1 == "IDV 1" & tabE$phase2 == "IDV 7", ]
okE <- e1$det == 14 && e1$fa == 6 && abs(e1$drift_SD - 0.82) < 0.01 &&
       e2$det == 20 && e2$fa == 0 && abs(e2$drift_SD - 0.06) < 0.01
cat(sprintf("   Anchor E (mean 14/20, 6 FA, drift 0.82; spatial median 20/20, 0 FA, drift 0.06): %s\n",
            if (okE) "OK" else "REVIEW"))

# =====================================================================
# Output
# =====================================================================
tabA <- do.call(rbind, c(list(L7), LD[setdiff(names(LD), "30")], list(L1))); rownames(tabA) <- NULL
cat("\n== A. Limits of every chart ==\n"); print(tabA, row.names = FALSE)
write.csv(tabA, file.path(DIR_OUT, "table_M4_final_A_limits.csv"), row.names = FALSE)
write.csv(tabB, file.path(DIR_OUT, "table_M4_final_B_idv7.csv"),   row.names = FALSE)
write.csv(tabC, file.path(DIR_OUT, "table_M4_final_C_faults.csv"), row.names = FALSE)
write.csv(tabD, file.path(DIR_OUT, "table_M4_final_D_step.csv"),   row.names = FALSE)
write.csv(tabE, file.path(DIR_OUT, "table_M4_final_E_center.csv"), row.names = FALSE)
cat(sprintf("\nAnchors: A OK | A2 OK | A3 %s | C %s | D %s | E %s   (%.1f min in total)\n",
            if (okA3) "OK" else "REVIEW", if (okC) "OK" else "REVIEW", if (okD) "OK" else "REVIEW",
            if (okE) "OK" else "REVIEW", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("Saved: table_M4_final_A_limits.csv, _B_idv7.csv, _C_faults.csv, _D_step.csv, _E_center.csv\n")
