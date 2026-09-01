# =====================================================================
# 25_spatial_median_center.R

# OUTPUT
#   02_resultados/spatial_median_center.csv
# =====================================================================

library(robustT2AFM)

VARS    <- c("xmv_9", "xmv_8", "xmeas_19", "xmeas_17")
I_LOTE  <- 20
PASO    <- 30
J       <- 4
ALPHA   <- 0.001
H_MCD   <- 0.67
DESDE_SANO  <- 1
DESDE_FALLO <- 161
N_SANOS_F1  <- 24
N_CONTA_F1  <- 6
N_SANOS_F2  <- 10
N_FALLO_F2  <- 20
DIR_DATOS <- "04_datos"
DIR_OUT   <- "02_resultados"

muestras_esp <- function(d) seq(d, by = PASO, length.out = I_LOTE)

spatial_median <- function(X, tol = 1e-10, maxit = 2000) {
  m <- colMeans(X)
  for (it in seq_len(maxit)) {
    d <- sqrt(rowSums((X - matrix(m, nrow(X), ncol(X), byrow = TRUE))^2))
    d[d < 1e-12] <- 1e-12
    w <- 1 / d
    m_new <- colSums(X * w) / sum(w)
    if (max(abs(m_new - m)) < tol) return(m_new)
    m <- m_new
  }
  warning("Weiszfeld did not converge")
  m
}

cat("Loading TEP data...\n")
if (!exists("ff")) ff <- read.csv(file.path(DIR_DATOS, "TEP_FaultFree_Testing.csv"))
if (!exists("fy")) fy <- read.csv(file.path(DIR_DATOS, "TEP_Faulty_Testing.csv"))

lote <- function(dat, run, fault, desde, nom) {
  sel <- dat$simulationRun == run & dat$faultNumber == fault &
    dat$sample %in% muestras_esp(desde)
  s <- dat[sel, VARS]
  if (nrow(s) != I_LOTE) return(NULL)
  data.frame(Batch = nom, s, row.names = NULL)
}

calibrar <- function(fallo_f1) {
  f1 <- do.call(rbind, c(
    lapply(1:N_SANOS_F1, function(r)
      lote(ff, r, 0, DESDE_SANO, sprintf("sano_%02d", r))),
    lapply(1:N_CONTA_F1, function(r)
      lote(fy, r, fallo_f1, DESDE_FALLO, sprintf("cont_%02d", r)))
  ))
  stopifnot(nrow(f1) == (N_SANOS_F1 + N_CONTA_F1) * I_LOTE)
  
  cal_r <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
  cal_c <- hotelling_classical_calibrate(f1, VARS)
  K     <- length(cal_r$weights)
  
  cal_s <- suppressMessages(
    calibrate_afm_mcd(f1[grepl("sano", f1$Batch), ], VARS, mcd_alpha = H_MCD))
  sd_ref <- apply(f1[grepl("sano", f1$Batch), VARS], 2, sd)
  
  C  <- do.call(rbind, cal_r$mcd_centers)[, VARS, drop = FALSE]
  Cs <- do.call(rbind, cal_s$mcd_centers)[, VARS, drop = FALSE]
  
  mu_mean <- cal_r$mu_r
  mu_sm   <- spatial_median(C)
  Z       <- sweep(C,  2, sd_ref, "/")
  Zs      <- sweep(Cs, 2, sd_ref, "/")
  mu_smz  <- spatial_median(Z) * sd_ref
  
  ref_mean <- cal_s$mu_r
  ref_sm   <- spatial_median(Cs)
  ref_smz  <- spatial_median(Zs) * sd_ref
  
  nm   <- names(cal_r$weights); cont <- grepl("cont", nm)
  w    <- as.numeric(cal_r$weights)
  
  list(cal_r = cal_r, cal_c = cal_c, K = K,
       ucl_r = ucl_F_adjusted(cal_r, I = I_LOTE, alpha = ALPHA)$UCL,
       ucl_c = hotelling_classical_ucl(K = K, I = I_LOTE, J = J,
                                       alpha = ALPHA, phase = "II")$UCL,
       ratio_w = mean(w[!cont]) / mean(w[cont]),
       infl    = det(cal_c$Sp) / det(cal_r$Sw),
       centers = list(mean = mu_mean, sm = mu_sm, sm_std = mu_smz),
       drift   = c(mean   = max(abs(mu_mean - ref_mean) / sd_ref),
                   sm     = max(abs(mu_sm   - ref_sm)   / sd_ref),
                   sm_std = max(abs(mu_smz  - ref_smz)  / sd_ref)))
}

monitorear <- function(cal, fallo_f2, mu) {
  p2 <- rbind(
    do.call(rbind, lapply(25:(24 + N_SANOS_F2), function(r)
      lote(ff, r, 0, DESDE_SANO, sprintf("F2sano_%02d", r)))),
    do.call(rbind, lapply(7:(6 + N_FALLO_F2), function(r)
      lote(fy, r, fallo_f2, DESDE_FALLO, sprintf("F2fallo_%03d", r))))
  )
  stopifnot(nrow(p2) == (N_SANOS_F2 + N_FALLO_F2) * I_LOTE)
  
  cal_mod <- cal$cal_r
  cal_mod$mu_r <- mu
  mr <- monitor_afm_mcd(p2, cal_mod, VARS)
  bf <- unique(p2$Batch); ef <- grepl("fallo", bf)
  tr <- mr$T2[match(bf, mr$Batch)]
  
  list(det = sum(tr[ef] > cal$ucl_r), fa = sum(tr[!ef] > cal$ucl_r),
       t2_healthy_med = median(tr[!ef]), t2_faulty_min = min(tr[ef]))
}

cals <- list(`7` = calibrar(7), `1` = calibrar(1))

cat("\n=== Anchor: IDV 7 -> IDV 7 with the simple average ===\n")
a <- cals[["7"]]; m <- monitorear(a, 7, a$centers$mean)
cat(sprintf("  UCL %.5f | ratio %.4f | inflation %.4f | %d of 20 | %d FA\n",
            a$ucl_r, a$ratio_w, a$infl, m$det, m$fa))
ok <- abs(a$ucl_r - 19.69285) < 1e-4 && abs(a$ratio_w - 6.4394) < 0.01 &&
  abs(a$infl - 16.1351) < 0.01 && m$det == 20 && m$fa == 0
cat(sprintf("  ANCHOR: %s\n", if (ok) "OK" else "*** FAILS, STOP ***"))
stopifnot(ok)

cat("\n=== Drift of the reference center (SD) ===\n")
for (f in c("7", "1")) {
  d <- cals[[f]]$drift
  cat(sprintf("  Phase 1 with IDV %s: mean %.4f | spatial median %.4f | spatial median (std) %.4f\n",
              f, d["mean"], d["sm"], d["sm_std"]))
}

cat("\n=== CROSS MATRIX, three center estimators ===\n")
filas <- list()
for (est in c("mean", "sm", "sm_std")) {
  for (f1 in c("7", "1")) for (f2 in c("7", "1")) {
    r <- monitorear(cals[[f1]], as.numeric(f2), cals[[f1]]$centers[[est]])
    cat(sprintf("  %-7s | Phase1 IDV %s -> Phase2 IDV %s : %2d of 20 detected, %d false alarms of 10\n",
                est, f1, f2, r$det, r$fa))
    filas[[length(filas) + 1]] <- data.frame(
      center_estimator = est, phase1 = paste0("IDV ", f1),
      phase2 = paste0("IDV ", f2), UCL_rob = cals[[f1]]$ucl_r,
      center_drift_SD = unname(cals[[f1]]$drift[est]),
      det_rob = r$det, fa_rob = r$fa,
      T2_median_healthy = r$t2_healthy_med,
      T2_min_faulty = r$t2_faulty_min)
  }
}
res <- do.call(rbind, filas)

cat("\n=== FULL TABLE ===\n")
print(res, row.names = FALSE)

if (!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)
write.csv(res, file.path(DIR_OUT, "spatial_median_center.csv"), row.names = FALSE)
cat("\nSaved: 02_resultados/spatial_median_center.csv\n")
cat("\n=== END. Copy the whole output. ===\n\n")