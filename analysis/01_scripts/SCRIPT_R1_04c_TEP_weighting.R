# =====================================================================
# SCRIPT_R1_04c_TEP_weighting.R
# ---------------------------------------------------------------------
# Revision round 1, Reviewer 1, comment 4 (and Reviewer 3, note 1):
# scale dependence of the AFM weighting. Real-data part (Tennessee
# Eastman, step 30, IDV 7), the same Phase 1 / Phase 2 as the paper.
#
# WHAT IT COMPUTES
#   The four TEP variables are in their own units (xmv_9, xmv_8 are valve
#   openings in %, xmeas_19 and xmeas_17 are process measurements), so
#   this is the case where the scale question matters in practice.
#   Three weighting schemes, all on the SAME MCD fits and the SAME
#   reference center, so only the weights differ:
#     (d) published: weights from lambda1 of each MCD covariance
#     (e) weights from lambda1 of each MCD correlation matrix
#     (f) global robust standardisation: each variable divided by
#         s_j = median over batches of the MCD standard deviation,
#         weights from lambda1 of D^-1 S_k D^-1
#   For each scheme it reports: healthy / contaminated weight ratio,
#   how many of the 6 IDV-7 batches are among the 6 lowest weights,
#   faulty batches detected out of 20 and false alarms out of 10 with
#   the Eq. (8) limit (which does not involve the weights).
#   It then repeats everything with each variable multiplied by 100
#   (a change of units), to show which scheme changes and which not.
#
# DESIGN DECISIONS
#   - Phase 1 / Phase 2 construction is PHASE1_PIPELINE_STEP30.R, sourced
#     unchanged (it prints its own seven anchors).
#   - covMcd() draws random subsets. To compare the original units with
#     the rescaled ones on equal terms, every calibration in this script
#     is preceded by set.seed(SEED_TEP): with the same draws, the MCD is
#     exactly equivariant, so (e) and (f) must not change at all.
#
# ANCHORS
#   C1  (d) rebuilt here equals the package: Sw and T2 identical to
#       calibrate_afm_mcd() / monitor_afm_mcd() (max abs diff < 1e-8).
#   C2  (e) and (f): weights identical (< 1e-10) under every rescaling.
#
# OUTPUT
#   05_revision_R1/02_resultados/table_R1_04c_TEP_weighting.csv
#   Run from C:/temp_paper. Single core; loading the TEP CSV files is the
#   slow part, the calculations themselves take seconds.
# =====================================================================
library(robustT2AFM)
dir.create("05_revision_R1/02_resultados", recursive = TRUE, showWarnings = FALSE)

source("01_scripts/PHASE1_PIPELINE_STEP30.R")   # f1, f2, VARS, I_LOTE, ALPHA, H_MCD, b, ef

SEED_TEP <- 2026
lam1 <- function(S) max(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
comb <- function(Sl, w) Reduce(`+`, Map(function(S, wi) wi * S, Sl, w))
t2_b <- function(X, mu, Si){ d <- colMeans(X) - mu; as.numeric(I_LOTE * t(d) %*% Si %*% d) }

esquemas <- function(cal){
  S  <- cal$mcd_covariances
  wd <- cal$weights[names(S)]
  le <- sapply(S, function(s) lam1(cov2cor(s)));          we <- (1 / le) / sum(1 / le)
  sg <- apply(sapply(S, function(s) sqrt(diag(s))), 1, median); Di <- diag(1 / sg)
  lf <- sapply(S, function(s) lam1(Di %*% s %*% Di));     wf <- (1 / lf) / sum(1 / lf)
  list(w = list(d = wd, e = we, f = wf), s_glob = sg)
}

evaluar <- function(f1x, f2x, etiqueta){
  set.seed(SEED_TEP)
  cal <- calibrate_afm_mcd(f1x, VARS, mcd_alpha = H_MCD)
  ucl <- ucl_F_adjusted(cal, I = I_LOTE, alpha = ALPHA)$UCL
  E   <- esquemas(cal); S <- cal$mcd_covariances
  filas <- lapply(names(E$w), function(m){
    w  <- E$w[[m]]
    Sw <- comb(S, w)
    T2 <- sapply(b, function(bb) t2_b(as.matrix(f2x[f2x$Batch == bb, VARS]), cal$mu_r, solve(Sw)))
    ic <- grepl("fallo", names(w))
    if (m == "d") {                                   # anchor C1
      mon <- monitor_afm_mcd(f2x, cal, VARS)
      c1  <- max(abs(Sw - cal$Sw), abs(T2 - mon$T2[match(b, mon$Batch)]))
      if (c1 > 1e-8) warning(sprintf("C1 FAIL (%s): max abs diff %.2e", etiqueta, c1))
    }
    data.frame(unidades = etiqueta, pesos = m,
               ratio_sano_cont   = round(mean(w[!ic]) / mean(w[ic]), 2),
               fallo_en_6_menores = sum(ic[order(w)][1:6]),
               deteccion_de_20   = sum(T2[ef]  > ucl),
               falsas_alarmas_de_10 = sum(T2[!ef] > ucl),
               T2_min_fallo = round(min(T2[ef]), 1), T2_max_sano = round(max(T2[!ef]), 1),
               UCL = round(ucl, 2))
  })
  list(tabla = do.call(rbind, filas), w = E$w, s_glob = E$s_glob)
}

# --- Original units + each variable multiplied by 100 ------------------
res0 <- evaluar(f1, f2, "original")
cat("\nRobust scale of each variable in the TEP Phase 1 (s_j):\n"); print(round(res0$s_glob, 3))

todas <- list(res0$tabla); ok_c2 <- TRUE
for (v in VARS) {
  f1x <- f1; f2x <- f2; f1x[[v]] <- 100 * f1x[[v]]; f2x[[v]] <- 100 * f2x[[v]]
  r <- evaluar(f1x, f2x, paste0(v, " x100"))
  todas[[length(todas) + 1]] <- r$tabla
  for (m in c("e", "f")) {
    dmax <- max(abs(r$w[[m]] - res0$w[[m]]))
    ok_c2 <- ok_c2 && dmax < 1e-10
    cat(sprintf("C2  %-12s scheme (%s): max weight change %.1e %s\n",
                paste0(v, " x100"), m, dmax, ifelse(dmax < 1e-10, "OK", "FAIL")))
  }
}
tab <- do.call(rbind, todas)

cat("\n== TEP, IDV 7: weighting scheme x units ==\n")
print(tab, row.names = FALSE)
write.csv(tab, "05_revision_R1/02_resultados/table_R1_04c_TEP_weighting.csv", row.names = FALSE)
cat(if (ok_c2) "\nC2 OK: (e) and (f) do not depend on the units.\n" else "\nC2 FAILURES: check before using.\n")
cat("C1: a warning above means the rebuilt (d) does not match the package.\n")
cat("Saved to 05_revision_R1/02_resultados/table_R1_04c_TEP_weighting.csv\n")
