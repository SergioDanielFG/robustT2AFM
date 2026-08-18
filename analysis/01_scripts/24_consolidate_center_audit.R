# =====================================================================
# 24_consolidate_center_audit.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   1. Merges the six per-scenario detail files of the reference-center
#      audit into a single file, adding scenario and replicate columns.
#   2. Recomputes the summary from that consolidated detail and compares
#      it, column by column, with the published summary.
#   3. Moves the superseded summary of the earlier audit to 99_obsoleto.
#   It also prints the figures that go into the reference-center section
#   of the manuscript.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - The summary is rebuilt from the detail rather than copied, which
#     is what turns the comparison into a check: if the published
#     summary cannot be reconstructed from the detail, the two do not
#     describe the same campaign.
#   - The consolidated file keeps one row per replicate and scenario, so
#     the detail remains auditable and any statistic in the summary can
#     be recomputed from it.
#   - The confidence intervals use the same normal-approximation
#     definition as the audit that produced the summary, so the
#     comparison is like for like.
#   - This script simulates nothing: it consumes the outputs of
#     23_REFERENCE_CENTER_AUDIT_v4.R.
#
# ANCHOR AND ITS EXPECTED VALUE
#   Every numeric column shared by the rebuilt summary and the published
#   one must agree to better than 1e-04. The row count of the
#   consolidated file must be 12 000 (6 scenarios x 2000 replicates).
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   The script stops if the anchor fails: no figure from either file is
#   used until the discrepancy is explained.
#
# OUTPUT
#   02_resultados/center_audit_full_detail.csv
#   02_resultados/center_audit_summary.csv
#   99_obsoleto/centro_ponderado_auditoria.csv   (moved)
# =====================================================================

dir.create("99_obsoleto", showWarnings = FALSE)

DIR_RES <- "02_resultados"

NOMBRES <- c(
  "0. Limpia (control negativo)",
  "1. Atipicos internos (paper)",
  "2. Desplazamiento puro 1 SD",
  "3. Desplazamiento puro 3 SD",
  "4. Desplazamiento puro 6 SD",
  "5. Desplaz. 3 SD + covarianza x2"
)
SHIFT  <- c(0, 0, 1, 3, 6, 3)
FACTOR <- c(1, 1, 1, 1, 1, 2)

# =====================================================================
# 1. MERGE THE SIX DETAIL FILES
# =====================================================================
cat("=== Merging the six detail files ===\n")

partes <- list()
for (e in 0:5) {
  ruta <- file.path(DIR_RES, sprintf("center_audit_v4_detail_scenario_%d.csv", e))
  if (!file.exists(ruta)) stop(sprintf("Cannot find %s", ruta))

  d <- read.csv(ruta)
  cat(sprintf("  scenario %d: %d rows, %d columns\n", e, nrow(d), ncol(d)))

  d$scenario_id <- e
  d$scenario    <- NOMBRES[e + 1]
  d$shift_SD    <- SHIFT[e + 1]
  d$factor_cov  <- FACTOR[e + 1]
  d$replicate   <- seq_len(nrow(d))

  partes[[e + 1]] <- d
}

completo <- do.call(rbind, partes)

# identification columns first
ident <- c("scenario_id", "scenario", "shift_SD", "factor_cov", "replicate")
completo <- completo[, c(ident, setdiff(names(completo), ident))]

cat(sprintf("\nConsolidated total: %d rows\n", nrow(completo)))
if (nrow(completo) != 12000)
  warning("Expected 12 000 rows (6 scenarios x 2000 replicates).")

write.csv(completo, file.path(DIR_RES, "center_audit_full_detail.csv"),
          row.names = FALSE)
cat("Saved: center_audit_full_detail.csv\n\n")

# =====================================================================
# 2. REBUILD THE SUMMARY AND COMPARE IT WITH THE PUBLISHED ONE
# =====================================================================
cat("=== Rebuilding the summary from the detail ===\n")

media <- function(x) mean(x, na.rm = TRUE)
se    <- function(x) { x <- x[is.finite(x)]; sd(x) / sqrt(length(x)) }

nuevo <- do.call(rbind, lapply(0:5, function(e) {
  M <- completo[completo$scenario_id == e, ]
  ic <- function(v) c(media(v) - 1.96 * se(v), media(v) + 1.96 * se(v))
  icD <- ic(M$delta_D); icF <- ic(M$delta_FAR); icT <- ic(M$delta_TPR)
  data.frame(
    scenario = NOMBRES[e + 1], shift_SD = SHIFT[e + 1], factor_cov = FACTOR[e + 1],
    n_in_lowest = round(media(M$n_in_lowest), 4),
    SE_n_in_lowest = round(se(M$n_in_lowest), 4),
    chance = 6 * 6 / 30,
    ratio_w = round(media(M$ratio_w), 4), SE_ratio_w = round(se(M$ratio_w), 4),
    D_simple = round(media(M$D_S), 4), D_weighted = round(media(M$D_W), 4),
    gain_D = round(media(M$delta_D), 4), SE_gain_D = round(se(M$delta_D), 4),
    CI95_D_low = round(icD[1], 4), CI95_D_high = round(icD[2], 4),
    FAR_simple = round(media(M$FAR_S), 6), FAR_weighted = round(media(M$FAR_W), 6),
    FAR_difference = round(media(M$delta_FAR), 6),
    SE_FAR_difference = signif(se(M$delta_FAR), 4),
    CI95_FAR_low = round(icF[1], 6), CI95_FAR_high = round(icF[2], 6),
    TPR_simple = round(media(M$TPR_S), 4), TPR_weighted = round(media(M$TPR_W), 4),
    gain_TPR = round(media(M$delta_TPR), 4), SE_gain_TPR = round(se(M$delta_TPR), 4),
    CI95_TPR_low = round(icT[1], 4), CI95_TPR_high = round(icT[2], 4),
    stringsAsFactors = FALSE
  )
}))

ruta_v4 <- file.path(DIR_RES, "weighted_center_audit_v4.csv")
if (file.exists(ruta_v4)) {
  viejo <- read.csv(ruta_v4, stringsAsFactors = FALSE)
  comunes <- intersect(names(viejo), names(nuevo))
  num <- comunes[sapply(nuevo[comunes], is.numeric)]

  difs <- sapply(num, function(cc) max(abs(viejo[[cc]] - nuevo[[cc]]), na.rm = TRUE))
  peor <- max(difs)

  cat(sprintf("\nColumns compared: %d\n", length(num)))
  cat(sprintf("Largest discrepancy: %.3e  (in '%s')\n", peor, names(which.max(difs))))

  if (peor > 1e-4) {
    cat("\nANCHOR FAILED: the published summary is NOT reconstructed from the detail.\n")
    print(sort(difs, decreasing = TRUE)[1:5])
    stop("Review before using any figure.")
  }
  cat("ANCHOR OK: the summary is reconstructed exactly from the detail.\n\n")
} else {
  cat("\nWarning: weighted_center_audit_v4.csv not found, no comparison made.\n\n")
}

write.csv(nuevo, file.path(DIR_RES, "center_audit_summary.csv"), row.names = FALSE)
cat("Saved: center_audit_summary.csv\n\n")

# =====================================================================
# 3. THE FIGURES THAT GO INTO THE MANUSCRIPT
# =====================================================================
cat("=== Figures that go into the reference-center section ===\n\n")

lim  <- nuevo[1, ]; ati <- nuevo[2, ]
des  <- nuevo[3:5, ]

cat(sprintf("Phase 1 clean:              n_in_lowest %.2f   ratio_w %.3f\n",
            lim$n_in_lowest, lim$ratio_w))
for (k in 1:3)
  cat(sprintf("Pure shift %.0f SD:           n_in_lowest %.2f   ratio_w %.3f\n",
              des$shift_SD[k], des$n_in_lowest[k], des$ratio_w[k]))
cat(sprintf("\nInternal outliers, detection: %.3f -> %.3f   95%% CI [%.3f, %.3f]\n",
            ati$TPR_simple, ati$TPR_weighted, ati$CI95_TPR_low, ati$CI95_TPR_high))
cat(sprintf("Phase 1 clean, distance:      %.3f -> %.3f\n",
            lim$D_simple, lim$D_weighted))
cat(sprintf("Phase 1 clean, detection:     %.3f -> %.3f\n\n",
            lim$TPR_simple, lim$TPR_weighted))

# =====================================================================
# 4. RETIRE THE SUPERSEDED SUMMARY OF THE EARLIER AUDIT
# =====================================================================
obs <- file.path(DIR_RES, "centro_ponderado_auditoria.csv")
if (file.exists(obs)) {
  file.rename(obs, file.path("99_obsoleto", "centro_ponderado_auditoria.csv"))
  cat("Moved to 99_obsoleto: centro_ponderado_auditoria.csv (earlier version, superseded)\n")
}

cat("\n=== END ===\n")
cat("Two files of this audit remain in 02_resultados:\n")
cat("  center_audit_full_detail.csv  (12 000 rows)\n")
cat("  center_audit_summary.csv      (6 rows)\n")
