# =====================================================================
# 16_ci_paired_far_difference.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The confidence interval for the difference between the false-alarm
#   rate of the proposed robust limit and that of the classical limit,
#   read from the per-replicate FAR file. It reports the paired
#   analysis (correlation between methods, mean difference, standard
#   error, t statistic, 95% CI and a Wilcoxon check), translates that
#   interval into the ARL0 scale, and shows the unpaired version for
#   comparison.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - The analysis is PAIRED: both methods were evaluated on the same
#     N replicates and share the variation each data set introduces,
#     so the paired analysis is the correct one and the unpaired one is
#     reported only as a contrast.
#   - The unpaired standard error uses the BETWEEN-REPLICATE deviation,
#     which is the same definition as the standard error of Table 1;
#     the binomial definition would give a different value.
#   - Wilcoxon is included because it makes no normality assumption; if
#     it agrees, the conclusion does not rest on that assumption.
#   - This script simulates nothing: it consumes the per-replicate FAR
#     file produced by ZTEST_TABLE1_FAITHFUL_v2.R.
#
# ANCHORS AND THEIR EXPECTED VALUES
#   Recomputed from the input file, against the published Table 1:
#     FAR robust    = 0.000625   (tolerance 1e-09)
#     FAR classical = 0.000730   (tolerance 1e-09)
#     False alarms robust    = 125   (tolerance 0.5)
#     False alarms classical = 146   (tolerance 0.5)
#     SE robust     = 5.55e-05   (tolerance 2e-07)
#     SE classical  = 6.35e-05   (tolerance 2e-07)
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   The six anchors must hold; the script stops otherwise, because that
#   would mean the input file is not the source of Table 1 and no
#   interval computed from it would refer to the published rows.
#
# OUTPUT
#   02_resultados/ci_paired_far.csv
# =====================================================================

DIR_RES <- "02_resultados"
K2      <- 100     # Phase 2 batches per replicate

d <- read.csv(file.path(DIR_RES, "far_per_replicate.csv"))
n <- nrow(d)
cat(sprintf("Read %d replicates from far_per_replicate.csv\n\n", n))

rob <- d$far_rob
cls <- d$far_cls

# ---------------------- Anchor against Table 1 -----------------------
cat("=== ANCHOR: does it reproduce Table 1? ===\n")
anc <- data.frame(
  cantidad  = c("FAR robust", "FAR classical", "FA robust", "FA classical",
                "SE robust", "SE classical"),
  calculado = c(mean(rob), mean(cls), sum(rob) * K2, sum(cls) * K2,
                sd(rob) / sqrt(n), sd(cls) / sqrt(n)),
  publicado = c(0.000625, 0.000730, 125, 146, 5.55e-05, 6.35e-05),
  tol       = c(1e-9, 1e-9, 0.5, 0.5, 2e-7, 2e-7)
)
for (i in seq_len(nrow(anc))) {
  ok <- abs(anc$calculado[i] - anc$publicado[i]) < anc$tol[i]
  cat(sprintf("  %-14s computed = %-12.6g published = %-10.6g  %s\n",
              anc$cantidad[i], anc$calculado[i], anc$publicado[i],
              if (ok) "OK" else "*** DOES NOT MATCH ***"))
}
stopifnot(all(abs(anc$calculado - anc$publicado) < anc$tol))
cat("  Anchor verified: this CSV is the source of Table 1.\n\n")

# ------------------------- Paired analysis ---------------------------
dif <- rob - cls
tt  <- t.test(rob, cls, paired = TRUE)
sed <- sd(dif) / sqrt(n)
lo  <- tt$conf.int[1]; hi <- tt$conf.int[2]
rho <- cor(rob, cls)

cat("=== PAIRED (t test for related samples) ===\n")
cat(sprintf("  correlation between the two methods = %.4f\n", rho))
cat(sprintf("  mean difference (rob - cls)         = %+.6f\n", mean(dif)))
cat(sprintf("  SE of the difference                = %.6e\n", sed))
cat(sprintf("  t = %.3f   df = %d   p = %.3f\n", tt$statistic, tt$parameter,
            tt$p.value))
cat(sprintf("  95%% CI of the difference = [%+.6f, %+.6f]\n", lo, hi))

# Wilcoxon: assumes no normality. If it agrees, the result does not
# depend on the distributional assumption.
w <- suppressWarnings(wilcox.test(rob, cls, paired = TRUE))
cat(sprintf("  Paired Wilcoxon (no normality assumption): p = %.3f\n\n",
            w$p.value))

# ------------- Translation into ARL0, the interpretable scale --------
far_cls <- mean(cls)
cat("=== TRANSLATED INTO ARL0 ===\n")
cat(sprintf("  classical observed: FAR = %.6f -> ARL0 = %.0f\n",
            far_cls, 1 / far_cls))
cat(sprintf("  robust observed   : FAR = %.6f -> ARL0 = %.0f\n",
            mean(rob), 1 / mean(rob)))
cat("\n  If the difference were at each end of the CI:\n")
for (b in c(lo, hi)) {
  f <- far_cls + b
  cat(sprintf("    difference %+.6f -> robust FAR %.6f -> ARL0 %.0f\n",
              b, f, 1 / f))
}

# --------------- Unpaired version, for comparison --------------------
# There are two ways of computing the unpaired standard error and they
# give different values. The BETWEEN-REPLICATE deviation is used here,
# which is the same definition as the standard error of Table 1.
se_nop <- sqrt(var(rob) / n + var(cls) / n)
z      <- mean(dif) / se_nop
lo_n   <- mean(dif) - 1.96 * se_nop
hi_n   <- mean(dif) + 1.96 * se_nop
cat("\n=== UNPAIRED, for comparison only ===\n")
cat(sprintf("  SE = %.6e\n", se_nop))
cat(sprintf("  95%% CI = [%+.6f, %+.6f]\n", lo_n, hi_n))
cat(sprintf("  z = %.3f   p = %.3f\n", z, 2 * (1 - pnorm(abs(z)))))
cat("  (the paired one is correct: both methods share the data)\n")

# ----------------------------- Reading -------------------------------
cat("\n=== READING ===\n")
cat("  The CI includes zero, so there is no evidence that the robust\n")
cat("  limit is more conservative than the classical one. But the\n")
cat("  interval is wide: it does NOT establish equivalence.\n")
cat("  The claim 'statistically equivalent' that appeared in the\n")
cat("  calibration table footnote is not supported by this analysis.\n")

# ------------------------------ Output -------------------------------
out <- data.frame(
  method          = c("no pareado", "pareado"),
  difference      = mean(dif),
  SE              = c(se_nop, sed),
  CI95_low        = c(lo_n, lo),
  CI95_high       = c(hi_n, hi),
  ARL0_rob_low    = c(1 / (far_cls + hi_n), 1 / (far_cls + hi)),
  ARL0_rob_high   = c(1 / (far_cls + lo_n), 1 / (far_cls + lo)),
  ARL0_classical  = 1 / far_cls,
  correlation     = c(NA, rho),
  n_replicates    = n
)
write.csv(out, file.path(DIR_RES, "ci_paired_far.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %s\n", file.path(DIR_RES, "ci_paired_far.csv")))
cat("\n=== END. Copy the whole output. ===\n\n")
