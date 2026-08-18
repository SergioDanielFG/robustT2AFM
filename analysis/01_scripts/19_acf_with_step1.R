# =====================================================================
# 19_acf_with_step1.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The lag-1 autocorrelation and the Ljung-Box p-value of the four
#   stripper variables on a healthy Tennessee Eastman run, for four
#   subsampling steps: 1, 10, 20 and 30. Step 1 is the un-subsampled
#   series and supplies the original autocorrelation figures quoted in
#   the autocorrelation section of the manuscript, which the previous
#   diagnostic did not cover because it started at step 10.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - A single healthy run (simulationRun 1, faultNumber 0) of the
#     fault-free test set is used, ordered by sample, so that the
#     series is a genuine time series with no run boundaries in it.
#   - The Ljung-Box test uses LB_LAG = 10 lags, the same setting as the
#     diagnostic it extends, so the twelve pre-existing rows must come
#     out identical.
#   - The steps 10, 20 and 30 are recomputed rather than copied, which
#     is what makes the comparison against the earlier file a check and
#     not a restatement.
#
# ANCHOR AND ITS EXPECTED VALUE
#   The twelve rows for steps 10, 20 and 30 must reproduce the earlier
#   Ljung-Box file exactly: maximum absolute difference below 1e-04 in
#   both ACF(1) and the p-value.
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   The script stops if any of those twelve rows moves; step 1 is only
#   trustworthy if the rest of the diagnostic is unchanged. If the
#   earlier file is not present the anchor is skipped with a warning.
#
# OUTPUT
#   02_resultados/step_sensitivity_ljungbox_with_step1.csv
# =====================================================================

VARS    <- c("xmv_9", "xmv_8", "xmeas_19", "xmeas_17")
PASOS   <- c(1, 10, 20, 30)      # step 1 is the addition
LB_LAG  <- 10                    # lags of the Ljung-Box test
DIR_DATOS <- "04_datos"
DIR_OUT   <- "02_resultados"

cat("Loading the healthy reference run...\n")
ff <- read.csv(file.path(DIR_DATOS, "TEP_FaultFree_Testing.csv"))
serie <- ff[ff$simulationRun == 1 & ff$faultNumber == 0, ]
serie <- serie[order(serie$sample), ]
cat(sprintf("Observations in the run: %d\n\n", nrow(serie)))

filas <- list()
for (paso in PASOS) {
  idx <- seq(1, nrow(serie), by = paso)
  sub <- serie[idx, VARS]

  cat(sprintf("--- STEP %d  (%d observations after subsampling) ---\n",
              paso, nrow(sub)))
  for (v in VARS) {
    x  <- sub[[v]]
    a1 <- as.numeric(acf(x, lag.max = 1, plot = FALSE)$acf[2])
    lb <- Box.test(x, lag = LB_LAG, type = "Ljung-Box")
    cat(sprintf("  %-10s ACF(1) = %7.4f   Ljung-Box p = %.4f   %s\n",
                v, a1, lb$p.value,
                if (lb$p.value > 0.05) "independent"
                else "*** AUTOCORRELATED ***"))
    filas[[length(filas) + 1]] <- data.frame(
      step = paso, variable = v, n_obs = nrow(sub),
      acf1 = round(a1, 4), p_ljungbox = round(lb$p.value, 4),
      independent = lb$p.value > 0.05
    )
  }
  cat("\n")
}
tabla <- do.call(rbind, filas)

# ------------------- Anchor against the earlier file ------------------
ruta_ant <- file.path(DIR_OUT, "step_sensitivity_ljungbox.csv")
if (file.exists(ruta_ant)) {
  cat("=== ANCHOR: do the twelve earlier rows match? ===\n")
  ant <- read.csv(ruta_ant)
  nue <- tabla[tabla$step %in% c(10, 20, 30), ]
  ant <- ant[order(ant$step, ant$variable), ]
  nue <- nue[order(nue$step, nue$variable), ]
  d_acf <- max(abs(ant$acf1 - nue$acf1))
  d_p   <- max(abs(ant$p_ljungbox - nue$p_ljungbox))
  cat(sprintf("  maximum difference in ACF(1)  : %.6f\n", d_acf))
  cat(sprintf("  maximum difference in p-value : %.6f\n", d_p))
  cat(sprintf("  ANCHOR: %s\n\n", if (d_acf < 1e-4 && d_p < 1e-4)
      "OK, the twelve rows are identical"
      else "*** THEY DIFFER, review before using step 1 ***"))
  stopifnot(d_acf < 1e-4, d_p < 1e-4)
} else {
  cat("WARNING: the earlier CSV was not found; the anchor could not be checked.\n\n")
}

# --------------------------- What matters ----------------------------
cat("=== THE TWO FIGURES QUOTED IN THE TEXT ===\n")
p1 <- tabla[tabla$step == 1, ]
for (v in c("xmv_9", "xmeas_19")) {
  a <- p1$acf1[p1$variable == v]
  cat(sprintf("  %-10s ACF(1) = %.4f  -> the text says %.2f\n",
              v, a, round(a, 2)))
}

cat("\n=== PROGRESSION AS THE STEP INCREASES ===\n")
for (v in VARS) {
  s <- tabla[tabla$variable == v, ]
  s <- s[order(s$step), ]
  cat(sprintf("  %-10s %s\n", v,
              paste(sprintf("step %2d: %6.3f", s$step, s$acf1),
                    collapse = " | ")))
}

if (!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)
ruta <- file.path(DIR_OUT, "step_sensitivity_ljungbox_with_step1.csv")
write.csv(tabla, ruta, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", ruta))
cat("(the earlier CSV is NOT overwritten; move it to 99_obsoleto)\n")
cat("\n=== END. Copy the whole output. ===\n\n")
