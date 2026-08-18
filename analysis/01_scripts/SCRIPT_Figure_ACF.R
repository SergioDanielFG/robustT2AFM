# =====================================================================
# SCRIPT_Figure_ACF.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The autocorrelation figure of the temporal-independence section:
#   eight panels, one per variable and series, showing the ACF of the
#   four stripper variables on the original series (top row) and on the
#   series after spaced subsampling with step 30 (bottom row). It also
#   prints the lag-1 ACF and the Ljung-Box p-value of the subsampled
#   series.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - A single healthy run (simulationRun 1, faultNumber 0) of the
#     fault-free test set is used, ordered by sample, so the series is a
#     genuine time series with no run boundaries in it.
#   - The two rows share the variables and differ only in the
#     subsampling, which is what makes the comparison read as the effect
#     of the step.
#   - The maximum lag of the bottom row is capped at the length of the
#     subsampled series, which is thirty times shorter.
#   - Labels are in English for consistency with the other figures.
#
# ANCHORS
#   None. The formal check of temporal independence is the Ljung-Box
#   test of SUBSAMPLING_STEP_SENSITIVITY.R and 19_acf_with_step1.R; this
#   script draws the figure and prints those values as a cross-check.
#
# OUTPUT
#   03_figuras/fig_acf.pdf
#   03_figuras/fig_acf.png   (600 dpi)
# =====================================================================
DIR_DATOS <- "04_datos"

VARS  <- c("xmv_9", "xmv_8", "xmeas_19", "xmeas_17")
ETIQ  <- c(expression(xmv[9]), expression(xmv[8]),
           expression(xmeas[19]), expression(xmeas[17]))
PASO  <- 30

if (!exists("ff_test")) ff_test <- read.csv(file.path(DIR_DATOS, "TEP_FaultFree_Testing.csv"))

serie <- ff_test[ff_test$simulationRun==1 & ff_test$faultNumber==0, ]
serie <- serie[order(serie$sample), ]

dibujar <- function() {
  par(mfrow=c(2,4), mar=c(4,4,3,1))
  # Row 1: original series
  for (k in seq_along(VARS)) {
    acf(serie[[VARS[k]]], main="", lag.max=30, ci.col="red",
        xlab="Lag", ylab="ACF")
    title(main=bquote(.(ETIQ[[k]]) ~ "(original)"), cex.main=1)
  }
  # Row 2: spaced series (step 30)
  for (k in seq_along(VARS)) {
    x_sub <- serie[[VARS[k]]][seq(1, nrow(serie), by=PASO)]
    acf(x_sub, main="", lag.max=min(15,length(x_sub)-1), ci.col="red",
        xlab="Lag", ylab="ACF")
    title(main=bquote(.(ETIQ[[k]]) ~ "(spaced, step 30)"), cex.main=1)
  }
}

dir.create("03_figuras", showWarnings=FALSE)
pdf("03_figuras/fig_acf.pdf", width=11, height=5.5); dibujar(); dev.off()
png("03_figuras/fig_acf.png", width=11, height=5.5, units="in", res=600); dibujar(); dev.off()

cat("Ljung-Box values (step 30):\n")
for (v in VARS) {
  x_sub <- serie[[v]][seq(1, nrow(serie), by=PASO)]
  acf1 <- acf(x_sub, plot=FALSE, lag.max=1)$acf[2]
  p_lb <- Box.test(x_sub, lag=min(10,length(x_sub)-1), type="Ljung-Box")$p.value
  cat(sprintf("  %-10s ACF1=%+.3f  Ljung-Box p=%.3f\n", v, acf1, p_lb))
}
cat("\nFigure saved to 03_figuras/fig_acf.pdf and .png\n")
