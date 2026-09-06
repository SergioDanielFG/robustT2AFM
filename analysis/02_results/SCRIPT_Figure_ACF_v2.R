# =====================================================================
# SCRIPT_Figure_ACF_v2.R
# ---------------------------------------------------------------------
# Same figure as SCRIPT_Figure_ACF.R with a layout change: the eight
# panels are arranged in four rows by two columns, one row per variable,
# with the original series on the left and the subsampled series on the
# right. The figure is exported at 6.5 inches wide so that it needs no
# reduction when placed at 16 cm in the manuscript, which keeps the axis
# labels legible.
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
  par(mfrow=c(4,2), mar=c(4,4,2.5,1), cex.axis=0.9, cex.lab=0.95)
  for (k in seq_along(VARS)) {
    # left panel: original series
    acf(serie[[VARS[k]]], main="", lag.max=30, ci.col="red",
        xlab="Lag", ylab="ACF")
    title(main=bquote(.(ETIQ[[k]]) ~ "(original)"), cex.main=1)
    # right panel: spaced series (step 30)
    x_sub <- serie[[VARS[k]]][seq(1, nrow(serie), by=PASO)]
    acf(x_sub, main="", lag.max=min(15,length(x_sub)-1), ci.col="red",
        xlab="Lag", ylab="ACF")
    title(main=bquote(.(ETIQ[[k]]) ~ "(spaced, step 30)"), cex.main=1)
  }
}

dir.create("03_figuras", showWarnings=FALSE)
pdf("03_figuras/fig_acf_v2.pdf", width=6.5, height=7); dibujar(); dev.off()
png("03_figuras/fig_acf_v2.png", width=6.5, height=7, units="in", res=600); dibujar(); dev.off()

cat("Ljung-Box values (step 30):\n")
for (v in VARS) {
  x_sub <- serie[[v]][seq(1, nrow(serie), by=PASO)]
  acf1 <- acf(x_sub, plot=FALSE, lag.max=1)$acf[2]
  p_lb <- Box.test(x_sub, lag=min(10,length(x_sub)-1), type="Ljung-Box")$p.value
  cat(sprintf("  %-10s ACF1=%+.3f  Ljung-Box p=%.3f\n", v, acf1, p_lb))
}
cat("\nFigure saved to 03_figuras/fig_acf_v2.pdf and .png\n")