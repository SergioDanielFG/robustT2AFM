# ============================================================
# Script_figures_detection_power.R
# ------------------------------------------------------------
# WHAT IT COMPUTES
#   The power-curve figure: the true positive rate against the mean
#   shift delta, for the classical T2 and the proposed method, in two
#   panels, one without contamination and one with six contaminated
#   Phase 1 batches, with shaded bands of plus or minus one standard
#   error.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - It plots, it does not compute: the rates and their standard
#     errors are read from the power table and are not recalculated.
#     Both methods were already equalised to the same ARL0 there, which
#     is what makes the two curves comparable.
#   - Both panels share the vertical scale from 0 to 1, so the gap
#     between methods can be compared across the two scenarios.
#   - The bands are one standard error, not a confidence interval.
#
# ANCHORS
#   None. The figure inherits the anchors of the table it reads.
#
# INPUT
#   02_resultados/table_3_3_power.csv, with columns contam_batches,
#   shift_F2, TPR_Hot, SE_Hot, TPR_Rob, SE_Rob and Advantage.
#
# OUTPUT
#   03_figuras/fig_power_curves.pdf
#   03_figuras/fig_power_curves.png   (600 dpi)
# ============================================================

dir.create("03_figuras", showWarnings = FALSE)

# --- Read results ---
dat <- read.csv("02_resultados/table_3_3_power.csv")

# --- Split the two scenarios ---
clean <- dat[dat$contam_batches == 0, ]
cont  <- dat[dat$contam_batches == 6, ]
clean <- clean[order(clean$shift_F2), ]
cont  <- cont[order(cont$shift_F2), ]

# --- Colors and style (sober, publication-ready) ---
col_hot <- "#C0392B"   # classical: red
col_rob <- "#185FA5"   # robust: blue
lwd_l   <- 2.2
cex_pt  <- 1.4

# --- Function that draws one panel ---
panel <- function(d, title_txt) {
  plot(d$shift_F2, d$TPR_Hot, type = "n",
       xlim = c(0.5, 2.0), ylim = c(0, 1),
       xlab = expression(paste("Mean shift ", delta, " (in ", sigma, ")")),
       ylab = "True positive rate (TPR)",
       main = title_txt, las = 1, cex.axis = 0.95)
  grid(col = "gray85", lty = 1)

  # standard-error band (classical)
  polygon(c(d$shift_F2, rev(d$shift_F2)),
          c(d$TPR_Hot - d$SE_Hot, rev(d$TPR_Hot + d$SE_Hot)),
          col = adjustcolor(col_hot, 0.15), border = NA)
  # standard-error band (robust)
  polygon(c(d$shift_F2, rev(d$shift_F2)),
          c(d$TPR_Rob - d$SE_Rob, rev(d$TPR_Rob + d$SE_Rob)),
          col = adjustcolor(col_rob, 0.15), border = NA)

  # curves
  lines(d$shift_F2, d$TPR_Hot, col = col_hot, lwd = lwd_l)
  points(d$shift_F2, d$TPR_Hot, col = col_hot, pch = 17, cex = cex_pt)
  lines(d$shift_F2, d$TPR_Rob, col = col_rob, lwd = lwd_l)
  points(d$shift_F2, d$TPR_Rob, col = col_rob, pch = 16, cex = cex_pt)

  legend("bottomright",
         legend = c(expression(paste("Classical ", T^2)), "Proposed method"),
         col = c(col_hot, col_rob), lwd = lwd_l, pch = c(17, 16),
         bty = "n", cex = 0.95)
}

# --- Function that draws the full figure (two panels) ---
draw_fig <- function() {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), mgp = c(2.6, 0.8, 0))
  panel(clean, "No contamination")
  panel(cont,  "Contamination (6 batches)")
}

# --- Save vector PDF ---
pdf("03_figuras/fig_power_curves.pdf", width = 10, height = 4.5)
draw_fig()
dev.off()

# --- Save 600 dpi PNG ---
png("03_figuras/fig_power_curves.png", width = 10, height = 4.5, units = "in", res = 600)
draw_fig()
dev.off()

# --- Also show on screen to review ---
draw_fig()

cat("Figure saved to: 03_figuras/fig_power_curves.pdf  and  03_figuras/fig_power_curves.png\n")
