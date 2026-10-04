# =====================================================================
# SCRIPT_Fig2_power_final.R
# ---------------------------------------------------------------------
# Revision round 1, P8 (Section 3.3, Figure 2).
#
# WHAT IT DOES
#   Redraws Figure 2 (power curves, TPR against the shift delta) with the
#   final method (robustT2AFM 0.3.0: MAD scaling, MCD of the centers,
#   m* = 14) instead of the published one. Same layout and style as
#   01_scripts/Script_figures_detection_power.R (two panels, SE bands).
#   Nothing is simulated: TPR and SE are read from the power table of M2,
#   where all charts were already equalised to the same ARL0.
#
# INPUT
#   05_revision_R1/02_resultados/table_M2_T4_power.csv
#   columns used: delta, ob, a_classical, e_NEW, SE_a_classical, SE_e_NEW
#
# ANCHORS (the script stops if they fail)
#   A1  delta = 1, 6 contaminated batches: classical 0.2700, NEW 0.8852
#   A2  delta = 1, clean Phase 1:          classical 0.9010, NEW 0.9152
#
# OUTPUT (the old fig_power_curves.* files are NOT overwritten)
#   03_figuras/fig_power_curves_v030.pdf
#   03_figuras/fig_power_curves_v030.png   (600 dpi)
#   Run from C:/temp_paper.
# =====================================================================

dat <- read.csv("05_revision_R1/02_resultados/table_M2_T4_power.csv")

# --- Anchors ---
g <- function(d, o, col) dat[dat$delta == d & dat$ob == o, col]
ok <- abs(g(1, 6, "a_classical") - 0.2700) < 1e-4 && abs(g(1, 6, "e_NEW") - 0.8852) < 1e-4 &&
      abs(g(1, 0, "a_classical") - 0.9010) < 1e-4 && abs(g(1, 0, "e_NEW") - 0.9152) < 1e-4
if (!ok) stop("ANCHOR FAILED: table_M2_T4_power.csv is not the expected one. Figure not drawn.")
cat("Anchors A1-A2 OK\n")

# --- Split the two scenarios ---
clean <- dat[dat$ob == 0, ]; clean <- clean[order(clean$delta), ]
cont  <- dat[dat$ob == 6, ]; cont  <- cont[order(cont$delta), ]

# --- Colors and style (same as the published figure) ---
col_hot <- "#C0392B"   # classical: red
col_rob <- "#185FA5"   # proposed: blue
lwd_l   <- 2.2
cex_pt  <- 1.4

panel <- function(d, title_txt) {
  plot(d$delta, d$a_classical, type = "n",
       xlim = c(0.5, 2.0), ylim = c(0, 1),
       xlab = expression(paste("Mean shift ", delta, " (in ", sigma, ")")),
       ylab = "True positive rate (TPR)",
       main = title_txt, las = 1, cex.axis = 0.95)
  grid(col = "gray85", lty = 1)
  polygon(c(d$delta, rev(d$delta)),
          c(d$a_classical - d$SE_a_classical, rev(d$a_classical + d$SE_a_classical)),
          col = adjustcolor(col_hot, 0.15), border = NA)
  polygon(c(d$delta, rev(d$delta)),
          c(d$e_NEW - d$SE_e_NEW, rev(d$e_NEW + d$SE_e_NEW)),
          col = adjustcolor(col_rob, 0.15), border = NA)
  lines(d$delta, d$a_classical, col = col_hot, lwd = lwd_l)
  points(d$delta, d$a_classical, col = col_hot, pch = 17, cex = cex_pt)
  lines(d$delta, d$e_NEW, col = col_rob, lwd = lwd_l)
  points(d$delta, d$e_NEW, col = col_rob, pch = 16, cex = cex_pt)
  legend("bottomright",
         legend = c(expression(paste("Classical ", T^2)), "MFA-MCD (proposed)"),
         col = c(col_hot, col_rob), lwd = lwd_l, pch = c(17, 16),
         bty = "n", cex = 0.95)
}

draw_fig <- function() {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), mgp = c(2.6, 0.8, 0))
  panel(clean, "No contamination")
  panel(cont,  "Contamination (6 batches)")
}

pdf("03_figuras/fig_power_curves_v030.pdf", width = 10, height = 4.5); draw_fig(); dev.off()
png("03_figuras/fig_power_curves_v030.png", width = 10, height = 4.5, units = "in", res = 600); draw_fig(); dev.off()
draw_fig()
cat("Figure saved to: 03_figuras/fig_power_curves_v030.pdf and 03_figuras/fig_power_curves_v030.png\n")
