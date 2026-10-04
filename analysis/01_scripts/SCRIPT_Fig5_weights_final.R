# =====================================================================
# SCRIPT_Fig5_weights_final.R
# ---------------------------------------------------------------------
# Revision round 1, P17 (Section 4.4, Figure 5; Reviewer 2, comment 1
# for the axis label).
#
# WHAT IT DOES
#   Redraws Figure 5 (MFA weight of each Phase 1 batch of the Tennessee
#   Eastman application, ordered from lowest to highest) with the final
#   method (robustT2AFM 0.3.0: MAD scaling of the weights, MCD of the
#   centers, m* = 14) instead of the published one, and with the axis
#   label "MFA weight" instead of "AFM weight". Same layout, colors and
#   size as 01_scripts/Script_Weights_v2.R (the published figure,
#   03_figuras/fig_pesos_afm_v2.png). Nothing is computed: the weights
#   are read from the CSV written by SCRIPT_M4c_TEP_center_weights.R.
#
# INPUT
#   05_revision_R1/02_resultados/table_M4c_weights.csv
#   rows phase1 == "IDV 7"; columns batch, contaminated, w_NEW
#
# ANCHORS (the script stops if they fail)
#   A1  30 batches, 6 contaminated, weights add up to 1
#   A2  ratio healthy / contaminated mean weight = 4.0226
#       (table_M4b_D2_D4.csv and table_M4_final_B_idv7.csv)
#   A3  the 6 contaminated batches are the 6 lowest weighted and all lie
#       below the uniform weight 1/30 (sentences of Section 4.4 and the
#       caption of Figure 5)
#
# OUTPUT (the published fig_pesos_afm_v2.* files are NOT overwritten)
#   03_figuras/fig_pesos_mfa_v030.pdf
#   03_figuras/fig_pesos_mfa_v030.png   (600 dpi, 3900 x 2400 px, as the
#                                        published one)
#   Run from C:/temp_paper.
# =====================================================================
library(ggplot2)

w <- read.csv("05_revision_R1/02_resultados/table_M4c_weights.csv", stringsAsFactors = FALSE)
w <- w[w$phase1 == "IDV 7", ]

# --- Anchors ---
cont  <- as.logical(w$contaminated)
ratio <- mean(w$w_NEW[!cont]) / mean(w$w_NEW[cont])
w_unif <- 1 / nrow(w)
okA1 <- nrow(w) == 30 && sum(cont) == 6 && abs(sum(w$w_NEW) - 1) < 1e-12
okA2 <- abs(round(ratio, 4) - 4.0226) < 1e-9
okA3 <- all(cont[order(w$w_NEW)][1:6]) && all(w$w_NEW[cont] < w_unif)
cat(sprintf("A1 30 batches, 6 contaminated, sum 1: %s\n", if (okA1) "OK" else "*** FAIL ***"))
cat(sprintf("A2 ratio healthy/contaminated %.4f (expected 4.0226): %s\n", ratio, if (okA2) "OK" else "*** FAIL ***"))
cat(sprintf("A3 6 lowest are the contaminated, all below 1/30: %s\n", if (okA3) "OK" else "*** FAIL ***"))
if (!(okA1 && okA2 && okA3)) stop("ANCHOR FAILED: table_M4c_weights.csv is not the expected one. Figure not drawn.")

# --- Data for the plot (as Script_Weights_v2.R) ---
df_pesos <- data.frame(Batch = w$batch, Weight = w$w_NEW,
                       Type = ifelse(cont, "Contaminated", "Healthy"))
df_pesos <- df_pesos[order(df_pesos$Weight), ]
df_pesos$Order <- seq_len(nrow(df_pesos))   # position, not the name

g_pesos <- ggplot(df_pesos, aes(x = Order, y = Weight, fill = Type)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = w_unif, linetype = "dashed", color = "grey40") +
  annotate("text", x = 1, y = w_unif,
           label = sprintf("Uniform weight = %.3f", w_unif),
           hjust = 0, vjust = -0.8, size = 3, color = "grey40") +
  scale_fill_manual(values = c("Healthy" = "#3FA9B6", "Contaminated" = "#A02D31"),
                    name = "Batch type") +
  scale_y_continuous(breaks = seq(0, 0.08, by = 0.02)) +
  labs(x = "Batch", y = "MFA weight") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

dir.create("03_figuras", showWarnings = FALSE)
ggsave("03_figuras/fig_pesos_mfa_v030.pdf", g_pesos, width = 6.5, height = 4)
ggsave("03_figuras/fig_pesos_mfa_v030.png", g_pesos, width = 6.5, height = 4, dpi = 600)

cat(sprintf("\nUniform weight: %.4f\n", w_unif))
cat(sprintf("Mean weight, healthy batches:      %.5f\n", mean(df_pesos$Weight[df_pesos$Type == "Healthy"])))
cat(sprintf("Mean weight, contaminated batches: %.5f\n", mean(df_pesos$Weight[df_pesos$Type == "Contaminated"])))
cat(sprintf("Ratio healthy / contaminated:      %.4f\n", ratio))
cat(sprintf("Largest weight: %.4f\n", max(df_pesos$Weight)))
cat("\nFigure saved to 03_figuras/fig_pesos_mfa_v030.pdf and .png\n")
