# =====================================================================
# Script_Weights_v2.R
# ---------------------------------------------------------------------
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - It plots, it does not compute: the weights come in from the
#     calibration performed by the base pipeline and are not
#     recalculated here.
#   - Batches are sorted by weight and placed on the x axis by position,
#     not by name, and the individual labels are removed: what the
#     figure has to show is where the contaminated batches fall in the
#     ordering, not which run each bar belongs to.
#   - The reference line is 1/K, the weight every batch would get if all
#     contributed equally.
#
#
# OUTPUT
#   03_figuras/fig_pesos_afm_v2.pdf
#   03_figuras/fig_pesos_afm_v2.png   (600 dpi)
# =====================================================================
library(ggplot2)

stopifnot(exists("cal_rob"))

df_pesos <- data.frame(Batch = names(cal_rob$weights),
                       Weight = as.numeric(cal_rob$weights))
df_pesos$Type <- ifelse(grepl("F1_fal", df_pesos$Batch), "Contaminated", "Healthy")
df_pesos <- df_pesos[order(df_pesos$Weight), ]
df_pesos$Order <- seq_len(nrow(df_pesos))   # position, not the name

w_unif <- 1 / nrow(df_pesos)

g_pesos <- ggplot(df_pesos, aes(x = Order, y = Weight, fill = Type)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = w_unif, linetype = "dashed", color = "grey40") +
  annotate("text", x = 1, y = w_unif,
           label = sprintf("Uniform weight = %.3f", w_unif),
           hjust = 0, vjust = -0.8, size = 3, color = "grey40") +
  scale_fill_manual(values = c("Healthy" = "#3FA9B6", "Contaminated" = "#A02D31"),
                    name = "Batch type") +
  scale_y_continuous(breaks = seq(0, 0.125, by = 0.025)) +
  labs(x = "Batch", y = "AFM weight") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_blank(),      # removes the individual labels
        axis.ticks.x = element_blank())

dir.create("03_figuras", showWarnings = FALSE)
ggsave("03_figuras/fig_pesos_afm_v2.pdf", g_pesos, width = 6.5, height = 4)
ggsave("03_figuras/fig_pesos_afm_v2.png", g_pesos, width = 6.5, height = 4, dpi = 600)

cat(sprintf("Uniform weight: %.4f\n", w_unif))
cat(sprintf("Mean weight, healthy batches:      %.5f\n",
            mean(df_pesos$Weight[df_pesos$Type == "Healthy"])))
cat(sprintf("Mean weight, contaminated batches: %.5f\n",
            mean(df_pesos$Weight[df_pesos$Type == "Contaminated"])))
cat(sprintf("Ratio healthy / contaminated:      %.4f\n",
            mean(df_pesos$Weight[df_pesos$Type == "Healthy"]) /
              mean(df_pesos$Weight[df_pesos$Type == "Contaminated"])))
cat("\nFigure saved to 03_figuras/fig_pesos_afm_v2.pdf and .png\n")