# =====================================================================
# Script_Weights.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The AFM weight figure of the Tennessee Eastman application: the
#   weight assigned to each of the thirty Phase 1 batches, sorted from
#   lowest to highest, with the six batches coming from faulty runs
#   marked and the uniform weight drawn as a reference line.
#
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
# REQUIRED INPUTS
#   cal_rob, the calibration object left in memory by
#   PHASE1_PIPELINE_STEP30.R. The script is not autonomous.
#
# OUTPUT
#   03_figuras/fig_pesos_afm.pdf
# =====================================================================
library(ggplot2)
df_pesos <- data.frame(Batch = names(cal_rob$weights),
                       Weight = as.numeric(cal_rob$weights))
df_pesos$Type <- ifelse(grepl("F1_fal", df_pesos$Batch), "Contaminated", "Healthy")
df_pesos <- df_pesos[order(df_pesos$Weight), ]
df_pesos$Order <- seq_len(nrow(df_pesos))   # position, not the name
w_unif <- 1 / nrow(df_pesos)

g_pesos <- ggplot(df_pesos, aes(x = Order, y = Weight, fill = Type)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = w_unif, linetype = "dashed", color = "grey40") +
  annotate("text", x = nrow(df_pesos), y = w_unif,
           label = "Uniform weight", hjust = 3.5, vjust = -2.5, size = 3, color = "grey40") +
  scale_fill_manual(values = c("Healthy" = "#3FA9B6", "Contaminated" = "#A02D31"),
                    name = "Batch type") +
  labs(x = "Batch", y = "AFM weight") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_blank(),      # removes the individual labels
        axis.ticks.x = element_blank())

ggsave("03_figuras/fig_pesos_afm.pdf", g_pesos, width = 8, height = 4.5)
