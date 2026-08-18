# =====================================================================
# Script_Control_Chart.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The two Phase 2 control charts of the Tennessee Eastman
#   application, side by side: the proposed AFM-MCD method and the
#   classical Hotelling method, each with its own upper control limit.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - It plots, it does not compute: the T2 values and the limits come
#     in from the calling script and are not recalculated here.
#   - The vertical scales of the two panels differ on purpose; what
#     matters is the position of each batch relative to its own UCL,
#     not the absolute value.
#   - The batch order on the x axis is the one of the incoming vectors:
#     ten healthy batches followed by twenty faulty ones.
#
# REQUIRED INPUTS
#   t2_f2_rob, t2_f2_hot, ucl_rob_F, ucl_hot_F, tipo_lote and FIG_OUT,
#   all of them supplied by BRIDGE_Script_Control_Chart.R. The script
#   is not autonomous and is only run through that bridge.
#
# OUTPUT
#   The PDF file named by FIG_OUT.
# =====================================================================
library(ggplot2)
library(patchwork)

df_ctrl <- data.frame(
  Order = seq_along(t2_f2_rob),
  Type  = factor(tipo_lote, levels = c("Healthy", "Faulty"),
                 labels = c("In-control batch", "Faulty batch")),
  T2_rob = t2_f2_rob,
  T2_hot = t2_f2_hot
)
col_map <- c("In-control batch" = "#3FA9B6", "Faulty batch" = "#A02D31")

tema_grafico <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")

# --- Robust chart ---
g_rob <- ggplot(df_ctrl, aes(x = Order, y = T2_rob)) +
  geom_hline(yintercept = ucl_rob_F, linetype = "dashed", color = "#A02D31", linewidth = 0.5) +
  geom_line(color = "grey70", linewidth = 0.3) +
  geom_point(aes(color = Type), size = 2) +
  scale_color_manual(values = col_map, name = "Batch origin") +
  annotate("text", x = max(df_ctrl$Order), y = ucl_rob_F,
           label = paste0("UCL = ", round(ucl_rob_F, 2)),
           vjust = -0.6, hjust = 1, size = 2.8, color = "#A02D31") +
  labs(title = "Proposed AFM-MCD method", x = "Batch", y = expression(T^2)) +
  tema_grafico

# --- Classical chart ---
g_hot <- ggplot(df_ctrl, aes(x = Order, y = T2_hot)) +
  geom_hline(yintercept = ucl_hot_F, linetype = "dashed", color = "#A02D31", linewidth = 0.5) +
  geom_line(color = "grey70", linewidth = 0.3) +
  geom_point(aes(color = Type), size = 2) +
  scale_color_manual(values = col_map, name = "Batch origin") +
  annotate("text", x = max(df_ctrl$Order), y = ucl_hot_F,
           label = paste0("UCL = ", round(ucl_hot_F, 2)),
           vjust = -0.6, hjust = 1, size = 2.8, color = "#A02D31") +
  labs(title = "Classical Hotelling method", x = "Batch", y = expression(T^2)) +
  tema_grafico

# --- Combine; collect the legend and force it to the bottom ---
combinado <- wrap_plots(g_rob, g_hot, ncol = 2, guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom"))

ggsave(FIG_OUT, combinado, width = 10, height = 4.4)
