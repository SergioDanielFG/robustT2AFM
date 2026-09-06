# =====================================================================
# Script_Control_Chart_v2.R
# --------------------------------------------------------------------
#
# CHANGES WITH RESPECT TO Script_Control_Chart.R
#   - Each T2 is divided by its own UCL. The decision rule does not
#     change: a batch is out of control when the ratio exceeds 1, which
#     is the same as the statistic exceeding its limit.
#   - The two panels now share the vertical axis, so the position of a
#     batch relative to its limit can be read across methods, which the
#     previous version required a caveat in the caption to explain.
#   - Adds a PNG export at 600 dpi and exports at 6.5 inches wide, so
#     the figure needs no reduction at 16 cm in the manuscript.
#   Nothing in the computation changes: the values and the limits still
#   come in from the calling script and are not recalculated here.
#
# REQUIRED INPUTS
#   t2_f2_rob, t2_f2_hot, ucl_rob_F, ucl_hot_F, tipo_lote and FIG_OUT,
#   all supplied by BRIDGE_Script_Control_Chart.R.
#
# OUTPUT
#   The PDF named by FIG_OUT and a PNG with the same name at 600 dpi.
# =====================================================================
library(ggplot2)
library(patchwork)

stopifnot(exists("t2_f2_rob"), exists("t2_f2_hot"),
          exists("ucl_rob_F"), exists("ucl_hot_F"), exists("tipo_lote"))

df_ctrl <- data.frame(
  Order = seq_along(t2_f2_rob),
  Type  = factor(tipo_lote, levels = c("Healthy", "Faulty"),
                 labels = c("In-control batch", "Faulty batch")),
  R_rob = t2_f2_rob / ucl_rob_F,
  R_hot = t2_f2_hot / ucl_hot_F
)

y_max <- max(df_ctrl$R_rob, df_ctrl$R_hot) * 1.05

col_map <- c("In-control batch" = "#3FA9B6", "Faulty batch" = "#A02D31")
tema_grafico <- theme_minimal(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")

panel <- function(y, titulo, ucl) {
  ggplot(df_ctrl, aes(x = Order, y = .data[[y]])) +
    geom_hline(yintercept = 1, linetype = "dashed",
               color = "#A02D31", linewidth = 0.5) +
    geom_line(color = "grey70", linewidth = 0.3) +
    geom_point(aes(color = Type), size = 2) +
    scale_color_manual(values = col_map, name = "Batch origin") +
    scale_y_continuous(limits = c(0, y_max)) +
    annotate("text", x = max(df_ctrl$Order), y = 1,
             label = paste0("UCL = ", round(ucl, 2)),
             vjust = -1.2, hjust = 1, size = 2.8, color = "#A02D31") +
    labs(title = titulo, x = "Batch", y = expression(T^2 / UCL)) +
    tema_grafico
}

g_rob <- panel("R_rob", "Proposed AFM-MCD method", ucl_rob_F)
g_hot <- panel("R_hot", "Classical Hotelling method", ucl_hot_F)

combinado <- wrap_plots(g_rob, g_hot, ncol = 2, guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom"))

ggsave(FIG_OUT, combinado, width = 6.5, height = 3.2)
ggsave(sub("\\.pdf$", ".png", FIG_OUT), combinado,
       width = 6.5, height = 3.2, dpi = 600)

cat(sprintf("Robust: %d of %d faulty batches above 1 | min ratio among faulty: %.2f\n",
            sum(df_ctrl$R_rob[df_ctrl$Type == "Faulty batch"] > 1),
            sum(df_ctrl$Type == "Faulty batch"),
            min(df_ctrl$R_rob[df_ctrl$Type == "Faulty batch"])))
cat(sprintf("Classical: %d of %d faulty batches above 1 | max ratio among healthy: %.2f\n",
            sum(df_ctrl$R_hot[df_ctrl$Type == "Faulty batch"] > 1),
            sum(df_ctrl$Type == "Faulty batch"),
            max(df_ctrl$R_hot[df_ctrl$Type == "In-control batch"])))
