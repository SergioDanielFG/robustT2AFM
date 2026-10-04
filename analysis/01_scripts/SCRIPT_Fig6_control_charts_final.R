# =====================================================================
# SCRIPT_Fig6_control_charts_final.R
# ---------------------------------------------------------------------
# Revision round 1, P18 (Section 4.5, Figure 6; Reviewer 2, comment 1
# for the panel title "MFA-MCD").
#
# WHAT IT DOES
#   Redraws Figure 6 (Phase 2 control charts of the Tennessee Eastman
#   application: ten normal batches followed by twenty IDV 7 batches,
#   each T2 divided by its own limit) with the final method
#   (robustT2AFM 0.3.0: MAD scaling of the weights, MCD of the centers,
#   Eq. (8) with m* = 14) instead of the published one (no scaling,
#   mean of the centers, m* = 13, UCL 19.69), and with the panel title
#   "Proposed MFA-MCD method" instead of "Proposed AFM-MCD method".
#   Same layout, colors and size as 01_scripts/Script_Control_Chart_v2.R
#   (the published figure, 03_figuras/fig_control_charts.png).
#
# DATA AND DESIGN (copied without changes from SCRIPT_M4b_TEP_D2_D4.R)
#   Phase 1: 24 normal runs (1-24) + 6 IDV 7 runs (1-6), step 30,
#            I = 20; set.seed(2026) before each calibration.
#   Phase 2: the published composition, normal runs 25-34 and IDV 7
#            runs 7-26, in that order (first ten normal).
#   NOTHING IS SIMULATED.
#
# ANCHORS (the script stops if any fails; values of
#          table_M4_final_B_idv7.csv and table_M4b_D2_D4.csv)
#   A1  published method (V7): UCL 19.69285, lowest faulty T2 43.3
#       (reproduces the published figure)
#   A2  final method (NEW): UCL 19.64480; 20 of 20 detected; 0 of 10
#       false alarms; lowest faulty T2 41.9; median T2 of normal 3.1
#   A3  classical: UCL 19.46440; 2 of 20; 0 of 10; lowest faulty 6.5;
#       median T2 of normal 1.9
#   A4  sentence of Section 4.5: every faulty batch exceeds the limit
#       of the proposed method by a factor of at least 2.1
#
# OUTPUT (the published fig_control_charts.* files are NOT overwritten)
#   03_figuras/fig_control_charts_v030.pdf
#   03_figuras/fig_control_charts_v030.png   (600 dpi)
#   05_revision_R1/02_resultados/table_Fig6_T2_batches.csv
#       (T2 and T2/UCL of each Phase 2 batch, for traceability)
#   Run from C:/temp_paper. Expected time: under one minute.
# =====================================================================
library(robustT2AFM)
library(ggplot2)
library(patchwork)
stopifnot(packageVersion("robustT2AFM") >= "0.3.0")

DIR_DATOS <- "04_datos"
DIR_OUT   <- "05_revision_R1/02_resultados"
VARS   <- c("xmv_9", "xmv_8", "xmeas_19", "xmeas_17")
I_LOTE <- 20; PASO <- 30; J <- 4
ALPHA  <- 0.001; H_MCD <- 0.67
DESDE_SANO <- 1; DESDE_FALLO <- 161
SEED_CAL <- 2026

if (!exists("ff_test")) { cat("Loading FaultFree_Testing...\n"); ff_test <- read.csv(file.path(DIR_DATOS, "TEP_FaultFree_Testing.csv")) }
if (!exists("fy"))      { cat("Loading Faulty_Testing...\n");    fy      <- read.csv(file.path(DIR_DATOS, "TEP_Faulty_Testing.csv")) }

# --- Batches (as in M4 final and M4b) --------------------------------
lote <- function(dat, run, fault, desde, nom) {
  sel <- dat$simulationRun == run & dat$faultNumber == fault &
         dat$sample %in% seq(desde, by = PASO, length.out = I_LOTE)
  s <- dat[sel, VARS]
  if (nrow(s) != I_LOTE) return(NULL)
  data.frame(Batch = nom, s, row.names = NULL)
}
f1 <- do.call(rbind, c(
  lapply(1:24, function(r) lote(ff_test, r, 0, DESDE_SANO,  sprintf("F1_sano_%02d",  r))),
  lapply(1:6,  function(r) lote(fy,      r, 7, DESDE_FALLO, sprintf("F1_fallo_%02d", r)))))
p2 <- do.call(rbind, c(
  lapply(25:34, function(r) lote(ff_test, r, 0, DESDE_SANO,  sprintf("F2_sano_%03d",  r))),
  lapply(7:26,  function(r) lote(fy,      r, 7, DESDE_FALLO, sprintf("F2_fallo_%03d", r)))))
stopifnot(length(unique(f1$Batch)) == 30, length(unique(p2$Batch)) == 30)

# --- Calibrations ----------------------------------------------------
set.seed(SEED_CAL)
K_v7  <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD, scaling = "none", center = "mean"))
set.seed(SEED_CAL)
K_new <- suppressMessages(calibrate_afm_mcd(f1, VARS, mcd_alpha = H_MCD))
K_cls <- hotelling_classical_calibrate(f1, VARS)
ucl_v7  <- ucl_F_adjusted(K_v7,  I = I_LOTE, alpha = ALPHA, m_star = "nominal")$UCL
ucl_new <- ucl_F_adjusted(K_new, I = I_LOTE, alpha = ALPHA)$UCL
ucl_cls <- hotelling_classical_ucl(K = K_cls$n_batches, I = I_LOTE, J = J, alpha = ALPHA, phase = "II")$UCL

b <- unique(p2$Batch)
t2_get <- function(m) setNames(m$T2[match(b, m$Batch)], b)
t2_v7  <- t2_get(monitor_afm_mcd(p2, K_v7,  VARS))
t2_new <- t2_get(monitor_afm_mcd(p2, K_new, VARS))
t2_cls <- t2_get(hotelling_classical_monitor(p2, K_cls, VARS))
ef <- grepl("fallo", b)
stopifnot(sum(ef) == 20, all(!ef[1:10]), all(ef[11:30]))

# --- Anchors ---------------------------------------------------------
r1 <- function(x) round(x, 1)
okA1 <- abs(ucl_v7 - 19.69285) < 1e-4 && r1(min(t2_v7[ef])) == 43.3
okA2 <- abs(ucl_new - 19.64480) < 1e-4 && sum(t2_new[ef] > ucl_new) == 20 &&
        sum(t2_new[!ef] > ucl_new) == 0 && r1(min(t2_new[ef])) == 41.9 &&
        r1(median(t2_new[!ef])) == 3.1
okA3 <- abs(ucl_cls - 19.46440) < 1e-4 && sum(t2_cls[ef] > ucl_cls) == 2 &&
        sum(t2_cls[!ef] > ucl_cls) == 0 && r1(min(t2_cls[ef])) == 6.5 &&
        r1(median(t2_cls[!ef])) == 1.9
fmin <- min(t2_new[ef]) / ucl_new
okA4 <- round(fmin, 1) == 2.1
cat("\n===== ANCHORS =====\n")
cat(sprintf("A1 V7       : UCL %.5f | min faulty T2 %.4f -> %s\n", ucl_v7, min(t2_v7[ef]), if (okA1) "OK" else "*** FAIL ***"))
cat(sprintf("A2 NEW      : UCL %.5f | %d/20 | FA %d | min faulty %.4f | median normal %.4f -> %s\n",
            ucl_new, sum(t2_new[ef] > ucl_new), sum(t2_new[!ef] > ucl_new), min(t2_new[ef]), median(t2_new[!ef]),
            if (okA2) "OK" else "*** FAIL ***"))
cat(sprintf("A3 classical: UCL %.5f | %d/20 | FA %d | min faulty %.4f | median normal %.4f -> %s\n",
            ucl_cls, sum(t2_cls[ef] > ucl_cls), sum(t2_cls[!ef] > ucl_cls), min(t2_cls[ef]), median(t2_cls[!ef]),
            if (okA3) "OK" else "*** FAIL ***"))
cat(sprintf("A4 NEW min faulty T2 / UCL = %.4f (text: at least 2.1) -> %s\n", fmin, if (okA4) "OK" else "*** FAIL ***"))
if (!(okA1 && okA2 && okA3 && okA4)) stop("ANCHOR FAILED: the figure is not drawn. STOP.")

# --- Traceability table ----------------------------------------------
tab <- data.frame(order = seq_along(b), batch = b, type = ifelse(ef, "Faulty", "Healthy"),
                  T2_NEW = t2_new, ratio_NEW = t2_new / ucl_new,
                  T2_classical = t2_cls, ratio_classical = t2_cls / ucl_cls,
                  UCL_NEW = ucl_new, UCL_classical = ucl_cls, row.names = NULL)
write.csv(tab, file.path(DIR_OUT, "table_Fig6_T2_batches.csv"), row.names = FALSE)

# --- Figure (as Script_Control_Chart_v2.R; only the method and title change)
df_ctrl <- data.frame(
  Order = tab$order,
  Type  = factor(tab$type, levels = c("Healthy", "Faulty"),
                 labels = c("In-control batch", "Faulty batch")),
  R_rob = tab$ratio_NEW,
  R_hot = tab$ratio_classical)
y_max <- max(df_ctrl$R_rob, df_ctrl$R_hot) * 1.05
col_map <- c("In-control batch" = "#3FA9B6", "Faulty batch" = "#A02D31")
tema_grafico <- theme_minimal(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal")
panel <- function(y, titulo, ucl) {
  ggplot(df_ctrl, aes(x = Order, y = .data[[y]])) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "#A02D31", linewidth = 0.5) +
    geom_line(color = "grey70", linewidth = 0.3) +
    geom_point(aes(color = Type), size = 2) +
    scale_color_manual(values = col_map, name = "Batch origin") +
    scale_y_continuous(limits = c(0, y_max)) +
    annotate("text", x = max(df_ctrl$Order), y = 1,
             label = sprintf("UCL = %.2f", ucl),
             vjust = -1.2, hjust = 1, size = 2.8, color = "#A02D31") +
    labs(title = titulo, x = "Batch", y = expression(T^2 / UCL)) +
    tema_grafico
}
g_rob <- panel("R_rob", "Proposed MFA-MCD method", ucl_new)
g_hot <- panel("R_hot", "Classical Hotelling method", ucl_cls)
combinado <- wrap_plots(g_rob, g_hot, ncol = 2, guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom"))

dir.create("03_figuras", showWarnings = FALSE)
ggsave("03_figuras/fig_control_charts_v030.pdf", combinado, width = 6.5, height = 3.2)
ggsave("03_figuras/fig_control_charts_v030.png", combinado, width = 6.5, height = 3.2, dpi = 600)

cat(sprintf("\nProposed : %d of 20 faulty above 1 | min ratio faulty %.2f | max ratio normal %.2f\n",
            sum(df_ctrl$R_rob[ef] > 1), min(df_ctrl$R_rob[ef]), max(df_ctrl$R_rob[!ef])))
cat(sprintf("Classical: %d of 20 faulty above 1 | min ratio faulty %.2f | max ratio normal %.2f\n",
            sum(df_ctrl$R_hot[ef] > 1), min(df_ctrl$R_hot[ef]), max(df_ctrl$R_hot[!ef])))
cat("\nSaved: 03_figuras/fig_control_charts_v030.pdf and .png\n")
cat("       05_revision_R1/02_resultados/table_Fig6_T2_batches.csv\n")
