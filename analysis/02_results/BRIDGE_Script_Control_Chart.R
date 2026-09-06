# =====================================================================
# BRIDGE_Script_Control_Chart.R
# --------------------------------------------------------------------
# Runs the Phase 1 pipeline (PHASE1_PIPELINE_STEP30.R), verifies the
# numerical anchors of the manuscript, and hands the statistics and
# limits to Script_Control_Chart_v2.R, which draws Figure 5.
# Nothing is recomputed here.
#
# OUTPUT
#   03_figuras/fig_control_charts.pdf and .png (600 dpi)
# =====================================================================


DESTINO <- "03_figuras"
dir.create(DESTINO, recursive = TRUE, showWarnings = FALSE)
FIG_OUT <- file.path(DESTINO, "fig_control_charts.pdf")

# --- 1. Base pipeline: computes tr, tc, ucl_rob, ucl_cls, ef ---------
source("01_scripts/PHASE1_PIPELINE_STEP30.R", echo = FALSE)

# --- 2. Check that the pipeline left what is expected ----------------
stopifnot(exists("tr"), exists("tc"), exists("ucl_rob"),
          exists("ucl_cls"), exists("ef"))
stopifnot(length(tr) == 30, length(tc) == 30, length(ef) == 30)
stopifnot(sum(ef) == 20, sum(!ef) == 10)

# --- 3. Anchors, before drawing anything -----------------------------
cat("\n===== BRIDGE ANCHORS =====\n")
cat(sprintf("  Robust UCL    : %.5f  (expected 19.69285)  %s\n",
            ucl_rob, if (abs(ucl_rob - 19.69285) < 1e-4) "OK" else "*** FAILS ***"))
cat(sprintf("  Classical UCL : %.5f  (expected 19.46440)  %s\n",
            ucl_cls, if (abs(ucl_cls - 19.46440) < 1e-4) "OK" else "*** FAILS ***"))
cat(sprintf("  Detected      : %d of 20 (expected 20)  %s\n",
            sum(tr[ef] > ucl_rob), if (sum(tr[ef] > ucl_rob) == 20) "OK" else "*** FAILS ***"))
cat(sprintf("  False alarms  : %d of 10 (expected 0)  %s\n",
            sum(tr[!ef] > ucl_rob), if (sum(tr[!ef] > ucl_rob) == 0) "OK" else "*** FAILS ***"))
stopifnot(abs(ucl_rob - 19.69285) < 1e-4,
          abs(ucl_cls - 19.46440) < 1e-4,
          sum(tr[ef] > ucl_rob) == 20,
          sum(tr[!ef] > ucl_rob) == 0)

# --- 4. THE BRIDGE: rename, recompute nothing ------------------------
t2_f2_rob <- tr
t2_f2_hot <- tc
ucl_rob_F <- ucl_rob
ucl_hot_F <- ucl_cls
tipo_lote <- ifelse(ef, "Faulty", "Healthy")

cat("\n  Batch order (the first 10 must be Healthy):\n    ")
cat(paste(substr(tipo_lote, 1, 1), collapse = ""), "\n")

# --- 5. Draw Figure 5 with the v2 script (T2 divided by UCL) --------
source("01_scripts/Script_Control_Chart_v2.R", echo = FALSE)

cat(sprintf("\nFigure written to: %s\n", FIG_OUT))
