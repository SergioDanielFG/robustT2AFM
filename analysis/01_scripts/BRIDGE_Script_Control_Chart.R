# =====================================================================
# BRIDGE_Script_Control_Chart.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   Nothing of its own. It runs the base pipeline, renames its objects
#   to the five names the control-chart script expects, checks the
#   anchors, and then runs that script without modifying it.
#
#   NAME EQUIVALENCE, which is all this file does
#     t2_f2_rob  <- tr        robust T2 of the 30 Phase 2 batches
#     t2_f2_hot  <- tc        classical T2 of the same batches
#     ucl_rob_F  <- ucl_rob   robust UCL
#     ucl_hot_F  <- ucl_cls   classical UCL
#     tipo_lote  <- ef        logical vector -> "Healthy" / "Faulty"
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - Nothing is recomputed: the five objects are renamed, not
#     re-estimated, so the figure shows exactly the numbers the pipeline
#     anchored.
#   - The batch order is the one fixed by b <- unique(f2$Batch) in the
#     pipeline: 10 healthy batches followed by 20 faulty ones. The
#     control-chart script uses seq_along() for the x axis, so that is
#     the order drawn.
#   - The figure is written to a separate working directory so that the
#     published figure in 03_figuras is NOT overwritten. Change FIG_OUT
#     to regenerate the published one.
#
# ANCHORS AND THEIR EXPECTED VALUES
#   robust UCL     = 19.69285   (tolerance 1e-04)
#   classical UCL  = 19.46440   (tolerance 1e-04)
#   detected       = 20 of 20   (exact)
#   false alarms   =  0 of 10   (exact)
#   Structure: 30 batches, 20 faulty and 10 healthy.
#
# SUCCESS CRITERION FIXED IN ADVANCE
#   All the anchors must hold before anything is drawn; the script stops
#   otherwise, since the figure would not correspond to the published
#   application.
#
# OUTPUT
#   99_obsoleto/_verificacion_logs/figura_regenerada/fig_control_charts.pdf
# =====================================================================

DESTINO <- "99_obsoleto/_verificacion_logs/figura_regenerada"
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

# --- 5. Run the original script, without touching it -----------------
source("01_scripts/Script_Control_Chart.R", echo = FALSE)

cat(sprintf("\nFigure written to: %s\n", FIG_OUT))
