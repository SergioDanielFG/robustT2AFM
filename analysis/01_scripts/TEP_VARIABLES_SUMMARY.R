# =====================================================================
# TEP_VARIABLES_SUMMARY.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The mean, standard deviation and coefficient of variation of the 52
#   Tennessee Eastman process variables under normal operating
#   conditions. Four of these rows appear as Table 6 of the manuscript:
#   the four stripper variables with the highest relative variability.
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - The summary is computed on the fault-free TRAINING set, not the
#     testing one. Variable selection has to characterise normal
#     operation, and the training set is the collection that does so;
#     the testing set is reserved for the batches that enter Phases 1
#     and 2. Computing it on the testing set instead shifts the CV of
#     xmv_9 from 0.056643 to 0.058679.
#   - No fault information is used at any point, so that the selection
#     cannot be biased towards the variables a particular fault
#     disturbs.
#   - No spaced subsampling is applied. The step-30 subsampling exists
#     to restore temporal independence within a batch; a mean and a
#     standard deviation over the whole series do not need it.
#   - The coefficient of variation is sd / |mean|, chosen because the
#     process variables live on very different scales: pressures run in
#     the thousands, flows and valve openings in the tens. Ranking by
#     standard deviation alone would simply rank by scale.
#
# HOW TO RUN
#   The working directory must be the project root, so that the
#   relative paths below resolve:  setwd("C:/temp_paper")
#
# EXPECTED OUTPUT
#   52 rows. The four stripper variables must come out at
#     xmv_9      cv = 0.056643
#     xmv_8      cv = 0.050528
#     xmeas_19   cv = 0.044633
#     xmeas_17   cv = 0.026879
#   which round to the 0.057 / 0.051 / 0.045 / 0.027 published in
#   Table 6.
# =====================================================================

DIR_DATA <- "04_datos"
DIR_OUT  <- "02_resultados"

cat("Working directory:", getwd(), "\n")

if (!exists("ff_train")) {
  cat("Reading TEP_FaultFree_Training.csv...\n")
  ff_train <- read.csv(file.path(DIR_DATA, "TEP_FaultFree_Training.csv"))
}
cat("  rows:", nrow(ff_train), "| columns:", ncol(ff_train), "\n")

vars_all <- c(paste0("xmeas_", 1:41), paste0("xmv_", 1:11))
healthy  <- ff_train[, vars_all]

res <- data.frame(
  variable  = vars_all,
  mean      = sapply(healthy, mean),
  sd        = sapply(healthy, sd),
  cv        = sapply(healthy, function(x) sd(x) / abs(mean(x))),
  row.names = NULL
)

# ---- Anchor check against the published table ----
anchors <- c(xmv_9 = 0.056643, xmv_8 = 0.050528,
             xmeas_19 = 0.044633, xmeas_17 = 0.026879)

cat("\n=== Anchor check, Table 6 ===\n")
ok <- TRUE
for (v in names(anchors)) {
  got <- res$cv[res$variable == v]
  good <- length(got) && abs(got - anchors[v]) < 5e-7
  ok <- ok && good
  cat(sprintf("  %-10s cv = %.6f | expected %.6f | %s\n",
              v, got, anchors[v], if (good) "OK" else "*** DIFFERS ***"))
}

if (ok) {
  write.csv(res, file.path(DIR_OUT, "tep_variables_summary.csv"),
            row.names = FALSE)
  cat("\nAll anchors match. Saved to",
      file.path(DIR_OUT, "tep_variables_summary.csv"), "\n")
} else {
  cat("\nAnchors do not match. Nothing written.\n")
}
