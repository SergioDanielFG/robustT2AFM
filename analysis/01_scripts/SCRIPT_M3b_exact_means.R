# =====================================================================
# SCRIPT_M3b_exact_means.R
# ---------------------------------------------------------------------
# Revision round 1, P12 (Section 3.7, Reviewer 1 comment 9, Reviewer 2
# comment 5).
#
# WHY
#   table_M3_all.csv rounds the mean TPR to 4 decimals, and 16 cells end in
#   5, so their 3-decimal value in the paper cannot be decided from the
#   CSV. NOTHING IS SIMULATED: this script re-reads the replicates saved
#   by SCRIPT_M3_adverse_scenarios.R and writes the means with 6 decimals.
#
# ANCHOR (the script stops if it fails)
#   The recomputed means of the 7 charts must reproduce table_M3_all.csv
#   (tolerance 5e-04) in the 29 configurations.
#
# OUTPUT
#   05_revision_R1/02_resultados/table_M3_all_6dec.csv
#   Run from C:/temp_paper. Takes seconds.
# =====================================================================
DIR_OUT <- "05_revision_R1/02_resultados"
DIR_M3  <- "05_revision_R1/03_T2_crudos/M3"
CHARTS  <- c("a_classical", "b_MCD_unif", "c_AFM_cls", "d_V7", "e_NEW", "f_MRCD", "g_RMCD_pooled")
ref <- read.csv(file.path(DIR_OUT, "table_M3_all.csv"))
ok <- TRUE
OUT <- do.call(rbind, lapply(ref$id, function(id) {
  f <- list.files(DIR_M3, paste0("^M3_", id, "_[0-9]{4}\\.rds$"), full.names = TRUE)
  M <- do.call(rbind, lapply(sort(f), readRDS)); stopifnot(nrow(M) == 2000)
  m <- colMeans(M[, CHARTS])
  okk <- all(abs(m - unlist(ref[ref$id == id, CHARTS])) < 5e-4); ok <<- ok && okk
  cat(sprintf("Anchor %-14s %s\n", id, if (okk) "OK" else "*** FAIL ***"))
  data.frame(id = id, t(round(m, 6)), check.names = FALSE)
}))
if (!ok) stop("ANCHOR FAILED: do not use.")
write.csv(OUT, file.path(DIR_OUT, "table_M3_all_6dec.csv"), row.names = FALSE)
print(OUT, row.names = FALSE)
cat("\nALL ANCHORS OK. File written: 05_revision_R1/02_resultados/table_M3_all_6dec.csv\n")
