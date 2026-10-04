# =====================================================================
# SCRIPT_M2c_directions_factorial.R
# ---------------------------------------------------------------------
# Revision round 1, P11 (Section 3.6, Reviewer 1 comments 5 and 10).
#
# WHY
#   table_M2_directions.csv reports charts a, d, e, f, g by direction, but
#   M2 also computed b (MCD only) and c (MFA weighting only) in every
#   replicate. They show which component produces the dependence on the
#   direction of the shift. This script tabulates them. NOTHING IS
#   SIMULATED: it only reads the replicates saved by SCRIPT_M2_power.R.
#
# WHAT IT DOES
#   For P_d1.0 (joint increase) and the 10 configurations of block D, with
#   a clean and a contaminated Phase 1: mean TPR (6 decimals) and SE of
#   the 7 charts, and the paired difference of b, c, d, e against the
#   classical chart (same |z| > 3 reading rule as SCRIPT_M2M3_diagnostics).
#
# ANCHOR (the script stops if it fails)
#   The means of a, d, e, f, g must reproduce table_M2_directions.csv
#   (tolerance 5e-04) in the 12 rows.
#
# OUTPUT (05_revision_R1/02_resultados/)
#   table_M2c_directions_factorial.csv   means and SE, 7 charts
#   table_M2c_directions_paired.csv      b, c, d, e minus classical, paired
#   Run from C:/temp_paper. Takes seconds.
# =====================================================================
DIR_OUT <- "05_revision_R1/02_resultados"
DIR_M2  <- "05_revision_R1/03_T2_crudos/M2"
CHARTS  <- c("a_classical", "b_MCD_unif", "c_AFM_cls", "d_V7", "e_NEW", "f_MRCD", "g_RMCD_pooled")
cfg <- data.frame(id = c("P_d1.0_c0", "P_d1.0_c6", paste0("D_", rep(c("mD1", "D4", "mD4", "D2", "D3"), each = 2), c("_c0", "_c6"))),
                  dir = c("D1", "D1", rep(c("mD1", "D4", "mD4", "D2", "D3"), each = 2)),
                  ob = rep(c(0, 6), 6), stringsAsFactors = FALSE)
read_m2 <- function(id) {
  f <- list.files(DIR_M2, paste0("^M2_", id, "_[0-9]{4}\\.rds$"), full.names = TRUE)
  M <- do.call(rbind, lapply(sort(f), readRDS)); stopifnot(nrow(M) == 2000); M
}
MM <- setNames(lapply(cfg$id, read_m2), cfg$id)

# --- Anchor ---
ref <- read.csv(file.path(DIR_OUT, "table_M2_directions.csv"))
ok <- TRUE
for (k in seq_len(nrow(cfg))) {
  r <- ref[ref$dir == cfg$dir[k] & ref$ob == cfg$ob[k], ]
  got <- colMeans(MM[[k]][, c("a_classical", "d_V7", "e_NEW", "f_MRCD", "g_RMCD_pooled")])
  exp <- unlist(r[, c("a_classical", "d_V7", "e_NEW", "f_MRCD", "g_RMCD_pooled")])
  okk <- all(abs(got - exp) < 5e-4); ok <- ok && okk
  cat(sprintf("Anchor %-10s %s\n", cfg$id[k], if (okk) "OK" else "*** FAIL ***"))
}
if (!ok) stop("ANCHOR FAILED: the M2 replicates do not reproduce table_M2_directions.csv. Do not use.")

# --- Tables ---
FACT <- do.call(rbind, lapply(seq_len(nrow(cfg)), function(k) { M <- MM[[k]]
  data.frame(cfg[k, c("dir", "ob")], t(round(colMeans(M[, CHARTS]), 6)),
             t(setNames(round(apply(M[, CHARTS], 2, function(x) sd(x) / sqrt(length(x))), 4), paste0("SE_", CHARTS))),
             row.names = NULL) }))
PAIR <- do.call(rbind, lapply(seq_len(nrow(cfg)), function(k) { M <- MM[[k]]
  do.call(rbind, lapply(c("b_MCD_unif", "c_AFM_cls", "d_V7", "e_NEW"), function(o) {
    dd <- M[, o] - M[, "a_classical"]; se <- sd(dd) / sqrt(length(dd)); z <- mean(dd) / se
    data.frame(dir = cfg$dir[k], ob = cfg$ob[k], chart = o, versus = "a_classical",
               diff = round(mean(dd), 6), SE_paired = round(se, 4), z = round(z, 1),
               verdict = if (abs(z) <= 3) "no difference" else if (z > 0) "better" else "worse")
  })) }))
write.csv(FACT, file.path(DIR_OUT, "table_M2c_directions_factorial.csv"), row.names = FALSE)
write.csv(PAIR, file.path(DIR_OUT, "table_M2c_directions_paired.csv"),    row.names = FALSE)
print(FACT[, c("dir", "ob", "a_classical", "b_MCD_unif", "c_AFM_cls", "d_V7", "e_NEW")], row.names = FALSE)
cat("\nALL ANCHORS OK. Files written in 05_revision_R1/02_resultados\n")
