# =====================================================================
# 28_console_pvalues.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The five p-values quoted in the text of Section 3.2 that until now
#   existed only as console output. Every one is recomputed from the
#   CSV files of 02_resultados, so that each number in the manuscript
#   traces back to a file.
#
#   1. Robust vs classical FAR, Table 1: two-proportion z test and
#      Fisher exact test on the false-alarm counts over 200,000 batches
#      each.                                  (text: p = 0.20 and 0.22)
#   2. Homogeneity of the FAR across K = 30 ... 500, Table 2:
#      chi-square test of equal proportions.  (text: p = 0.65)
#   3. Trend of the FAR with K, Table 2: Cochran-Armitage trend test
#      with K as the score.                   (text: p = 0.95)
#   4. Clean vs contaminated Phase 1, Table 3: two-proportion z test on
#      the pooled counts over 1,000,000 batches each.   (text: p = 0.54)
#
# INPUT
#   02_resultados/table1_ztest_faithful.csv
#   02_resultados/arl0_convergence_K_v3.csv
#   02_resultados/arl0_clean_vs_contaminated_v2.csv
# OUTPUT
#   02_resultados/console_pvalues.csv
# =====================================================================
DIR_OUT <- "02_resultados"

t1 <- read.csv(file.path(DIR_OUT, "table1_ztest_faithful.csv"))
t2 <- read.csv(file.path(DIR_OUT, "arl0_convergence_K_v3.csv"))
t3 <- read.csv(file.path(DIR_OUT, "arl0_clean_vs_contaminated_v2.csv"))

# --- 1. Table 1: robust vs classical ---------------------------------
fa  <- t1$FA_total
n   <- t1$n_batches
z_test  <- prop.test(fa, n, correct = FALSE)
fisher  <- fisher.test(cbind(fa, n - fa))

# --- 2 and 3. Table 2: homogeneity and trend across K ----------------
homog <- prop.test(t2$false_alarms, t2$total_batches, correct = FALSE)
trend <- prop.trend.test(t2$false_alarms, t2$total_batches, score = t2$K)

# --- 4. Table 3: clean vs contaminated Phase 1 ------------------------
cond_col <- if ("condition" %in% names(t3)) "condition" else "condicion"
fa3 <- tapply(t3$false_alarms, t3[[cond_col]], sum)
n3  <- tapply(t3$total_batches, t3[[cond_col]], sum)
cvc <- prop.test(as.numeric(fa3), as.numeric(n3), correct = FALSE)

# --- Collect ----------------------------------------------------------
out <- data.frame(
  test = c("Table 1: two-proportion z, robust vs classical",
           "Table 1: Fisher exact, robust vs classical",
           "Table 2: chi-square homogeneity across K",
           "Table 2: Cochran-Armitage trend in K",
           "Table 3: two-proportion z, clean vs contaminated"),
  statistic = c(unname(z_test$statistic), NA,
                unname(homog$statistic), unname(trend$statistic),
                unname(cvc$statistic)),
  df = c(unname(z_test$parameter), NA, unname(homog$parameter),
         unname(trend$parameter), unname(cvc$parameter)),
  p_value = c(z_test$p.value, fisher$p.value, homog$p.value,
              trend$p.value, cvc$p.value),
  text_value = c(0.20, 0.22, 0.65, 0.95, 0.54)
)
out$p_rounded <- round(out$p_value, 2)
out$match     <- ifelse(out$p_rounded == out$text_value, "OK", "REVIEW")

cat("\n===== P-VALUES QUOTED IN SECTION 3.2 =====\n")
print(out[, c("test", "p_value", "p_rounded", "text_value", "match")],
      row.names = FALSE, digits = 4)

write.csv(out, file.path(DIR_OUT, "console_pvalues.csv"), row.names = FALSE)
cat(sprintf("\nSaved to %s/console_pvalues.csv\n", DIR_OUT))