# =====================================================================
# SCRIPT_M1b_limit_posthoc.R
# ---------------------------------------------------------------------
# Revision round 1, module M1, second part: analysis of the control
# limit on the T2 saved by SCRIPT_M1_limit.R (Reviewer 1, comments 1,
# 2 and 3). Nothing is simulated: it only recounts false alarms.
#
# WHAT IT COMPUTES
#   A. ARL0 of Eq. (8) with different effective batch sizes m*:
#        13    = round(I * h), the value of the manuscript
#        14    = raw MCD subset (quan), measured in M1
#        16.7  = observations kept after reweighting, measured in M1
#                (mean over clean batches; Eq. (8) accepts a
#                non-integer m*)
#        17    = its rounded value
#        20    = the whole batch (no trimming)
#   B. Empirical limit: the 0.999 quantile of the in-control T2, and
#      the correction factor c = empirical limit / Eq. (8) limit with
#      m* = 13. With 1 000 000 batches per condition the 0.999 quantile
#      rests on about 1000 exceedances.
#   Both for the published method (V7) and for the final one (NEW), and
#   for the isolated center (MCDc_V7w) and the spatial median (SM_std).
#
# GROUPS
#   clean_pooled / contaminated_pooled : the five T3 cells of each
#     condition (K = 30), 1 000 000 batches each
#   T1_K30_contam                       : Table 1 cell
#   K = 50, 100, 200, 500 (clean)       : Table 2 cells
#
# ANCHOR
#   With m* = 13 the false alarms must equal those of SCRIPT_M1_limit.R:
#   clean_pooled V7 729 / NEW 974, contaminated_pooled V7 706 / NEW 795.
#
# OUTPUT
#   05_revision_R1/02_resultados/table_M1b_mstar_arl0.csv
#   05_revision_R1/02_resultados/table_M1b_empirical_limit.csv
#   Run from C:/temp_paper.
# =====================================================================
DIR_RAW <- "05_revision_R1/03_T2_crudos"
DIR_OUT <- "05_revision_R1/02_resultados"
I <- 20; J <- 4; ALPHA <- 0.001
M_GRID  <- c(13, 14, 16.7, 17, 20)
METHODS <- c("V7", "NEW", "MCDc_V7w", "SM_std")

ucl_m <- function(K, m) { d2 <- K * m - K - J + 1
  (J * (K + 1) * (m - 1) / d2) * qf(1 - ALPHA, J, d2) }
arl_ci <- function(fa, n) c(ARL0 = if (fa > 0) round(n / fa) else NA,
  CI95_low = round(n / (qchisq(0.975, 2 * fa + 2) / 2)),
  CI95_high = if (fa > 0) round(n / (qchisq(0.025, 2 * fa) / 2)) else Inf)

files <- list.files(DIR_RAW, pattern = "^M1_.*\\.rds$", full.names = TRUE)
raw   <- setNames(lapply(files, readRDS), sub("^M1_(.*)\\.rds$", "\\1", basename(files)))
cat("Cells read:", length(raw), "\n")
stopifnot(length(raw) == 15)

grupos <- c(
  list(clean_pooled        = grep("^T3_.*_clean$",  names(raw), value = TRUE),
       contaminated_pooled = grep("^T3_.*_contam$", names(raw), value = TRUE),
       T1_K30_contam       = "T1_K30_contam"),
  setNames(as.list(grep("^T2_", names(raw), value = TRUE)), grep("^T2_", names(raw), value = TRUE)))

t2_de <- function(cells, m) unlist(lapply(cells, function(cc) as.vector(raw[[cc]]$T2[, , m])))

filasA <- list(); filasB <- list()
for (g in names(grupos)) {
  cc <- grupos[[g]]; K <- raw[[cc[1]]]$cell$K
  for (m in METHODS) {
    x <- t2_de(cc, m); n <- length(x)
    for (ms in M_GRID) {
      u <- ucl_m(K, ms); fa <- sum(x > u)
      filasA[[length(filasA) + 1]] <- data.frame(group = g, K = K, method = m, m_star = ms,
        UCL = round(u, 4), false_alarms = fa, batches = n, t(arl_ci(fa, n)))
    }
    q  <- unname(quantile(x, 1 - ALPHA, type = 8))
    u13 <- ucl_m(K, 13)
    filasB[[length(filasB) + 1]] <- data.frame(group = g, K = K, method = m, batches = n,
      UCL_eq8_m13 = round(u13, 4), UCL_empirical = round(q, 4), factor_c = round(q / u13, 4),
      m_star_equivalent = tryCatch(round(uniroot(function(mm) ucl_m(K, mm) - q, c(5.5, 1000))$root, 2),
                                    error = function(e) NA))   # NA: no m* reaches it
  }
}
A <- do.call(rbind, filasA); B <- do.call(rbind, filasB)

# --- Anchor ----------------------------------------------------------------
get_fa <- function(g, m) A$false_alarms[A$group == g & A$method == m & A$m_star == 13]
anc <- c(get_fa("clean_pooled", "V7"), get_fa("clean_pooled", "NEW"),
         get_fa("contaminated_pooled", "V7"), get_fa("contaminated_pooled", "NEW"))
ok <- all(anc == c(729, 974, 706, 795))
cat(sprintf("Anchor (m* = 13): %s  expected 729/974/706/795 -> %s\n\n",
            paste(anc, collapse = "/"), if (ok) "OK" else "*** FAIL ***"))

write.csv(A, file.path(DIR_OUT, "table_M1b_mstar_arl0.csv"),      row.names = FALSE)
write.csv(B, file.path(DIR_OUT, "table_M1b_empirical_limit.csv"), row.names = FALSE)

cat("=== A. ARL0 of Eq. (8) by m* (V7 and NEW) ===\n")
w <- reshape(A[A$method %in% c("V7", "NEW"), c("group", "method", "m_star", "ARL0")],
             idvar = c("group", "method"), timevar = "m_star", direction = "wide")
names(w) <- sub("ARL0\\.", "m*=", names(w)); print(w, row.names = FALSE)
cat("\n=== B. Empirical limit and correction factor ===\n")
print(B[B$method %in% c("V7", "NEW"), ], row.names = FALSE)
cat(if (ok) "\nAnchor OK: valid.\n" else "\nANCHOR FAILED: do not use.\n")
