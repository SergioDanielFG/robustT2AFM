# =====================================================================
# SCRIPT_M1c_recount_mstar14.R
# ---------------------------------------------------------------------
# Revision round 1, module M1, third part (Reviewer 1, comments 2 and 3;
# Tables 1, 2, 3 and 13; decision D12 of REGISTRO_CAMBIOS_v9.md).
# Nothing is simulated: it recounts the false alarms of the T2 saved by
# SCRIPT_M1_limit.R in 03_T2_crudos with the limit of Eq. (8) at
# m* = 14 (size of the MCD subset), the value adopted in the revision.
#
# WHAT IT COMPUTES
#   A. Every cell of M1 (Tables 1, 2, 3 and 13), every chart: the four
#      MFA-MCD variants with Eq. (8) at m* = 13 (as published) and at
#      m* = 14 (final); the classical chart with its own limit. False
#      alarms, FAR, SE of the FAR across replicates, ARL0 and exact
#      Poisson 95% interval, with the same formulas as M1.
#   B. The tests quoted in Section 3.2, with the same code as
#      28_console_pvalues.R and 16_ci_paired_far_difference.R, run for
#      V7 with m* = 13 (must reproduce the manuscript) and for NEW with
#      m* = 14 (values for the revision).
#   C. Empirical 0.999 limit of NEW and its ratio to Eq. (8) with
#      m* = 14 (the factor c of M1b was computed against m* = 13).
#
# ANCHORS (nothing is written if any fails)
#   A1 m* = 13: false alarms of every chart and cell equal
#      table_M1_T1_calibration, table_M1_T2_K and table_M1_T3_seeds.
#   A2 m* = 14: NEW equals table_M1b_mstar_arl0 (T1 164; pooled clean
#      995, contaminated 804; K = 50/100/200/500: 128/132/120/123).
#   A3 Tests with V7 m* = 13 equal console_pvalues.csv and
#      ci_paired_far.csv.
#   A4 Empirical limit of NEW equals table_M1b_empirical_limit.
#
# OUTPUT (05_revision_R1/02_resultados/)
#   table_M1c_cells_mstar.csv       one row per cell x chart x m*
#   table_M1c_tests.csv             tests of B
#   table_M1c_empirical_limit.csv   part C
#   Run from C:/temp_paper. Expected time: under one minute.
# =====================================================================
DIR_RAW <- "05_revision_R1/03_T2_crudos"
DIR_OUT <- "05_revision_R1/02_resultados"
I <- 20; J <- 4; ALPHA <- 0.001; N_F2 <- 100; N_REP <- 2000
AFM <- c("V7", "SM_std", "MCDc_V7w", "NEW")

ucl_m <- function(K, m) { d2 <- K * m - K - J + 1
  (J * (K + 1) * (m - 1) / d2) * qf(1 - ALPHA, J, d2) }

files <- list.files(DIR_RAW, pattern = "^M1_.*\\.rds$", full.names = TRUE)
raw   <- setNames(lapply(files, readRDS), sub("^M1_(.*)\\.rds$", "\\1", basename(files)))
cat("Cells read:", length(raw), "\n")
stopifnot(length(raw) == 15)
cellK  <- function(K) if (K == 30) "T3_b2026_clean" else sprintf("T2_K%d_clean", K)
ucl_of <- function(cc, m, ms) if (m == "CLASSICAL") unname(raw[[cc]]$ucl["CLASSICAL"]) else ucl_m(raw[[cc]]$cell$K, ms)
fa_rep <- function(cc, m, ms) rowSums(raw[[cc]]$T2[, , m] > ucl_of(cc, m, ms))  # false alarms per replicate

# --- A. Recount per cell ---------------------------------------------------
summ <- function(cc, m, ms) {
  f <- fa_rep(cc, m, ms); fa <- sum(f); n <- length(raw[[cc]]$T2[, , m])
  data.frame(cell = cc, table = raw[[cc]]$cell$table, K = raw[[cc]]$cell$K,
             seed_base = raw[[cc]]$cell$base,
             phase1 = if (raw[[cc]]$cell$ob > 0) "contaminated" else "clean",
             method = m, m_star = ms, UCL = round(ucl_of(cc, m, ms), 4),
             false_alarms = fa, batches = n, FAR = fa / n,
             SE_FAR = sd(f / N_F2) / sqrt(length(f)),
             ARL0 = if (fa > 0) round(n / fa) else NA,
             CI95_low = round(n / (qchisq(0.975, 2 * fa + 2) / 2)),
             CI95_high = if (fa > 0) round(n / (qchisq(0.025, 2 * fa) / 2)) else Inf)
}
A <- do.call(rbind, lapply(names(raw), function(cc) rbind(
  do.call(rbind, lapply(AFM, function(m) rbind(summ(cc, m, 13), summ(cc, m, 14)))),
  summ(cc, "CLASSICAL", NA))))
rownames(A) <- NULL

# --- Anchor A1: m* = 13 reproduces M1 ----------------------------------------
get13 <- function(cc, m) A$false_alarms[A$cell == cc & A$method == m & (is.na(A$m_star) | A$m_star == 13)]
t1 <- read.csv(file.path(DIR_OUT, "table_M1_T1_calibration.csv"))
t2 <- read.csv(file.path(DIR_OUT, "table_M1_T2_K.csv"))
t3 <- read.csv(file.path(DIR_OUT, "table_M1_T3_seeds.csv"))
a1 <- c(
  mapply(function(m, fa) get13("T1_K30_contam", m) == fa, t1$method, t1$false_alarms),
  mapply(function(K, m, fa) get13(cellK(K), m) == fa, t2$K, t2$method, t2$false_alarms),
  mapply(function(b, ph, m, fa) get13(sprintf("T3_b%d_%s", b, if (ph == "clean") "clean" else "contam"), m) == fa,
         t3$seed_base, t3$phase1, t3$method, t3$false_alarms))
cat(sprintf("A1 (m* = 13 reproduces M1): %d of %d counts match -> %s\n",
            sum(a1), length(a1), if (all(a1)) "OK" else "*** FAIL ***"))

# --- Anchor A2: m* = 14 reproduces M1b ---------------------------------------
get14  <- function(cc) A$false_alarms[A$cell == cc & A$method == "NEW" & A$m_star %in% 14]
pool14 <- function(ph) sum(sapply(grep(sprintf("^T3_.*_%s$", ph), names(raw), value = TRUE), get14))
g2 <- c(get14("T1_K30_contam"), pool14("clean"), pool14("contam"), get14("T2_K50_clean"),
        get14("T2_K100_clean"), get14("T2_K200_clean"), get14("T2_K500_clean"))
a2 <- all(g2 == c(164, 995, 804, 128, 132, 120, 123))
cat(sprintf("A2 (NEW m* = 14 vs M1b): %s  expected 164/995/804/128/132/120/123 -> %s\n",
            paste(g2, collapse = "/"), if (a2) "OK" else "*** FAIL ***"))

# --- B. Tests of Section 3.2 -------------------------------------------------
tests <- function(meth, ms) {
  n1 <- N_REP * N_F2
  fr <- fa_rep("T1_K30_contam", meth, ms); fc <- fa_rep("T1_K30_contam", "CLASSICAL", NA)
  z  <- prop.test(c(sum(fr), sum(fc)), c(n1, n1), correct = FALSE)
  fi <- fisher.test(cbind(c(sum(fr), sum(fc)), n1 - c(sum(fr), sum(fc))))
  rob <- fr / N_F2; cls <- fc / N_F2
  tt <- t.test(rob, cls, paired = TRUE)
  Ks  <- c(30, 50, 100, 200, 500)
  faK <- sapply(Ks, function(K) sum(fa_rep(cellK(K), meth, ms)))
  ho  <- prop.test(faK, rep(n1, 5), correct = FALSE)
  tr  <- prop.trend.test(faK, rep(n1, 5), score = Ks)
  fa3 <- sapply(c("clean", "contam"), function(ph)
    sum(sapply(grep(sprintf("^T3_.*_%s$", ph), names(raw), value = TRUE), function(cc) sum(fa_rep(cc, meth, ms)))))
  cv  <- prop.test(fa3, rep(5 * n1, 2), correct = FALSE)
  data.frame(chart = meth, m_star = ms,
    test = c("Table 1: two-proportion z, robust vs classical",
             "Table 1: Fisher exact, robust vs classical",
             "Table 1: paired t, FAR robust - classical",
             "Table 2: chi-square homogeneity across K",
             "Table 2: Cochran-Armitage trend in K",
             "Table 3: two-proportion z, clean vs contaminated"),
    counts = c(rep(paste(sum(fr), sum(fc), sep = "/"), 3), rep(paste(faK, collapse = "/"), 2),
               paste(fa3, collapse = "/")),
    statistic = unname(c(z$statistic, NA, tt$statistic, ho$statistic, tr$statistic, cv$statistic)),
    p_value   = c(z$p.value, fi$p.value, tt$p.value, ho$p.value, tr$p.value, cv$p.value),
    difference  = c(NA, NA, mean(rob - cls), NA, NA, NA),
    SE          = c(NA, NA, sd(rob - cls) / sqrt(length(rob)), NA, NA, NA),
    CI95_low    = c(NA, NA, tt$conf.int[1], NA, NA, NA),
    CI95_high   = c(NA, NA, tt$conf.int[2], NA, NA, NA),
    correlation = c(NA, NA, cor(rob, cls), NA, NA, NA))
}
B <- rbind(tests("V7", 13), tests("NEW", 14))

# --- Anchor A3: V7 m* = 13 reproduces the published tests ---------------------
v7 <- B[B$chart == "V7", ]
a3 <- c(abs(v7$p_value[1] - 0.201923175404219) < 1e-9,
        abs(v7$p_value[2] - 0.224174725743826) < 1e-9,
        abs(v7$p_value[4] - 0.650720392664498) < 1e-9,
        abs(v7$p_value[5] - 0.949659623924169) < 1e-9,
        abs(v7$p_value[6] - 0.543601303103876) < 1e-9,
        abs(v7$difference[3] - (-0.000105)) < 1e-10,
        abs(v7$SE[3] - 5.97601299912265e-05) < 1e-10)
cat(sprintf("A3 (V7 tests vs console_pvalues and ci_paired_far): %d of 7 -> %s\n",
            sum(a3), if (all(a3)) "OK" else "*** FAIL ***"))

# --- C. Empirical limit of NEW against Eq. (8) with m* = 14 -------------------
groups <- list(clean_pooled = grep("^T3_.*_clean$", names(raw), value = TRUE),
               contaminated_pooled = grep("^T3_.*_contam$", names(raw), value = TRUE),
               T1_K30_contam = "T1_K30_contam", T2_K50_clean = "T2_K50_clean",
               T2_K100_clean = "T2_K100_clean", T2_K200_clean = "T2_K200_clean",
               T2_K500_clean = "T2_K500_clean")
C <- do.call(rbind, lapply(names(groups), function(g) {
  K <- raw[[groups[[g]][1]]]$cell$K
  x <- unlist(lapply(groups[[g]], function(cc) as.vector(raw[[cc]]$T2[, , "NEW"])))
  q <- unname(quantile(x, 1 - ALPHA, type = 8)); u <- ucl_m(K, 14)
  data.frame(group = g, K = K, batches = length(x), UCL_eq8_m14 = round(u, 4),
             UCL_empirical = round(q, 4), factor_c14 = round(q / u, 4))
}))
eb <- read.csv(file.path(DIR_OUT, "table_M1b_empirical_limit.csv")); eb <- eb[eb$method == "NEW", ]
a4 <- all(sapply(seq_len(nrow(C)), function(i)
  abs(C$UCL_empirical[i] - eb$UCL_empirical[eb$group == C$group[i]]) < 1e-4))
cat(sprintf("A4 (empirical limit of NEW vs M1b): %s\n", if (a4) "OK" else "*** FAIL ***"))

# --- Output --------------------------------------------------------------------
ok <- all(a1) && a2 && all(a3) && a4
if (!ok) stop("ANCHOR FAILURE: nothing written. Copy the whole output.")
write.csv(A, file.path(DIR_OUT, "table_M1c_cells_mstar.csv"),     row.names = FALSE)
write.csv(B, file.path(DIR_OUT, "table_M1c_tests.csv"),           row.names = FALSE)
write.csv(C, file.path(DIR_OUT, "table_M1c_empirical_limit.csv"), row.names = FALSE)

vista <- function(x) print(x[, c("cell", "method", "m_star", "UCL", "false_alarms", "FAR",
                                 "SE_FAR", "ARL0", "CI95_low", "CI95_high")], row.names = FALSE, digits = 4)
keep <- function(x) x[(x$method == "V7" & x$m_star %in% 13) | (x$method == "NEW" & x$m_star %in% 14) |
                      x$method == "CLASSICAL", ]
cat("\n=== Table 1 (T1 cell): V7 m*=13, NEW m*=14, classical ===\n"); vista(keep(A[A$cell == "T1_K30_contam", ]))
cat("\n=== Table 2 (clean, K = 30 ... 500) ===\n")
vista(keep(A[A$cell %in% sapply(c(30, 50, 100, 200, 500), cellK), ]))
cat("\n=== Table 3 (K = 30, five seed bases) ===\n"); vista(keep(A[A$table == "T3", ]))
cat("\n=== Table 13 (base 2026): all MFA-MCD variants with m* = 14 ===\n")
vista(A[grepl("^T3_b2026_", A$cell) & (A$m_star %in% 14 | A$method == "CLASSICAL"), ])
cat("\n=== Tests of Section 3.2 ===\n")
print(B[, c("chart", "m_star", "test", "counts", "statistic", "p_value", "difference", "SE",
            "CI95_low", "CI95_high")], row.names = FALSE, digits = 4)
cat("\n=== Empirical limit of NEW and factor c with m* = 14 ===\n"); print(C, row.names = FALSE)
cat("\nALL ANCHORS OK: valid. Files written in", DIR_OUT, "\n")
