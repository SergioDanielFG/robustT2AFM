# Generates the two datasets shipped with the package, afm_phase1 and
# afm_phase2. Run from the package root:
#
#   source("data-raw/make_datasets.R")
#
# The call below is the base configuration of Frutos-Galarza et al. (2026)
# and is reproduced verbatim in the @source of R/data.R and in
# tests/testthat/test-data.R, which regenerates the data and compares it
# with what is shipped. If you change anything here, change it there too.

devtools::load_all(".", quiet = TRUE)

sim <- simulate_batch_process(
  K1 = 30, K2 = 20, I = 20, J = 4, rho = 0.6,
  outlier_batches_F1 = 6, outlier_rate = 0.20, outlier_shift = 4,
  prop_ooc_F2 = 0.5, shift_ooc = 1.0,
  seed = 20260425
)

# Split by phase, drop the now-unused factor levels and renumber the rows,
# so each object looks like a data set a user would load rather than like a
# slice of a larger one.
afm_phase1 <- droplevels(subset(sim, Phase == "Phase 1"))
afm_phase2 <- droplevels(subset(sim, Phase == "Phase 2"))
rownames(afm_phase1) <- NULL
rownames(afm_phase2) <- NULL

save(afm_phase1, file = "data/afm_phase1.rda", compress = "xz", version = 2)
save(afm_phase2, file = "data/afm_phase2.rda", compress = "xz", version = 2)

message("afm_phase1: ", nrow(afm_phase1), " rows, ",
        round(file.size("data/afm_phase1.rda") / 1024, 1), " KB")
message("afm_phase2: ", nrow(afm_phase2), " rows, ",
        round(file.size("data/afm_phase2.rda") / 1024, 1), " KB")
