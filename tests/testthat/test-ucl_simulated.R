# The simulated limit recalibrates the method B times, so these tests use a
# small B and alpha = 0.01: 5 x 200 = 1000 simulated T2, 10 beyond the limit.

cal_sim <- function() {
  data(afm_phase1, package = "robustT2AFM", envir = environment())
  calibrate_afm_mcd(afm_phase1, paste0("Var", 1:4))
}

test_that("ucl_simulated returns a limit with the documented structure", {
  cal <- cal_sim()
  us  <- ucl_simulated(cal, I = 20, alpha = 0.01, B = 5, seed = 1)

  expect_named(us, c("UCL", "method", "parameters"))
  expect_equal(us$method, "Simulated (parametric bootstrap)")
  expect_true(is.finite(us$UCL))
  # A 0.99 quantile of a T2 with J = 4: above the chi-square value (13.3),
  # well below the classical 0.999 limit territory.
  expect_gt(us$UCL, 10)
  expect_lt(us$UCL, 30)
  expect_equal(us$parameters$n_T2, 5 * 200)
  expect_equal(us$parameters$K, 30)
  expect_equal(us$parameters$J, 4)
  # The synthetic samples are calibrated with the options of the real one.
  expect_equal(us$parameters$scaling, cal$scaling)
  expect_equal(us$parameters$center, cal$center)
  expect_equal(us$parameters$mcd_alpha, cal$mcd_alpha)
})

test_that("ucl_simulated is reproducible with a seed and leaves the RNG alone", {
  cal <- cal_sim()
  u1 <- ucl_simulated(cal, I = 20, alpha = 0.01, B = 5, seed = 7)$UCL
  u2 <- ucl_simulated(cal, I = 20, alpha = 0.01, B = 5, seed = 7)$UCL
  u3 <- ucl_simulated(cal, I = 20, alpha = 0.01, B = 5, seed = 8)$UCL
  expect_identical(u1, u2)
  expect_false(identical(u1, u3))

  # With a seed, the user's random number stream is restored afterwards.
  set.seed(123); a <- stats::runif(3)
  set.seed(123); invisible(ucl_simulated(cal, I = 20, alpha = 0.01, B = 5, seed = 7))
  b <- stats::runif(3)
  expect_identical(a, b)
})

test_that("ucl_simulated warns when few simulated values lie beyond the limit", {
  cal <- cal_sim()
  expect_warning(ucl_simulated(cal, I = 20, alpha = 0.001, B = 5, seed = 1),
                 "imprecise")
})

test_that("ucl_simulated validates its inputs", {
  cal <- cal_sim()
  expect_error(ucl_simulated("nope", I = 20), "calibrate_afm_mcd")
  expect_error(ucl_simulated(cal, I = 1), "I' must be")
  expect_error(ucl_simulated(cal, I = 4), "larger than the number of variables")
  expect_error(ucl_simulated(cal, I = 20, alpha = 1.5), "alpha' must be")
  expect_error(ucl_simulated(cal, I = 20, B = 0), "'B'")
  expect_error(ucl_simulated(cal, I = 20, n_new = 2.5), "'n_new'")
  expect_error(ucl_simulated(cal, I = 20, seed = "a"), "'seed'")
})
