test_that("ucl_F_adjusted matches the closed-form value (m* = MCD subset size)", {
  # Minimal calibration stub: only the pieces ucl_F_adjusted needs.
  cal <- list(
    Sw        = diag(4),
    mcd_alpha = 0.67,
    weights   = setNames(rep(1 / 30, 30), paste0("B", 1:30))
  )

  ucl <- ucl_F_adjusted(cal, I = 20)   # alpha default = 0.001, m_star = "subset"

  # Reconstruct the paper's formula:
  # UCL = [J (K+1) (m*-1) / (K m* - K - J + 1)] * F_{J, df2, 1-alpha}
  J      <- 4L
  K      <- 30L
  alpha  <- 0.001
  # MCD subset size for I = 20, J = 4, alpha = 0.67:
  #   n2 = (20 + 4 + 1) %/% 2 = 12;  floor(2*12 - 20 + 2*(20 - 12)*0.67) = floor(14.72) = 14
  m_star <- 14
  df2    <- K * m_star - K - J + 1                   # 387
  scale  <- J * (K + 1) * (m_star - 1) / df2         # 1612 / 387
  Fq     <- stats::qf(1 - alpha, df1 = J, df2 = df2)
  expected_UCL <- scale * Fq

  expect_equal(ucl$UCL, expected_UCL)
  expect_equal(ucl$parameters$m_star, m_star)
  expect_equal(ucl$parameters$m_star_rule, "subset")
  expect_equal(ucl$parameters$df2,    df2)
  expect_equal(ucl$parameters$alpha,  alpha)
  expect_gt(ucl$UCL, 0)
})

test_that("m* = 'subset' is the subset size covMcd really uses", {
  # The point of the default: m* is the number of observations the MCD of a
  # batch is computed from, which robustbase reports as $quan.
  set.seed(1)
  X <- matrix(stats::rnorm(20 * 4), 20, 4)
  quan <- robustbase::covMcd(X, alpha = 0.67)$quan
  cal <- list(Sw = diag(4), mcd_alpha = 0.67,
              weights = setNames(rep(1 / 30, 30), paste0("B", 1:30)))
  expect_equal(ucl_F_adjusted(cal, I = 20)$parameters$m_star, as.numeric(quan))
  expect_equal(as.numeric(quan), 14)
})

test_that("m_star = 'nominal' reproduces the limit of version 0.2.0", {
  cal <- list(Sw = diag(4), mcd_alpha = 0.67,
              weights = setNames(rep(1 / 30, 30), paste0("B", 1:30)))
  ucl <- ucl_F_adjusted(cal, I = 20, m_star = "nominal")
  m_star <- round(20 * 0.67)                         # 13
  df2    <- 30 * m_star - 30 - 4 + 1                 # 357
  expected <- 4 * 31 * (m_star - 1) / df2 * stats::qf(0.999, 4, df2)
  expect_equal(ucl$UCL, expected)
  expect_equal(ucl$parameters$m_star, 13)
  expect_equal(ucl$parameters$df2, 357)
  expect_equal(ucl$parameters$m_star_rule, "nominal")
  # The subset rule gives a lower limit (more degrees of freedom).
  expect_lt(ucl_F_adjusted(cal, I = 20)$UCL, ucl$UCL)
})

test_that("ucl_F_adjusted validates its inputs", {
  cal <- list(
    Sw = diag(4), mcd_alpha = 0.67,
    weights = setNames(rep(1 / 30, 30), paste0("B", 1:30))
  )
  expect_error(ucl_F_adjusted(cal, I = 1),                   "I' must be")
  expect_error(ucl_F_adjusted(cal, I = 20, alpha = 1.5),     "alpha' must be")
  expect_error(ucl_F_adjusted(cal, I = 20, m_star = "mean"), "should be one of")
  expect_error(ucl_F_adjusted("nope", I = 20),
               "calibrate_afm_mcd")
})
