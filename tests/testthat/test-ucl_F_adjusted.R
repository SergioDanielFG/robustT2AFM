test_that("ucl_F_adjusted matches the closed-form value", {
  # Minimal calibration stub: only the pieces ucl_F_adjusted needs.
  cal <- list(
    Sw        = diag(4),
    mcd_alpha = 0.67,
    weights   = setNames(rep(1 / 30, 30), paste0("B", 1:30))
  )

  ucl <- ucl_F_adjusted(cal, I = 20)   # alpha default = 0.001

  # Reconstruct the paper's formula:
  # UCL = [J (K+1) (m*-1) / (K m* - K - J + 1)] * F_{J, df2, 1-alpha}
  J      <- 4L
  K      <- 30L
  I      <- 20L
  h      <- 0.67
  alpha  <- 0.001
  m_star <- round(I * h)                             # 13
  df2    <- K * m_star - K - J + 1                   # 357
  scale  <- J * (K + 1) * (m_star - 1) / df2         # 1488 / 357
  Fq     <- stats::qf(1 - alpha, df1 = J, df2 = df2)
  expected_UCL <- scale * Fq

  expect_equal(ucl$UCL, expected_UCL)
  expect_equal(ucl$parameters$m_star, m_star)
  expect_equal(ucl$parameters$df2,    df2)
  expect_equal(ucl$parameters$alpha,  alpha)
  expect_gt(ucl$UCL, 0)
})

test_that("ucl_F_adjusted validates its inputs", {
  cal <- list(
    Sw = diag(4), mcd_alpha = 0.67,
    weights = setNames(rep(1 / 30, 30), paste0("B", 1:30))
  )
  expect_error(ucl_F_adjusted(cal, I = 1),                 "I' must be")
  expect_error(ucl_F_adjusted(cal, I = 20, alpha = 1.5),   "alpha' must be")
  expect_error(ucl_F_adjusted("nope", I = 20),
               "calibrate_afm_mcd")
})
