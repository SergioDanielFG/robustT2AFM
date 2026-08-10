vars <- paste0("Var", 1:4)

test_that("hotelling_classical_monitor adds is_ooc only when ucl is supplied", {

  cal <- hotelling_classical_calibrate(afm_phase1, vars)

  mon_plain <- hotelling_classical_monitor(afm_phase2, cal, vars)
  expect_named(mon_plain, c("Batch", "I", "T2"))
  expect_false("is_ooc" %in% names(mon_plain))

  u       <- hotelling_classical_ucl(K = 30, I = 20, J = 4)$UCL
  mon_ucl <- hotelling_classical_monitor(afm_phase2, cal, vars, ucl = u)
  expect_named(mon_ucl, c("Batch", "I", "T2", "is_ooc"))
  expect_type(mon_ucl$is_ooc, "logical")

  # El limite solo anade la columna: los T2 no se tocan.
  expect_identical(mon_ucl$T2, mon_plain$T2)
  expect_identical(mon_ucl[c("Batch", "I", "T2")], mon_plain)

  # ucl = NULL reproduce la salida de tres columnas exactamente.
  expect_identical(hotelling_classical_monitor(afm_phase2, cal, vars,
                                               ucl = NULL),
                   mon_plain)
})

test_that("hotelling_classical_monitor flags with a strict greater-than", {

  cal <- hotelling_classical_calibrate(afm_phase1, vars)
  mon <- hotelling_classical_monitor(afm_phase2, cal, vars)

  # Un lote justo sobre el limite queda bajo control; hay que marcarlo solo
  # por encima. Se toma el T2 de un lote real como limite exacto.
  on_limit <- mon$T2[1]
  flags    <- hotelling_classical_monitor(afm_phase2, cal, vars,
                                          ucl = on_limit)$is_ooc
  expect_false(flags[1])

  # Un pelo por debajo del mismo valor si lo marca.
  just_under <- hotelling_classical_monitor(afm_phase2, cal, vars,
                                            ucl = on_limit * (1 - 1e-9))$is_ooc
  expect_true(just_under[1])

  # Y la columna es exactamente la comparacion estricta.
  u <- hotelling_classical_ucl(K = 30, I = 20, J = 4)$UCL
  expect_identical(
    hotelling_classical_monitor(afm_phase2, cal, vars, ucl = u)$is_ooc,
    mon$T2 > u
  )
})

test_that("hotelling_classical_monitor validates ucl", {

  cal <- hotelling_classical_calibrate(afm_phase1, vars)

  expect_error(hotelling_classical_monitor(afm_phase2, cal, vars, ucl = -1),
               "must be a single positive")
  expect_error(hotelling_classical_monitor(afm_phase2, cal, vars, ucl = 0),
               "must be a single positive")
  expect_error(hotelling_classical_monitor(afm_phase2, cal, vars,
                                           ucl = c(10, 20)),
               "must be a single positive")
  expect_error(hotelling_classical_monitor(afm_phase2, cal, vars,
                                           ucl = "19.5"),
               "must be a single positive")
  expect_error(hotelling_classical_monitor(afm_phase2, cal, vars,
                                           ucl = NA_real_),
               "must be a single positive")

  # Pasar el objeto entero en vez de su $UCL es el error previsible, y el
  # mensaje tiene que nombrar la salida.
  expect_error(
    hotelling_classical_monitor(afm_phase2, cal, vars,
                                ucl = hotelling_classical_ucl(K = 30, I = 20,
                                                              J = 4)),
    "hotelling_classical_ucl\\(K, I, J\\)\\$UCL"
  )
})

test_that("vectorising the classical loop returns an identical object", {

  # Reimplementacion literal del bucle con rbind() anterior.
  old_rbind_loop <- function(new_data, calibration, variables,
                             batch_col = "Batch") {
    mu_global <- calibration$mu_global
    Sp_inv    <- solve(calibration$Sp)
    batch_id  <- new_data[[batch_col]]
    batches   <- unique(batch_id)

    results <- data.frame(Batch = character(), I = integer(), T2 = numeric(),
                          stringsAsFactors = FALSE)
    for (batch in batches) {
      subset_batch <- new_data[batch_id == batch, variables, drop = FALSE]
      I     <- nrow(subset_batch)
      x_bar <- if (I == 1) {
        as.numeric(subset_batch[1, ])
      } else {
        colMeans(as.matrix(subset_batch))
      }
      diff <- x_bar - mu_global
      T2   <- as.numeric(I * t(diff) %*% Sp_inv %*% diff)
      results <- rbind(results, data.frame(Batch = as.character(batch), I = I,
                                           T2 = T2, stringsAsFactors = FALSE))
    }
    results
  }

  cal <- hotelling_classical_calibrate(afm_phase1, vars)
  expect_identical(hotelling_classical_monitor(afm_phase2, cal, vars),
                   old_rbind_loop(afm_phase2, cal, vars))

  u   <- hotelling_classical_ucl(K = 30, I = 20, J = 4)$UCL
  new <- hotelling_classical_monitor(afm_phase2, cal, vars, ucl = u)
  old <- old_rbind_loop(afm_phase2, cal, vars)
  old$is_ooc <- old$T2 > u
  expect_identical(new, old)

  # El caso I == 1, que la gemela clasica trata aparte y la robusta no.
  one <- afm_phase2[c(1, which(afm_phase2$Batch == "F2_B02")[1]), ]
  expect_identical(hotelling_classical_monitor(one, cal, vars),
                   old_rbind_loop(one, cal, vars))
  expect_identical(hotelling_classical_monitor(one, cal, vars)$I, c(1L, 1L))
})

test_that("an invalid classical batch aborts at the same point too", {
  cal <- hotelling_classical_calibrate(afm_phase1, vars)

  d_last <- afm_phase2
  d_last$Var1[which(d_last$Batch == "F2_B20")[1]] <- NA
  expect_error(hotelling_classical_monitor(d_last, cal, vars),
               "Batch 'F2_B20' contains non-finite values")

  d_two <- afm_phase2
  d_two$Var1[which(d_two$Batch == "F2_B15")[1]] <- NA
  d_two$Var1[which(d_two$Batch == "F2_B04")[1]] <- NA
  expect_error(hotelling_classical_monitor(d_two, cal, vars),
               "Batch 'F2_B04' contains non-finite values")
})

test_that("the classical count can no longer come out NA either", {

  cal <- hotelling_classical_calibrate(afm_phase1, vars)
  u   <- hotelling_classical_ucl(K = 30, I = 20, J = 4)$UCL

  mon <- hotelling_classical_monitor(afm_phase2, cal, vars, ucl = u)
  expect_false(anyNA(mon$T2))
  expect_false(anyNA(mon$is_ooc))
  expect_false(is.na(sum(mon$is_ooc)))

  # El mismo agujero que la gemela robusta, y por la misma via: un NA en
  # Fase 2 convertia sum(is_ooc) entero en NA. Ahora aborta.
  d <- afm_phase2
  d$Var1[1] <- NA
  expect_error(hotelling_classical_monitor(d, cal, vars, ucl = u),
               "contains non-finite values")
  expect_error(hotelling_classical_monitor(d, cal, vars),
               "contains non-finite values")

  d_chr <- afm_phase2; d_chr$Var2 <- as.character(d_chr$Var2)
  expect_error(hotelling_classical_monitor(d_chr, cal, vars, ucl = u),
               "not numeric: Var2")
})

test_that("the classical block of summary() reports a number, not NA", {
  # El camino real: run_afm_mcd -> hotelling_classical_monitor -> summary.
  # Es donde el NA se imprimia al usuario.
  s <- summary(run_afm_mcd(afm_phase1, afm_phase2, vars, plot = FALSE,
                           compare_classical = TRUE))
  expect_false(is.na(s$classical$n_ooc))
  expect_type(s$classical$n_ooc, "integer")
})

test_that("run_afm_mcd delegates the classical flag instead of computing it", {

  study <- run_afm_mcd(afm_phase1, afm_phase2, vars, plot = FALSE,
                       compare_classical = TRUE)

  # La forma antigua: el envoltorio calculaba is_ooc a mano. Delegarla en
  # hotelling_classical_monitor() debe devolver el mismo objeto, columna a
  # columna y en el mismo orden.
  cal_c <- hotelling_classical_calibrate(afm_phase1, vars)
  ucl_c <- hotelling_classical_ucl(K = cal_c$n_batches, I = 20, J = 4,
                                   alpha = 0.001)
  mon_c <- hotelling_classical_monitor(afm_phase2, cal_c, vars)
  mon_c$is_ooc <- mon_c$T2 > ucl_c$UCL

  expect_identical(study$classical$monitoring, mon_c)
})
