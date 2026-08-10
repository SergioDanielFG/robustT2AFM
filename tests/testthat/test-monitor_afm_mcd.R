test_that("monitor_afm_mcd returns non-negative T^2 for every batch", {
  set.seed(20260417)
  sim  <- simulate_batch_process(K1 = 12, K2 = 5, I = 20, J = 4,
                                 prop_ooc_F2 = 0.4, shift_ooc = 2,
                                 seed = 20260417)
  vars <- paste0("Var", 1:4)
  cal  <- suppressMessages(
    calibrate_afm_mcd(subset(sim, Phase == "Phase 1"), vars)
  )
  mon  <- monitor_afm_mcd(subset(sim, Phase == "Phase 2"), cal, vars)

  expect_s3_class(mon, "data.frame")
  expect_named(mon, c("Batch", "I", "T2"))
  expect_true(all(mon$T2 >= 0))
  expect_true(all(mon$I == 20L))
})

test_that("monitor_afm_mcd yields T^2 = 0 when batch mean equals mu_r", {
  # Trivial calibration: mu_r = 0, Sw = I; batch with zero mean must return 0.
  cal <- list(
    mu_r    = c(0, 0, 0, 0),
    Sw      = diag(4),
    weights = setNames(rep(1 / 30, 30), paste0("B", 1:30))
  )
  new_batch <- data.frame(
    Batch = "F2_new",
    Var1  = c(-1, 1),
    Var2  = c(-1, 1),
    Var3  = c(-1, 1),
    Var4  = c(-1, 1)
  )
  mon <- monitor_afm_mcd(new_batch, cal, paste0("Var", 1:4))
  expect_equal(mon$T2, 0)
  expect_equal(mon$I, 2L)
})

test_that("monitor_afm_mcd adds is_ooc only when ucl is supplied", {
  sim  <- simulate_batch_process(K1 = 12, K2 = 5, I = 20, J = 4,
                                 prop_ooc_F2 = 0.4, shift_ooc = 2,
                                 seed = 20260417)
  vars <- paste0("Var", 1:4)
  cal  <- calibrate_afm_mcd(subset(sim, Phase == "Phase 1"), vars)
  ph2  <- subset(sim, Phase == "Phase 2")

  mon_plain <- monitor_afm_mcd(ph2, cal, vars)
  mon_flag  <- monitor_afm_mcd(ph2, cal, vars,
                               ucl = ucl_F_adjusted(cal, I = 20)$UCL)

  # Default output is untouched: three columns, not one more.
  expect_named(mon_plain, c("Batch", "I", "T2"))

  # ucl = NULL reproduces the default output exactly.
  expect_identical(monitor_afm_mcd(ph2, cal, vars, ucl = NULL), mon_plain)

  # The flag is additive: the first three columns are unchanged.
  expect_named(mon_flag, c("Batch", "I", "T2", "is_ooc"))
  expect_identical(mon_flag[c("Batch", "I", "T2")], mon_plain)
  expect_type(mon_flag$is_ooc, "logical")
})

test_that("monitor_afm_mcd flags strictly above the limit, never on it", {
  # mu_r = 0, Sw = I, batch of 2 observations with x_bar = (2, 2, 2, 2):
  # T2 = I * ||x_bar||^2 = 2 * 16 = 32, an exact value with no rounding.
  cal  <- list(mu_r = c(0, 0, 0, 0), Sw = diag(4))
  vars <- paste0("Var", 1:4)
  new_batch <- data.frame(
    Batch = "F2_new",
    Var1  = c(1, 3), Var2 = c(1, 3),
    Var3  = c(1, 3), Var4 = c(1, 3)
  )
  expect_equal(monitor_afm_mcd(new_batch, cal, vars)$T2, 32)

  expect_false(monitor_afm_mcd(new_batch, cal, vars, ucl = 33)$is_ooc)
  expect_true(monitor_afm_mcd(new_batch, cal, vars, ucl = 31)$is_ooc)

  # Exactly on the limit: NOT flagged (strictly greater, not >=).
  expect_false(monitor_afm_mcd(new_batch, cal, vars, ucl = 32)$is_ooc)
})

# Reimplementacion literal del bucle con rbind() que habia antes de
# vectorizar. Se conserva aqui para poder comparar objeto contra objeto: los
# T2 se calculan igual, asi que la salida debe ser identica columna por
# columna, tipos y row.names incluidos.
old_rbind_loop <- function(new_data, calibration, variables,
                           batch_col = "Batch") {
  mu_r     <- calibration$mu_r
  Sw_inv   <- solve(calibration$Sw)
  batch_id <- new_data[[batch_col]]
  batches  <- unique(batch_id)

  results <- data.frame(Batch = character(), I = integer(), T2 = numeric(),
                        stringsAsFactors = FALSE)
  for (batch in batches) {
    subset_batch <- new_data[batch_id == batch, variables, drop = FALSE]
    I     <- nrow(subset_batch)
    x_bar <- colMeans(subset_batch)
    diff  <- x_bar - mu_r
    T2    <- as.numeric(I * t(diff) %*% Sw_inv %*% diff)
    results <- rbind(results, data.frame(Batch = as.character(batch), I = I,
                                         T2 = T2, stringsAsFactors = FALSE))
  }
  results
}

test_that("vectorising the loop returns an object identical to the old rbind", {
  vars <- paste0("Var", 1:4)
  cal  <- calibrate_afm_mcd(afm_phase1, vars)

  expect_identical(monitor_afm_mcd(afm_phase2, cal, vars),
                   old_rbind_loop(afm_phase2, cal, vars))

  # Y con el limite puesto, la columna is_ooc se apoya en los mismos T2.
  u   <- ucl_F_adjusted(cal, I = 20)$UCL
  new <- monitor_afm_mcd(afm_phase2, cal, vars, ucl = u)
  old <- old_rbind_loop(afm_phase2, cal, vars)
  old$is_ooc <- old$T2 > u
  expect_identical(new, old)

  # Tambien con lotes de tamano desigual y en orden no alfabetico, que es
  # donde una vectorizacion descuidada reordena las filas.
  shuffled <- afm_phase2[afm_phase2$Batch %in% c("F2_B07", "F2_B02", "F2_B11"), ]
  shuffled <- shuffled[order(match(shuffled$Batch,
                                   c("F2_B07", "F2_B02", "F2_B11"))), ]
  shuffled <- shuffled[-c(1, 2), ]   # F2_B07 se queda con 18 observaciones
  expect_identical(monitor_afm_mcd(shuffled, cal, vars),
                   old_rbind_loop(shuffled, cal, vars))
  expect_identical(monitor_afm_mcd(shuffled, cal, vars)$Batch,
                   c("F2_B07", "F2_B02", "F2_B11"))
})

test_that("an invalid batch still aborts at the same point in the loop", {
  vars <- paste0("Var", 1:4)
  cal  <- calibrate_afm_mcd(afm_phase1, vars)

  # El ultimo lote: el error solo aparece si el bucle llega hasta el.
  d_last <- afm_phase2
  d_last$Var1[which(d_last$Batch == "F2_B20")[1]] <- NA
  expect_error(monitor_afm_mcd(d_last, cal, vars),
               "Batch 'F2_B20' contains non-finite values")

  # Con dos lotes invalidos se nombra el primero en orden de aparicion, no el
  # primero por fila ni el ultimo. Si la guarda se sacara del bucle a una
  # pasada previa sobre todo el data frame, este es el test que lo detecta.
  d_two <- afm_phase2
  d_two$Var1[which(d_two$Batch == "F2_B15")[1]] <- NA
  d_two$Var1[which(d_two$Batch == "F2_B04")[1]] <- NA
  expect_error(monitor_afm_mcd(d_two, cal, vars),
               "Batch 'F2_B04' contains non-finite values")

  # Y si el orden de aparicion cambia, cambia el lote nombrado: se recorre por
  # aparicion, no por nombre.
  d_rev <- d_two[order(match(d_two$Batch,
                             c("F2_B15", "F2_B04",
                               setdiff(levels(d_two$Batch),
                                       c("F2_B15", "F2_B04"))))), ]
  expect_error(monitor_afm_mcd(d_rev, cal, vars),
               "Batch 'F2_B15' contains non-finite values")
})

test_that("sum(mon$is_ooc) can no longer come out NA", {
  vars <- paste0("Var", 1:4)
  cal  <- calibrate_afm_mcd(afm_phase1, vars)
  u    <- ucl_F_adjusted(cal, I = 20)$UCL

  # 1. En el camino normal el conteo es un numero, nunca NA.
  mon <- monitor_afm_mcd(afm_phase2, cal, vars, ucl = u)
  expect_false(anyNA(mon$T2))
  expect_false(anyNA(mon$is_ooc))
  expect_false(is.na(sum(mon$is_ooc)))
  expect_identical(sum(mon$is_ooc), 8L)

  # 2. La entrada que antes producia el NA ahora aborta, asi que no queda
  #    ningun camino por el que monitor_afm_mcd() devuelva un is_ooc con NA.
  #    Antes de la guarda esto daba T2 = NA, is_ooc = NA y suma = NA.
  d <- afm_phase2
  d$Var1[1] <- NA
  expect_error(monitor_afm_mcd(d, cal, vars, ucl = u),
               "contains non-finite values")

  # 3. Y lo mismo para NaN e Inf, que propagan igual.
  d_nan <- afm_phase2; d_nan$Var2[3] <- NaN
  expect_error(monitor_afm_mcd(d_nan, cal, vars, ucl = u),
               "contains non-finite values")
  d_inf <- afm_phase2; d_inf$Var3[7] <- Inf
  expect_error(monitor_afm_mcd(d_inf, cal, vars, ucl = u),
               "contains non-finite values")

  # 4. La guarda no depende de que se pase el limite: sin ucl tambien aborta,
  #    porque el NA ya esta en T2 antes de que exista is_ooc.
  expect_error(monitor_afm_mcd(d, cal, vars), "contains non-finite values")

  # 5. Una columna de texto en Fase 2 da un mensaje que nombra la variable,
  #    en vez del "'x' debe ser numerico" de colMeans().
  d_chr <- afm_phase2; d_chr$Var2 <- as.character(d_chr$Var2)
  expect_error(monitor_afm_mcd(d_chr, cal, vars, ucl = u),
               "not numeric: Var2")
})

test_that("monitor_afm_mcd validates ucl", {
  cal  <- list(mu_r = c(0, 0, 0, 0), Sw = diag(4))
  vars <- paste0("Var", 1:4)
  new_batch <- data.frame(
    Batch = "F2_new",
    Var1  = c(-1, 1), Var2 = c(-1, 1),
    Var3  = c(-1, 1), Var4 = c(-1, 1)
  )

  expect_error(monitor_afm_mcd(new_batch, cal, vars, ucl = -1),
               "single positive, finite number")
  expect_error(monitor_afm_mcd(new_batch, cal, vars, ucl = 0),
               "single positive, finite number")
  expect_error(monitor_afm_mcd(new_batch, cal, vars, ucl = c(10, 20)),
               "single positive, finite number")
  expect_error(monitor_afm_mcd(new_batch, cal, vars, ucl = "19.7"),
               "single positive, finite number")
  expect_error(monitor_afm_mcd(new_batch, cal, vars, ucl = NA_real_),
               "single positive, finite number")

  # Passing the whole ucl_F_adjusted() object is the predictable mistake:
  # the error must name the fix.
  expect_error(monitor_afm_mcd(new_batch, cal, vars,
                               ucl = list(UCL = 19.7, method = "F-adjusted")),
               "ucl_F_adjusted\\(calibration, I\\)\\$UCL")
})

test_that("monitor_afm_mcd reports a singular Sw with a usable message", {
  # Var4 is an exact copy of Var3, so Sw is rank-deficient and not invertible.
  Sw <- matrix(c(1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 1, 1,
                 0, 0, 1, 1), nrow = 4, byrow = TRUE)
  cal <- list(mu_r = c(0, 0, 0, 0), Sw = Sw)
  new_batch <- data.frame(
    Batch = "F2_new",
    Var1  = c(-1, 1), Var2 = c(-1, 1),
    Var3  = c(-1, 1), Var4 = c(-1, 1)
  )

  expect_error(monitor_afm_mcd(new_batch, cal, paste0("Var", 1:4)),
               "singular and cannot be inverted")
  # The message must tell the engineer what to do, not only what failed.
  expect_error(monitor_afm_mcd(new_batch, cal, paste0("Var", 1:4)),
               "redundant")
})
