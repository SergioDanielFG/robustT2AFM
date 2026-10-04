# Semilla interna del MCD de los centros (version 0.3.0). Un valor fijo y
# documentado: no se ajusto a ningun resultado.
CENTER_SEED <- 20260927L

# Ejecuta `expr` con una semilla fija y devuelve el generador al estado en que
# estaba, exista o no una semilla previa en la sesion del usuario.
with_local_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  set.seed(seed)
  expr
}
