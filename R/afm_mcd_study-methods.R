#' Print and Summarise an AFM-MCD Study
#'
#' @param x An \code{afm_mcd_study} object (or, for the summary method, a
#'   \code{summary.afm_mcd_study} object) to print.
#' @param object An \code{afm_mcd_study} object to summarise.
#' @param ... Not used; present for compatibility with the generics.
#'
#' @return \code{summary} returns an object of class
#'   \code{summary.afm_mcd_study} holding the reported quantities, so they can
#'   be reused programmatically. The \code{print} methods return their
#'   argument invisibly and are called for their console output.
#'
#' @seealso \code{\link{run_afm_mcd}}
#'
#' @name afm_mcd_study-methods
NULL


#' @rdname afm_mcd_study-methods
#' @export
print.afm_mcd_study <- function(x, ...) {

  mon    <- x$monitoring
  n      <- nrow(mon)
  n_ooc  <- sum(mon$is_ooc)
  pct    <- if (n > 0) 100 * n_ooc / n else 0
  K1     <- length(x$calibration$weights)
  top    <- which.max(mon$T2)

  cat("AFM-MCD robust T2 study\n")
  cat(strrep("-", 62), "\n", sep = "")
  cat(sprintf("Variables (%d):   %s\n", length(x$variables),
              paste(x$variables, collapse = ", ")))
  cat(sprintf("Phase 1:         %d calibration batches, %d observations each\n",
              K1, x$calibration$I_phase1))
  cat(sprintf("Phase 2:         %d monitored batches\n", n))
  cat(sprintf("Control limit:   UCL = %.2f  (alpha = %s)\n\n",
              x$ucl$UCL, format(x$ucl$parameters$alpha, scientific = FALSE)))

  if (n_ooc == 0) {
    cat(sprintf("No Phase 2 batch is out of control (0 of %d).\n", n))
    cat(sprintf("  Highest: %s, T2 = %.2f (%.2fx the limit)\n",
                mon$Batch[top], mon$T2[top], mon$T2[top] / x$ucl$UCL))
  } else {
    cat(sprintf("%d of %d Phase 2 batches are out of control (%.1f%%):\n",
                n_ooc, n, pct))
    ooc_names <- mon$Batch[mon$is_ooc]
    shown     <- ooc_names[seq_len(min(10L, length(ooc_names)))]
    cat("  ", paste(shown, collapse = ", "),
        if (length(ooc_names) > length(shown)) {
          paste0(", ... and ", length(ooc_names) - length(shown), " more")
        } else "", "\n", sep = "")
    # Dos decimales, como la columna "x limit" del summary: es el mismo
    # cociente y no puede salir con distinta precision en dos sitios.
    cat(sprintf("  Highest: %s, T2 = %.2f (%.2fx the limit)\n",
                mon$Batch[top], mon$T2[top], mon$T2[top] / x$ucl$UCL))
  }

  cat("\nsummary() for the full report.  $chart for the control chart.\n")
  invisible(x)
}


#' @rdname afm_mcd_study-methods
#' @export
summary.afm_mcd_study <- function(object, ...) {

  mon   <- object$monitoring
  cal   <- object$calibration
  w     <- cal$weights
  K1    <- length(w)
  n_ooc <- sum(mon$is_ooc)

  ooc_tab <- mon[mon$is_ooc, c("Batch", "T2"), drop = FALSE]
  if (nrow(ooc_tab) > 0) {
    ooc_tab       <- ooc_tab[order(-ooc_tab$T2), , drop = FALSE]
    ooc_tab$ratio <- ooc_tab$T2 / object$ucl$UCL
    rownames(ooc_tab) <- NULL
  }

  # Lote mas alto que sigue bajo control: el candidato a vigilar.
  in_ctrl <- mon[!mon$is_ooc, c("Batch", "T2"), drop = FALSE]
  closest <- if (nrow(in_ctrl) > 0) {
    row <- in_ctrl[which.max(in_ctrl$T2), , drop = FALSE]
    list(Batch = row$Batch, T2 = row$T2, ratio = row$T2 / object$ucl$UCL)
  } else {
    NULL
  }

  classical <- if (is.null(object$classical)) {
    NULL
  } else {
    list(
      UCL   = object$classical$ucl$UCL,
      n_ooc = sum(object$classical$monitoring$is_ooc),
      det_ratio   = det(object$classical$calibration$Sp) / det(cal$Sw),
      trace_ratio = sum(diag(object$classical$calibration$Sp)) /
        sum(diag(cal$Sw)),
      centre_gap  = sqrt(sum((object$classical$calibration$mu_global -
                                cal$mu_r)^2))
    )
  }

  structure(
    list(
      call        = object$call,
      variables   = object$variables,
      K1          = K1,
      I1          = cal$I_phase1,
      K2          = nrow(mon),
      I2          = if (length(unique(mon$I)) == 1) mon$I[1] else NA_integer_,
      mcd_alpha   = cal$mcd_alpha,
      mu_r        = cal$mu_r,
      w_min       = min(w),
      w_min_batch = names(w)[which.min(w)],
      w_max       = max(w),
      w_max_batch = names(w)[which.max(w)],
      uniform     = 1 / K1,
      lowest      = sort(w)[seq_len(min(3L, K1))],
      ucl         = object$ucl,
      n_ooc       = n_ooc,
      pct_ooc     = if (nrow(mon) > 0) 100 * n_ooc / nrow(mon) else 0,
      ooc_table   = ooc_tab,
      closest     = closest,
      classical   = classical
    ),
    class = "summary.afm_mcd_study"
  )
}


#' @rdname afm_mcd_study-methods
#' @export
print.summary.afm_mcd_study <- function(x, ...) {

  cat("AFM-MCD robust T2 study\n")
  cat(strrep("=", 62), "\n", sep = "")
  # deparse() parte la llamada larga en varias lineas sangradas; al unirlas
  # hay que comprimir los espacios o quedan huecos en medio.
  call_txt <- gsub("\\s+", " ", paste(deparse(x$call), collapse = " "))
  cat("Call:  ", call_txt, "\n\n", sep = "")

  cat("DATA\n")
  cat(sprintf("  Variables (%d):     %s\n", length(x$variables),
              paste(x$variables, collapse = ", ")))
  cat(sprintf("  Phase 1:           %d batches, %d observations each\n",
              x$K1, x$I1))
  cat(sprintf("  Phase 2:           %d batches%s\n", x$K2,
              if (is.na(x$I2)) "" else sprintf(", %d observations each", x$I2)))

  cat("\nPHASE 1 CALIBRATION\n")

  cat("  Reference centre:  ",
      paste(format(round(x$mu_r, 3), nsmall = 3), collapse = "  "), "\n", sep = "")
  cat(sprintf("  AFM weights:       lowest  %.4f  (%s)\n",
              x$w_min, x$w_min_batch))
  cat(sprintf("                     highest %.4f  (%s)\n",
              x$w_max, x$w_max_batch))
  cat(sprintf("                     even weighting would be 1/%d = %.4f\n",
              x$K1, x$uniform))
  cat("  Most internally dispersed batches:\n")
  for (i in seq_along(x$lowest)) {
    cat(sprintf("                     %-8s %.4f\n",
                names(x$lowest)[i], x$lowest[i]))
  }
  cat("  Weight ranks internal dispersion, not contamination: see\n")
  cat(sprintf("  ?plot_afm_weights.  $weight_plot shows all %d.\n", x$K1))

  cat("\nCONTROL LIMIT\n")
  cat(sprintf("  UCL = %.4f      %s, alpha = %s\n",
              x$ucl$UCL, x$ucl$method,
              format(x$ucl$parameters$alpha, scientific = FALSE)))
  cat(sprintf("                     (J = %d, K = %d, m* = %d, df2 = %d)\n",
              x$ucl$parameters$J, x$ucl$parameters$K,
              x$ucl$parameters$m_star, x$ucl$parameters$df2))

  cat(sprintf("\nPHASE 2 MONITORING            %d of %d out of control (%.1f%%)\n",
              x$n_ooc, x$K2, x$pct_ooc))
  if (x$n_ooc > 0) {
    cat(sprintf("  %-8s %8s %9s\n", "Batch", "T2", "x limit"))
    shown <- x$ooc_table[seq_len(min(15L, nrow(x$ooc_table))), , drop = FALSE]
    for (i in seq_len(nrow(shown))) {
      cat(sprintf("  %-8s %8.2f %9.2f\n",
                  shown$Batch[i], shown$T2[i], shown$ratio[i]))
    }
    if (nrow(x$ooc_table) > nrow(shown)) {
      cat(sprintf("  ... and %d more, see $monitoring\n",
                  nrow(x$ooc_table) - nrow(shown)))
    }
  }
  if (!is.null(x$closest)) {
    cat(sprintf("  Closest batch still in control:  %s, T2 = %.2f (%.2fx)\n",
                x$closest$Batch, x$closest$T2, x$closest$ratio))
  }
  cat("  Full table in $monitoring.\n")

  if (!is.null(x$classical)) {
    cat("\nCLASSICAL BASELINE (Hotelling T2, non-robust)\n")
    cat(sprintf("  Classical:  UCL = %8.4f    %d of %d out of control\n",
                x$classical$UCL, x$classical$n_ooc, x$K2))
    cat(sprintf("  AFM-MCD:    UCL = %8.4f    %d of %d out of control\n",
                x$ucl$UCL, x$n_ooc, x$K2))
    cat(sprintf("  trace(Sp) / trace(Sw) = %.2f     det ratio = %.2f\n",
                x$classical$trace_ratio, x$classical$det_ratio))
    cat(sprintf("  reference centres %.3f apart\n", x$classical$centre_gap))
  }

  cat("\nNEXT STEP\n")
  if (x$n_ooc > 0) {
    cat(sprintf(
      "  Pull the process records for the %d flagged batch(es) and check\n",
      x$n_ooc))
    cat("  what they have in common: materials, shift, operator, equipment.\n")
  } else {
    cat("  Nothing flagged. Keep monitoring, and recalibrate when the\n")
    cat("  process itself changes (new equipment, new formulation).\n")
  }

  invisible(x)
}
