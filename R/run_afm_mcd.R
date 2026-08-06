#' Run the Complete AFM-MCD Study in One Call
#'
#' Chains the four steps of the method - Phase 1 calibration, operational
#' control limit, Phase 2 monitoring and control chart - and returns a single
#' object that knows how to describe itself. This is the entry point for
#' someone who wants a control chart without first learning which function
#' takes which argument; the individual steps remain available separately as
#' \code{\link{calibrate_afm_mcd}}, \code{\link{ucl_F_adjusted}},
#' \code{\link{monitor_afm_mcd}} and \code{\link{plot_control_chart}}.
#'
#' @param phase1 A data frame of historical, in-control batches used for
#'   calibration. Must contain a column identifying each batch; see
#'   \code{batch_col}.
#' @param phase2 A data frame of the new batches to monitor, with the same
#'   batch column.
#' @param variables Character vector with the names of the process variables.
#'   Default \code{NULL} auto-detects every numeric column of \code{phase1}
#'   other than 'Batch' and 'Phase', and reports the choice with a
#'   \code{message()}. Read that message: identifier or index columns that
#'   happen to be numeric would be picked up as if they were process
#'   variables, which yields plausible and wrong results. Name the variables
#'   explicitly whenever the data frame carries anything but process
#'   measurements.
#' @param mcd_alpha Numeric in (0.60, 0.90). Proportion of observations kept
#'   by MCD within each Phase 1 batch. Default 0.67, as in
#'   Frutos-Galarza et al. (2026).
#' @param alpha Numeric in (0, 1). Nominal false-alarm rate for the control
#'   limit. Default 0.001.
#' @param plot Logical. If \code{TRUE} (default) the control chart is drawn on
#'   the active graphics device. Both plots are built and stored either way.
#' @param compare_classical Logical. If \code{TRUE}, the classical
#'   (non-robust) Hotelling chart is computed on the same data and reported
#'   alongside in \code{summary()}. Default \code{FALSE}.
#' @param save_path Optional character. Path stem (without extension) for
#'   exporting the figures. When supplied, four files are written:
#'   \code{<stem>_control_chart.png/.pdf} and
#'   \code{<stem>_afm_weights.png/.pdf}. Default \code{NULL} writes nothing.
#' @param batch_col Character. Name of the column identifying the batch in
#'   both phases. Default \code{"Batch"}, and propagated to every step. It is
#'   also excluded from the automatic variable detection, together with
#'   \code{"Phase"}.
#' @param x An \code{afm_mcd_study} object (or, for the summary method, a
#'   \code{summary.afm_mcd_study} object) to print.
#' @param object An \code{afm_mcd_study} object to summarise.
#' @param ... Not used; present for compatibility with the generics.
#'
#' @return \code{run_afm_mcd} returns an object of class
#'   \code{afm_mcd_study}, a list with components:
#' \describe{
#'   \item{calibration}{The list returned by \code{\link{calibrate_afm_mcd}}.}
#'   \item{ucl}{The list returned by \code{\link{ucl_F_adjusted}}.}
#'   \item{monitoring}{The data frame returned by
#'         \code{\link{monitor_afm_mcd}}, including the \code{is_ooc} column.}
#'   \item{chart}{The \code{ggplot} control chart.}
#'   \item{weight_plot}{The \code{ggplot} of the Phase 1 AFM weights.}
#'   \item{classical}{\code{NULL}, or a list with the classical
#'         \code{calibration}, \code{ucl} and \code{monitoring} when
#'         \code{compare_classical = TRUE}.}
#'   \item{variables}{The process variables actually used.}
#'   \item{call}{The matched call that produced the object.}
#' }
#'   \code{summary} returns an object of class
#'   \code{summary.afm_mcd_study} holding the reported quantities, so they can
#'   be reused programmatically. The \code{print} methods return their
#'   argument invisibly and are called for their console output.
#'
#' @details
#' Nothing is recomputed here: the wrapper calls the same four functions in
#' order and stores their results untouched, so any figure it reports is the
#' figure those functions produce. The number of observations per batch is not
#' an argument because \code{\link{calibrate_afm_mcd}} already records it as
#' \code{I_phase1}, and that is what is passed to the control limit.
#'
#' The classical comparison, when requested, prints the two sets of numbers
#' side by side and draws no conclusion. Which method does better depends on
#' how contaminated the calibration set is, and that is not something the
#' package can determine from the data alone.
#'
#' @section Reading the example, and what it does not show:
#' The example below runs the base configuration of Frutos-Galarza et al.
#' (2026) on one seed. With that seed the proposed method flags 8 of the 10
#' off-target batches and the classical chart flags 2, neither raising a false
#' alarm on the 10 in-control batches.
#'
#' Those are the counts of a single replication, not detection rates. The
#' paper reports 0.863 for the proposed method and 0.270 for the classical one
#' averaged over 2000 replications, and it measures them differently: there
#' both methods are first equalised to a common in-control ARL0 using the
#' empirical quantile of each method's own null distribution over 5000
#' in-control batches, whereas this example applies each method's operational
#' limit. The two figures describe similar but not identical quantities and
#' should not be read as one estimating the other.
#'
#' The example also illustrates only half of the paper's argument. With a
#' shift of 1 sigma and alpha = 0.001, what separates the two methods here is
#' detection power: across twelve seeds and 200 replications neither method
#' produced a single false alarm. Nothing in this example speaks to false
#' alarm behaviour.
#'
#' @references
#' Frutos-Galarza, S. D., Ruiz-Barzola, O., Ramirez, J., &
#' Galindo-Villardon, P. (2026). A robust Hotelling-type T2 control chart
#' combining the minimum covariance determinant estimator with multiple
#' factor analysis weighting. Under review.
#'
#' @export
#'
#' @examples
#' # The whole method in one call. The variables are detected automatically
#' # because the four Var columns are the only numeric ones.
#' data(afm_phase1)
#' data(afm_phase2)
#' study <- run_afm_mcd(afm_phase1, afm_phase2, plot = FALSE)
#'
#' study                      # the short answer
#' summary(study)             # the full report
#' study$monitoring           # the T2 of every Phase 2 batch
#'
#' # The same run with the classical chart alongside. On this data the
#' # proposed method flags 8 of the 10 off-target batches and the classical
#' # chart flags 2, with no false alarm from either. One replication is not
#' # a detection rate: see the section above before quoting these counts.
#' summary(run_afm_mcd(afm_phase1, afm_phase2, plot = FALSE,
#'                     compare_classical = TRUE))
run_afm_mcd <- function(phase1,
                        phase2,
                        variables = NULL,
                        mcd_alpha = 0.67,
                        alpha = 0.001,
                        plot = TRUE,
                        compare_classical = FALSE,
                        save_path = NULL,
                        batch_col = "Batch") {

  this_call <- match.call()

  # --- Input validation ---
  if (!is.data.frame(phase1)) {
    stop("'phase1' must be a data frame of historical batches.")
  }
  if (!is.data.frame(phase2)) {
    stop("'phase2' must be a data frame of the batches to monitor.")
  }
  check_batch_col(phase1, batch_col, "phase1",
                  "run_afm_mcd(phase1, phase2, batch_col = \"<name>\")")
  check_batch_col(phase2, batch_col, "phase2",
                  "run_afm_mcd(phase1, phase2, batch_col = \"<name>\")")
  if (!is.logical(plot) || length(plot) != 1 || is.na(plot)) {
    stop("'plot' must be TRUE or FALSE.")
  }
  if (!is.logical(compare_classical) || length(compare_classical) != 1 ||
      is.na(compare_classical)) {
    stop("'compare_classical' must be TRUE or FALSE.")
  }
  if (!is.null(save_path)) {
    if (!is.character(save_path) || length(save_path) != 1 ||
        nchar(save_path) == 0) {
      stop("'save_path' must be a single non-empty character string, or NULL.")
    }
  }

  # --- Variable selection ---
  if (is.null(variables)) {
    candidates <- setdiff(colnames(phase1), c(batch_col, "Phase"))
    variables <- candidates[vapply(phase1[candidates], is.numeric, logical(1))]
    if (length(variables) < 1) {
      stop("No numeric process variables found in 'phase1' besides 'Batch' ",
           "and 'Phase'. Pass 'variables' explicitly with the column names ",
           "to monitor.")
    }
    # Siempre se avisa: una columna identificadora numerica (run, sample,
    # faultNumber) entraria como variable de proceso y daria un resultado
    # verosimil y falso. El aviso no depende de que el usuario imprima nada.
    message("Auto-detected ", length(variables), " process variable(s): ",
            paste(variables, collapse = ", "),
            "\n  Check this list. Numeric identifier or index columns ",
            "(run number, sample number, fault code) would be monitored as ",
            "if they were measurements, producing plausible but wrong ",
            "results. Pass 'variables' explicitly to choose them yourself.")
  }
  check_variables(phase1, variables, "phase1",
                  "run_afm_mcd(phase1, phase2, variables), or NULL to detect them")
  missing_p2 <- setdiff(variables, colnames(phase2))
  if (length(missing_p2) > 0) {
    stop("Variable(s) not found in 'phase2': ",
         paste(missing_p2, collapse = ", "),
         ". They are present in 'phase1', so the two data frames do not carry ",
         "the same process variables. Check the Phase 2 export.")
  }

  # --- The four steps, in order, with nothing recomputed ---
  calibration <- calibrate_afm_mcd(phase1, variables, mcd_alpha = mcd_alpha,
                                   batch_col = batch_col)
  limit       <- ucl_F_adjusted(calibration, I = calibration$I_phase1,
                                alpha = alpha)
  monitoring  <- monitor_afm_mcd(phase2, calibration, variables,
                                 ucl = limit$UCL, batch_col = batch_col)

  chart <- plot_control_chart(
    monitoring,
    UCL          = limit$UCL,
    alpha        = alpha,
    method_label = "AFM-MCD Control Chart - Phase 2",
    save_path    = if (is.null(save_path)) NULL else {
      paste0(save_path, "_control_chart")
    }
  )
  weight_plot <- plot_afm_weights(
    calibration,
    save_path = if (is.null(save_path)) NULL else {
      paste0(save_path, "_afm_weights")
    }
  )

  # --- Optional classical baseline on the same data ---
  classical <- NULL
  if (isTRUE(compare_classical)) {
    cal_c <- hotelling_classical_calibrate(phase1, variables,
                                           batch_col = batch_col)
    ucl_c <- hotelling_classical_ucl(K = cal_c$n_batches,
                                     I = calibration$I_phase1,
                                     J = length(variables),
                                     alpha = alpha)
    mon_c <- hotelling_classical_monitor(phase2, cal_c, variables,
                                         batch_col = batch_col)
    mon_c$is_ooc <- mon_c$T2 > ucl_c$UCL
    classical <- list(calibration = cal_c, ucl = ucl_c, monitoring = mon_c)
  }

  if (isTRUE(plot)) {
    print(chart)
  }

  structure(
    list(
      calibration = calibration,
      ucl         = limit,
      monitoring  = monitoring,
      chart       = chart,
      weight_plot = weight_plot,
      classical   = classical,
      variables   = variables,
      call        = this_call
    ),
    class = "afm_mcd_study"
  )
}


#' @rdname run_afm_mcd
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


#' @rdname run_afm_mcd
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
      det_ratio = det(object$classical$calibration$Sp) / det(cal$Sw)
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


#' @rdname run_afm_mcd
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
  cat(sprintf(
    "  Each batch is summarised from the least-scattered %.0f%% of its own\n",
    100 * x$mcd_alpha))
  cat("  observations, which limits how far stray points inside a batch can\n")
  cat("  distort the reference covariance. It does not protect the reference\n")
  cat("  centre: a batch shifted bodily off-target keeps its internal shape,\n")
  cat("  so no setting of that percentage makes it stand out here.\n")
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
  cat("  A low weight means the batch is more spread out than the rest.\n")
  cat("  It does not mean the batch is off-target, and it does not mean\n")
  cat(sprintf("  the batch is contaminated.  $weight_plot shows all %d.\n", x$K1))

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
    cat(sprintf("  det(Sp) / det(Sw) = %.2f\n", x$classical$det_ratio))
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
