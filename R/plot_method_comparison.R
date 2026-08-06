#' Compare the Robust and Classical Charts on the Same Phase 2 Batches
#'
#' Draws both control charts side by side over the same Phase 2 batches in the
#' same order, each against its own control limit. By default the figure is
#' bare: two panels, a labelled limit line and a legend, with no explanatory
#' text of any kind. Diagnostic marks are available on request.
#'
#' @param robust Either an \code{afm_mcd_study} returned by
#'   \code{\link{run_afm_mcd}} with \code{compare_classical = TRUE}, or the
#'   data frame returned by \code{\link{monitor_afm_mcd}}. In the second form
#'   the next three arguments are required.
#' @param classical The data frame returned by
#'   \code{\link{hotelling_classical_monitor}} on the same batches. Ignored
#'   when \code{robust} is a study.
#' @param UCL_robust Numeric scalar, the robust limit, normally
#'   \code{ucl_F_adjusted(cal, I)$UCL}. Ignored when \code{robust} is a study.
#' @param UCL_classical Numeric scalar, the classical limit, normally
#'   \code{hotelling_classical_ucl(K, I, J)$UCL}. Ignored when \code{robust}
#'   is a study.
#' @param faulty Optional character vector with the identifiers of the batches
#'   known to be faulty. Only supply it when that is genuinely known: a
#'   simulation you generated, or a labelled benchmark such as Tennessee
#'   Eastman. Default \code{NULL}. See Details.
#' @param colour_by One of \code{"auto"} (default), \code{"truth"} or
#'   \code{"decision"}. \code{"auto"} colours by truth when \code{faulty} is
#'   supplied and by the chart's own verdict otherwise. \code{"truth"} without
#'   \code{faulty} is an error, since there is no truth to colour by.
#' @param diagnostics Logical. If \code{TRUE}, adds the marks that compare the
#'   two verdicts: a hollow ring on every batch signalled by one method and not
#'   the other, drawn in the panel that stayed silent, a dotted vertical guide
#'   through both panels at those batches, and a count in the corner of each
#'   panel. Default \code{FALSE}, which leaves the figure bare.
#' @param near_miss Numeric in (0, 1), or \code{NULL} (default). When set,
#'   shades the band of that width immediately below the classical limit, as a
#'   fraction of it, in the classical panel only. See Details before reading
#'   anything into it.
#' @param scale Either \code{"ratio"} (default) to plot each statistic divided
#'   by its own method's limit, or \code{"T2"} to plot the raw statistics with
#'   independent vertical scales. See Details.
#' @param labels Character vector of length 2 with the panel titles. Default
#'   \code{c("AFM-MCD (robust)", "Classical Hotelling")}.
#' @param title Character. Title of the whole figure, or \code{NULL} for none.
#'   Default "Same Phase 2 batches, two methods". A figure exported to PNG or
#'   PDF and pasted into a report carries no journal caption, so it has to say
#'   what it is on its own.
#' @param save_path Optional character. Path stem (without extension) for
#'   export as PNG (300 dpi) and PDF. Default \code{NULL} writes nothing.
#'
#' @return A \code{ggplot} object with two panels.
#'
#' @details
#' \strong{What this figure is about.} When Phase 1 contains contaminated
#' batches, masking inflates the classical pooled covariance. A larger
#' covariance produces a wider control limit, and a wider limit signals less
#' often. The classical chart does not become nervous under contamination, it
#' becomes deaf, and the failure mode to illustrate is lost detection, not
#' false alarms. That sentence is the whole point of the figure and it is
#' deliberately not printed under it: a caption nobody reads is worse than no
#' caption, because the lines that do matter get skipped with it.
#'
#' \strong{Colouring by truth versus by verdict.} With \code{faulty} the points
#' are coloured "Faulty batch" and "Fault-free batch" under a legend titled
#' "Batch condition", and the reader sees faulty batches sitting below a limit
#' that never fired. Without it the colours can only show each chart's own
#' verdict, "Out of control" and "In control", under a legend titled "Chart
#' verdict". The two pairs share no word and never appear together, and the
#' legend title names which of the two the colours encode: a reader looking at
#' this figure next to \code{\link{plot_control_chart}} should never have to
#' wonder whether the same colour means the same thing.
#' Which batches are genuinely faulty is knowable in a simulation, because we
#' generate it, and in a labelled benchmark such as Tennessee Eastman, because
#' the data carries it. It is never knowable on a plant floor: that is the very
#' thing the chart is there to find out. Colouring by truth therefore requires
#' the caller to supply the truth, and is a validation device, not a production
#' one. See \code{\link{plot_control_chart}}, which never colours by truth.
#'
#' \strong{Disagreement versus missed detection.} The diagnostic marks report
#' which batches one method signals and the other does not, which is computable
#' anywhere. Only when \code{faulty} is supplied do the counts become
#' detections and the wording speak of batches missed.
#'
#' \strong{Why \code{scale = "ratio"} is the default.} Raw statistics force
#' independent vertical axes, because the robust and classical T-squared values
#' are of different magnitudes, and independent axes make the panel heights
#' meaningless to compare. That needs a warning, and this figure has no caption
#' to put one in. Dividing each statistic by its own limit removes the need for
#' the warning instead of restating it: both panels then share one scale whose
#' reference line is 1.0, the axis says so, and heights are directly
#' comparable. \code{scale = "T2"} reproduces the look of the published figure
#' and is the one case where a caption line appears, because there the figure
#' has to defend itself against its own natural reading.
#'
#' \strong{The near-miss band is a display device, not a verdict.} It marks
#' batches whose classical statistic falls within \code{near_miss} of the limit
#' from below, which is how masking shows up: the signal is not absent, it is
#' pushed under a limit that moved. The threshold is an argument because it is
#' a choice. A statistic at 85 percent of its limit was not "almost detected":
#' the relation between T-squared and the probability of signalling is not
#' linear, and the package makes no such claim.
#'
#' @references
#' Frutos-Galarza, S. D., Ruiz-Barzola, O., Ramirez, J., &
#' Galindo-Villardon, P. (2026). A robust Hotelling-type T2 control chart
#' combining the minimum covariance determinant estimator with multiple
#' factor analysis weighting. Under review.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline geom_rect
#'   geom_text facet_wrap vars labs scale_color_manual scale_x_continuous
#'   scale_y_continuous expansion theme_minimal theme element_text
#'   element_blank element_line margin ggsave
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @export
#'
#' @examples
#' data(afm_phase1)
#' data(afm_phase2)
#' vars <- paste0("Var", 1:4)
#'
#' # Form 1: from a study that already carries the classical baseline.
#' study <- run_afm_mcd(afm_phase1, afm_phase2, vars, plot = FALSE,
#'                      compare_classical = TRUE)
#' plot_method_comparison(study)
#'
#' # Form 2: from the pieces, which is how the paper's own scripts chain.
#' cal   <- calibrate_afm_mcd(afm_phase1, vars)
#' ucl   <- ucl_F_adjusted(cal, I = 20)
#' mon   <- monitor_afm_mcd(afm_phase2, cal, vars, ucl = ucl$UCL)
#' cal_c <- hotelling_classical_calibrate(afm_phase1, vars)
#' ucl_c <- hotelling_classical_ucl(K = 30, I = 20, J = 4)
#' mon_c <- hotelling_classical_monitor(afm_phase2, cal_c, vars)
#'
#' plot_method_comparison(mon, mon_c, ucl$UCL, ucl_c$UCL)
#'
#' # These are simulated data, so the faulty batches are known and the
#' # points can be coloured by truth: faulty batches sitting below a limit
#' # that never fired. Real monitoring data carries no such column.
#' truth <- as.character(unique(
#'   afm_phase2$Batch[afm_phase2$Status == "Faulty"]
#' ))
#' plot_method_comparison(study, faulty = truth)
#'
#' # Everything the bare figure leaves out, for a closer look.
#' plot_method_comparison(study, faulty = truth,
#'                        diagnostics = TRUE, near_miss = 0.20)
plot_method_comparison <- function(robust,
                                   classical = NULL,
                                   UCL_robust = NULL,
                                   UCL_classical = NULL,
                                   faulty = NULL,
                                   colour_by = c("auto", "truth", "decision"),
                                   diagnostics = FALSE,
                                   near_miss = NULL,
                                   scale = c("ratio", "T2"),
                                   labels = c("AFM-MCD (robust)",
                                              "Classical Hotelling"),
                                   title = "Same Phase 2 batches, two methods",
                                   save_path = NULL) {

  scale     <- match.arg(scale)
  colour_by <- match.arg(colour_by)

  # --- Resolve the two input forms ---
  if (inherits(robust, "afm_mcd_study")) {
    if (is.null(robust$classical)) {
      stop("This study has no classical baseline to compare against. ",
           "Re-run it with compare_classical = TRUE:\n",
           "  run_afm_mcd(phase1, phase2, variables, compare_classical = TRUE)")
    }
    mon_r <- robust$monitoring
    mon_c <- robust$classical$monitoring
    ucl_r <- robust$ucl$UCL
    ucl_c <- robust$classical$ucl$UCL
  } else {
    if (!is.data.frame(robust)) {
      stop("'robust' must be an afm_mcd_study from run_afm_mcd(), or the ",
           "data frame from monitor_afm_mcd().")
    }
    if (!is.data.frame(classical)) {
      stop("'classical' must be the data frame from ",
           "hotelling_classical_monitor(). It is required when 'robust' is ",
           "a data frame rather than a study object.")
    }
    if (!is.numeric(UCL_robust) || length(UCL_robust) != 1 ||
        !is.finite(UCL_robust) || UCL_robust <= 0) {
      stop("'UCL_robust' must be a single positive finite number, ",
           "normally ucl_F_adjusted(calibration, I)$UCL.")
    }
    if (!is.numeric(UCL_classical) || length(UCL_classical) != 1 ||
        !is.finite(UCL_classical) || UCL_classical <= 0) {
      stop("'UCL_classical' must be a single positive finite number, ",
           "normally hotelling_classical_ucl(K, I, J)$UCL.")
    }
    mon_r <- robust
    mon_c <- classical
    ucl_r <- UCL_robust
    ucl_c <- UCL_classical
  }

  for (nm in c("Batch", "T2")) {
    if (!nm %in% colnames(mon_r) || !nm %in% colnames(mon_c)) {
      stop("Both monitoring results must contain columns 'Batch' and 'T2'.")
    }
  }
  if (nrow(mon_r) < 1) stop("The monitoring results have zero rows.")

  b_r <- as.character(mon_r$Batch)
  b_c <- as.character(mon_c$Batch)
  if (!setequal(b_r, b_c)) {
    stop("The two monitoring results cover different batches, so they ",
         "cannot be compared batch by batch. Run both methods on the same ",
         "Phase 2 data.")
  }
  # Se reordena el clasico al orden del robusto: cada lote ocupa la misma
  # posicion en los dos paneles.
  mon_c <- mon_c[match(b_r, b_c), , drop = FALSE]

  # --- Remaining validation ---
  if (!is.logical(diagnostics) || length(diagnostics) != 1 ||
      is.na(diagnostics)) {
    stop("'diagnostics' must be TRUE or FALSE.")
  }
  if (!is.null(near_miss)) {
    if (!is.numeric(near_miss) || length(near_miss) != 1 ||
        !is.finite(near_miss) || near_miss <= 0 || near_miss >= 1) {
      stop("'near_miss' must be a single number in (0, 1) - the width of the ",
           "band below the classical limit, as a fraction of it - or NULL to ",
           "omit the band.")
    }
  }
  if (!is.character(labels) || length(labels) != 2 || anyDuplicated(labels)) {
    stop("'labels' must be two distinct character strings.")
  }
  if (!is.null(title) && (!is.character(title) || length(title) != 1)) {
    stop("'title' must be a single character string, or NULL for no title.")
  }
  if (!is.null(faulty)) {
    faulty <- as.character(faulty)
    unknown <- setdiff(faulty, b_r)
    if (length(unknown) > 0) {
      stop("These 'faulty' batches are not among the monitored batches: ",
           paste(unknown, collapse = ", "),
           ". 'faulty' must list identifiers that appear in the Phase 2 data.")
    }
  }
  if (colour_by == "truth" && is.null(faulty)) {
    stop("colour_by = \"truth\" needs the truth: pass 'faulty' with the ",
         "identifiers of the batches known to be faulty. That is knowable in ",
         "a simulation or a labelled benchmark, never in production, which is ",
         "why the package will not guess it.")
  }
  if (!is.null(save_path)) {
    if (!is.character(save_path) || length(save_path) != 1 ||
        nchar(save_path) == 0) {
      stop("'save_path' must be a single non-empty character string, or NULL.")
    }
  }

  # --- Assemble the plotting frame ---
  K <- length(b_r)
  flag_r <- mon_r$T2 > ucl_r
  flag_c <- mon_c$T2 > ucl_c

  plot_df <- data.frame(
    Method   = factor(rep(labels, each = K), levels = labels),
    Batch    = rep(b_r, 2),
    BatchIdx = rep(seq_len(K), 2),
    T2       = c(mon_r$T2, mon_c$T2),
    UCL      = rep(c(ucl_r, ucl_c), each = K),
    flagged  = c(flag_r, flag_c),
    stringsAsFactors = FALSE
  )
  plot_df$y <- if (scale == "ratio") plot_df$T2 / plot_df$UCL else plot_df$T2

  # --- Colour: truth when it was supplied, verdict otherwise ---
  by_truth <- colour_by == "truth" ||
    (colour_by == "auto" && !is.null(faulty))
  # Dos vocabularios, uno por idea, sin ninguna palabra en comun, y el titulo
  # de la leyenda dice cual de los dos se esta mirando.
  if (by_truth) {
    lv <- c("Fault-free batch", "Faulty batch")
    legend_title <- "Batch condition"
    plot_df$Status <- factor(
      ifelse(plot_df$Batch %in% faulty, lv[2], lv[1]), levels = lv
    )
  } else {
    lv <- c("In control", "Out of control")
    legend_title <- "Chart verdict"
    plot_df$Status <- factor(ifelse(plot_df$flagged, lv[2], lv[1]), levels = lv)
  }
  color_map <- stats::setNames(c("#3FA9B6", "#A02D31"), lv)

  # --- Reference line, one per panel ---
  hline_df <- data.frame(
    Method = factor(labels, levels = labels),
    y      = if (scale == "ratio") c(1, 1) else c(ucl_r, ucl_c),
    label  = if (scale == "ratio") {
      c("own limit", "own limit")
    } else {
      c(paste0("UCL = ", round(ucl_r, 2)), paste0("UCL = ", round(ucl_c, 2)))
    },
    stringsAsFactors = FALSE
  )

  # --- Build the plot: bare unless something is asked for ---
  p <- ggplot2::ggplot()

  if (!is.null(near_miss)) {
    lo <- (1 - near_miss) * ucl_c
    band_df <- data.frame(
      Method = factor(labels[2], levels = labels),
      ymin   = if (scale == "ratio") 1 - near_miss else lo,
      ymax   = if (scale == "ratio") 1 else ucl_c,
      stringsAsFactors = FALSE
    )
    band_df$label <- sprintf("within %.0f%% of the limit", 100 * near_miss)
    p <- p +
      ggplot2::geom_rect(
        data = band_df,
        mapping = ggplot2::aes(xmin = -Inf, xmax = Inf,
                               ymin = .data$ymin, ymax = .data$ymax),
        fill = "#A02D31", alpha = 0.08, inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data = band_df,
        mapping = ggplot2::aes(x = Inf,
                               y = (.data$ymin + .data$ymax) / 2,
                               label = .data$label),
        hjust = 1.03, vjust = 0.5, size = 2.6, color = "#A02D31",
        inherit.aes = FALSE
      )
  }

  if (isTRUE(diagnostics)) {
    only_r <- b_r[flag_r & !flag_c]
    only_c <- b_r[flag_c & !flag_r]
    disagree <- c(only_r, only_c)

    if (length(disagree) > 0) {
      p <- p + ggplot2::geom_vline(
        xintercept = which(b_r %in% disagree),
        color = "#A02D31", linetype = "dotted", linewidth = 0.4
      )
    }
  }

  p <- p +
    ggplot2::geom_hline(
      data = hline_df,
      mapping = ggplot2::aes(yintercept = .data$y),
      color = "#A02D31", linetype = "dashed", linewidth = 0.6
    ) +
    ggplot2::geom_point(
      data = plot_df,
      mapping = ggplot2::aes(x = .data$BatchIdx, y = .data$y,
                             color = .data$Status),
      size = 2.4
    ) +
    ggplot2::scale_color_manual(values = color_map, drop = FALSE,
                                name = legend_title)

  if (isTRUE(diagnostics)) {
    only_r <- b_r[flag_r & !flag_c]
    only_c <- b_r[flag_c & !flag_r]

    # El anillo va en el panel donde el lote NO suena.
    ring_df <- rbind(
      plot_df[plot_df$Method == labels[2] & plot_df$Batch %in% only_r, ],
      plot_df[plot_df$Method == labels[1] & plot_df$Batch %in% only_c, ]
    )
    if (nrow(ring_df) > 0) {
      p <- p + ggplot2::geom_point(
        data = ring_df,
        mapping = ggplot2::aes(x = .data$BatchIdx, y = .data$y),
        shape = 21, size = 4.4, stroke = 0.9,
        color = "#A02D31", fill = NA, inherit.aes = FALSE
      )
    }

    count_label <- function(flag) {
      if (is.null(faulty)) {
        sprintf("signalled %d of %d", sum(flag), K)
      } else {
        is_f <- b_r %in% faulty
        sprintf("detected %d of %d faulty, %d false alarm(s)",
                sum(flag & is_f), length(faulty), sum(flag & !is_f))
      }
    }
    ann_df <- data.frame(
      Method = factor(labels, levels = labels),
      label  = c(count_label(flag_r), count_label(flag_c)),
      stringsAsFactors = FALSE
    )
    p <- p + ggplot2::geom_text(
      data = ann_df,
      mapping = ggplot2::aes(x = -Inf, y = Inf, label = .data$label),
      hjust = -0.04, vjust = 1.6, size = 3, fontface = "bold",
      color = "grey25", inherit.aes = FALSE
    )
  }

  # Etiquetas del eje x: identificadores reales, adelgazados si son muchos.
  # Nunca posiciones: quien ve F2_B01, F2_B03 ubica el resto contando; quien
  # ve 1, 3 pierde la trazabilidad a su lote.
  brk <- seq_len(K)
  if (K > 15) brk <- brk[seq(1L, K, by = 2L)]

  # El pie existe solo cuando la figura necesita defenderse de su lectura
  # natural, es decir con escalas libres.
  cap <- if (scale == "T2") {
    paste0("Panels use independent vertical scales: heights are NOT ",
           "comparable between panels.")
  } else {
    NULL
  }

  p <- p +
    ggplot2::geom_text(
      data = hline_df,
      mapping = ggplot2::aes(x = Inf, y = .data$y, label = .data$label),
      hjust = 1.03, vjust = -0.6, size = 3, color = "#A02D31",
      inherit.aes = FALSE
    ) +
    ggplot2::scale_x_continuous(breaks = brk, labels = b_r[brk],
                                expand = c(0.04, 0.04)) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.15))
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$Method), ncol = 2,
                        scales = if (scale == "ratio") "fixed" else "free_y") +
    ggplot2::labs(
      title = title,
      x = "Batch",
      y = if (scale == "ratio") {
        "T2 / own control limit"
      } else {
        expression(T^2 ~ "statistic")
      },
      caption = cap
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.caption       = ggplot2::element_text(size = 7.5, hjust = 0,
                                                 color = "grey35",
                                                 margin = ggplot2::margin(t = 8)),
      strip.text         = ggplot2::element_text(face = "bold", size = 10),
      axis.text.x        = ggplot2::element_text(angle = 45, hjust = 1,
                                                 size = 7),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey90",
                                                 linewidth = 0.3),
      legend.position    = "bottom",
      plot.margin        = ggplot2::margin(t = 8, r = 12, b = 8, l = 8)
    )

  if (!is.null(save_path)) {
    ggplot2::ggsave(paste0(save_path, ".png"), plot = p,
                    width = 9, height = 4.5, units = "in", dpi = 300)
    ggplot2::ggsave(paste0(save_path, ".pdf"), plot = p,
                    width = 9, height = 4.5, units = "in")
  }

  return(p)
}
