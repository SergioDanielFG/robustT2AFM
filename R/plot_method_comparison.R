#' Compare the Robust and Classical Charts on the Same Phase 2 Batches
#'
#' Draws both control charts as two stacked panels over the same Phase 2
#' batches in the same order, each against its own control limit, so that a
#' batch is a single vertical line crossing both panels. Batches that one
#' method flags and the other does not are marked in both panels, which is
#' what the figure exists to show.
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
#'   Eastman. Default \code{NULL}, in which case the figure reports
#'   disagreement between the methods and never uses the word "missed". See
#'   Details.
#' @param near_miss Numeric in (0, 1), or \code{NULL}. Width, as a fraction of
#'   the classical limit, of a shaded band drawn immediately below that limit
#'   in the classical panel. Default 0.20. \code{NULL} removes the band.
#' @param scale Either \code{"ratio"} (default) to plot each statistic divided
#'   by its own method's limit, which puts both panels on a common scale with
#'   a single reference line at 1.0, or \code{"T2"} to plot the raw
#'   statistics with independent vertical scales. See Details.
#' @param labels Character vector of length 2 with the panel titles. Default
#'   \code{c("AFM-MCD (robust)", "Classical Hotelling")}.
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
#' becomes deaf: the failure mode to illustrate is lost detection, not false
#' alarms.
#'
#' \strong{Disagreement versus missed detection.} Without \code{faulty}, the
#' figure can only report which batches one method flags and the other does
#' not. That is computable anywhere, including in production, where nobody
#' knows which batches are truly faulty. Supplying \code{faulty} is a claim of
#' ground truth, and only then does the figure count detections and describe a
#' batch as missed. The distinction is kept in the wording: "flagged by one
#' method only" without it, "missed" with it.
#'
#' \strong{Why \code{scale = "ratio"} is the default.} With raw statistics the
#' two panels need independent vertical axes, because the robust and classical
#' T-squared values are of different magnitudes. Independent axes make the
#' panel heights meaningless to compare, and a reader who does not notice the
#' caveat will read the shorter classical bars as if they meant something.
#' Dividing each statistic by its own limit removes the problem: both panels
#' then share one scale whose reference line is 1.0, and heights are directly
#' comparable. It is the same device as the \code{x limit} column of
#' \code{summary.afm_mcd_study}. \code{scale = "T2"} remains available and
#' states the caveat in its own caption.
#'
#' \strong{The near-miss band is a display device, not a verdict.} The band
#' marks batches whose classical statistic falls within \code{near_miss} of
#' the limit from below, which is how masking shows up: the signal is not
#' absent, it is pushed under a limit that moved. The threshold is an argument
#' because it is a choice, its definition is written into the legend, and the
#' band is drawn only in the classical panel. The figure reports how many
#' batches fall in the band and nothing more. A statistic at 85 percent of its
#' limit was not "almost detected": the relation between T-squared and the
#' probability of signalling is not linear, and the package makes no such
#' claim.
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
#' @export
#'
#' @examples
#' # Base configuration of Frutos-Galarza et al. (2026).
#' sim <- simulate_batch_process(
#'   K1 = 30, K2 = 20, I = 20, J = 4, rho = 0.6,
#'   outlier_batches_F1 = 6, outlier_rate = 0.20, outlier_shift = 4,
#'   prop_ooc_F2 = 0.5, shift_ooc = 1.0,
#'   seed = 20260425
#' )
#' p1   <- subset(sim, Phase == "Phase 1")
#' p2   <- subset(sim, Phase == "Phase 2")
#' vars <- paste0("Var", 1:4)
#'
#' # Form 1: from a study that already carries the classical baseline.
#' study <- run_afm_mcd(p1, p2, vars, plot = FALSE, compare_classical = TRUE)
#' plot_method_comparison(study)
#'
#' # Form 2: from the pieces, which is how the paper's own scripts chain.
#' cal   <- calibrate_afm_mcd(p1, vars)
#' ucl   <- ucl_F_adjusted(cal, I = 20)
#' mon   <- monitor_afm_mcd(p2, cal, vars, ucl = ucl$UCL)
#' cal_c <- hotelling_classical_calibrate(p1, vars)
#' ucl_c <- hotelling_classical_ucl(K = 30, I = 20, J = 4)
#' mon_c <- hotelling_classical_monitor(p2, cal_c, vars)
#'
#' plot_method_comparison(mon, mon_c, ucl$UCL, ucl_c$UCL)
#'
#' # In a simulation the faulty batches are known, so the figure may
#' # count missed detections instead of mere disagreement.
#' truth <- as.character(unique(
#'   p2$Batch[p2$Status == "Out of Control"]
#' ))
#' plot_method_comparison(study, faulty = truth)
#'
#' # Raw statistics, with the scale caveat stated in the caption.
#' plot_method_comparison(study, scale = "T2")
plot_method_comparison <- function(robust,
                                   classical = NULL,
                                   UCL_robust = NULL,
                                   UCL_classical = NULL,
                                   faulty = NULL,
                                   near_miss = 0.20,
                                   scale = c("ratio", "T2"),
                                   labels = c("AFM-MCD (robust)",
                                              "Classical Hotelling"),
                                   save_path = NULL) {

  scale <- match.arg(scale)

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
  # Se reordena el clasico al orden del robusto: cada lote es una vertical.
  mon_c <- mon_c[match(b_r, b_c), , drop = FALSE]

  # --- Remaining validation ---
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
  if (!is.null(faulty)) {
    faulty <- as.character(faulty)
    unknown <- setdiff(faulty, b_r)
    if (length(unknown) > 0) {
      stop("These 'faulty' batches are not among the monitored batches: ",
           paste(unknown, collapse = ", "),
           ". 'faulty' must list identifiers that appear in the Phase 2 data.")
    }
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
  plot_df$Status <- factor(ifelse(plot_df$flagged, "Signalled", "No signal"),
                           levels = c("No signal", "Signalled"))

  # --- Disagreement: flagged by one method only, in either direction ---
  only_r <- b_r[flag_r & !flag_c]
  only_c <- b_r[flag_c & !flag_r]
  disagree <- c(only_r, only_c)

  # El anillo va en el panel donde el lote NO suena.
  ring_df <- rbind(
    plot_df[plot_df$Method == labels[2] & plot_df$Batch %in% only_r, ],
    plot_df[plot_df$Method == labels[1] & plot_df$Batch %in% only_c, ]
  )

  # --- Per-panel counts ---
  n_faulty <- if (is.null(faulty)) NA_integer_ else length(faulty)
  count_label <- function(flag) {
    if (is.null(faulty)) {
      sprintf("signalled %d of %d batches", sum(flag), K)
    } else {
      is_f <- b_r %in% faulty
      sprintf("detected %d of %d faulty  |  %d false alarm(s)",
              sum(flag & is_f), n_faulty, sum(flag & !is_f))
    }
  }
  ann_df <- data.frame(
    Method = factor(labels, levels = labels),
    label  = c(count_label(flag_r), count_label(flag_c)),
    stringsAsFactors = FALSE
  )

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

  # --- Near-miss band, classical panel only ---
  band_df <- NULL
  n_band  <- 0L
  if (!is.null(near_miss)) {
    lo <- (1 - near_miss) * ucl_c
    n_band <- sum(mon_c$T2 >= lo & mon_c$T2 <= ucl_c)
    band_df <- data.frame(
      Method = factor(labels[2], levels = labels),
      ymin   = if (scale == "ratio") 1 - near_miss else lo,
      ymax   = if (scale == "ratio") 1 else ucl_c,
      stringsAsFactors = FALSE
    )
  }

  # --- Caption: mechanism first, then what each device means ---
  cap <- paste0(
    "Masking inflates the Phase 1 covariance, which widens the limit; a ",
    "wider limit signals less often.\n",
    if (scale == "ratio") {
      paste0("Both panels share one scale: T2 divided by each method's own ",
             "limit, so 1.0 is that limit and heights are comparable.\n")
    } else {
      paste0("Panels use independent vertical scales: the heights are NOT ",
             "comparable between panels.\n")
    },
    if (length(disagree) > 0) {
      paste0("Ringed and marked by a vertical guide: ",
             if (is.null(faulty)) {
               paste0(length(disagree),
                      " batch(es) signalled by one method only.\n")
             } else {
               paste0(sum(only_r %in% faulty),
                      " faulty batch(es) missed by the classical chart and ",
                      "caught by the robust one.\n")
             })
    } else "",
    if (!is.null(near_miss)) {
      sprintf(paste0("Shaded band: within %.0f%% below the classical limit ",
                     "(%d batch(es)). This counts batches, it does not mean ",
                     "they were almost detected."),
              100 * near_miss, n_band)
    } else ""
  )

  color_map <- c("No signal" = "#3FA9B6", "Signalled" = "#A02D31")

  # --- Build the plot ---
  p <- ggplot2::ggplot()

  if (!is.null(band_df)) {
    p <- p + ggplot2::geom_rect(
      data = band_df,
      mapping = ggplot2::aes(xmin = -Inf, xmax = Inf,
                             ymin = .data$ymin, ymax = .data$ymax),
      fill = "#A02D31", alpha = 0.08, inherit.aes = FALSE
    )
  }

  if (length(disagree) > 0) {
    p <- p + ggplot2::geom_vline(
      xintercept = which(b_r %in% disagree),
      color = "#A02D31", linetype = "dotted", linewidth = 0.4
    )
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
                                name = NULL)

  # Anillo sin etiqueta: el eje x ya lleva el identificador y la guia vertical
  # baja hasta el. Etiquetar cada anillo produce solapes cuando los lotes en
  # desacuerdo son consecutivos.
  if (nrow(ring_df) > 0) {
    p <- p + ggplot2::geom_point(
      data = ring_df,
      mapping = ggplot2::aes(x = .data$BatchIdx, y = .data$y),
      shape = 21, size = 4.4, stroke = 0.9,
      color = "#A02D31", fill = NA, inherit.aes = FALSE
    )
  }

  p <- p +
    ggplot2::geom_text(
      data = ann_df,
      mapping = ggplot2::aes(x = -Inf, y = Inf, label = .data$label),
      hjust = -0.04, vjust = 1.6, size = 3.2, fontface = "bold",
      color = "grey25", inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = hline_df,
      mapping = ggplot2::aes(x = Inf, y = .data$y, label = .data$label),
      hjust = 1.03, vjust = -0.6, size = 3, color = "#A02D31",
      inherit.aes = FALSE
    ) +
    # Aire extra arriba: el recuento se dibuja en la esquina superior y si no
    # se amplia el rango choca con los puntos mas altos del panel.
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.20))
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$Method), ncol = 1,
                        scales = if (scale == "ratio") "fixed" else "free_y") +
    ggplot2::labs(
      title = "Same Phase 2 batches, two methods",
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
      plot.title         = ggplot2::element_text(face = "bold", size = 13,
                                                 hjust = 0.5,
                                                 margin = ggplot2::margin(b = 8)),
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

  # Etiquetas del eje x: identificadores, no indices.
  p <- p + ggplot2::scale_x_continuous(breaks = seq_len(K), labels = b_r,
                                       expand = c(0.03, 0.03))

  if (!is.null(save_path)) {
    ggplot2::ggsave(paste0(save_path, ".png"), plot = p,
                    width = 7, height = 6.5, units = "in", dpi = 300)
    ggplot2::ggsave(paste0(save_path, ".pdf"), plot = p,
                    width = 7, height = 6.5, units = "in")
  }

  return(p)
}
