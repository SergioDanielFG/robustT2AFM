#' Plot a Hotelling T-squared Control Chart
#'
#' Produces a publication-quality control chart from the output of
#' \code{\link{monitor_afm_mcd}} or \code{\link{hotelling_classical_monitor}},
#' overlaying the Upper Control Limit (UCL) as a horizontal reference line.
#' In-control batches are drawn in cyan and out-of-control batches in
#' burgundy red, with the connecting line segments colored by the status of
#' their endpoints. A status legend appears on the right margin and the UCL
#' is labeled inline with its numerical value.
#'
#' @param monitor_result A data frame returned by
#'   \code{\link{monitor_afm_mcd}} or
#'   \code{\link{hotelling_classical_monitor}}. Must contain columns
#'   \code{Batch} and \code{T2}.
#' @param UCL Numeric scalar. The Upper Control Limit, typically the
#'   \code{$UCL} element returned by \code{\link{ucl_F_adjusted}} (the
#'   operational default) or by \code{\link{hotelling_classical_ucl}} when
#'   producing a classical benchmark chart.
#' @param method_label Character. Method name shown in the plot title (e.g.
#'   "Classical Hotelling T2 Control Chart - Phase 2", "AFM-MCD Control
#'   Chart - Phase 2"). Default: "Control chart".
#' @param alpha Optional numeric in (0, 1). Nominal false-alarm rate used to
#'   compute the UCL. If provided, it is appended to the UCL inline label
#'   as "UCL = X (alpha = Y)". Default \code{NULL} (not displayed).
#' @param y_max Optional numeric. Upper limit of the y-axis. Default
#'   \code{NULL} = auto-scale to \code{max(1.20 * UCL, 1.18 * max(T2))} so
#'   that OOC batch labels have head-room.
#' @param show_values Logical. If \code{TRUE}, prints the T-squared numerical
#'   value above every batch. Default \code{FALSE}: with more than a handful of
#'   batches those labels collide with each other and with the out-of-control
#'   identifiers, and the exact value of an in-control batch supports no
#'   decision. The numbers are all in the \code{T2} column of the data frame
#'   you passed in.
#' @param save_path Optional character. Path stem (without extension) to
#'   which the chart is exported as both PNG (300 dpi) and PDF at
#'   publication dimensions (7 x 4.5 in). Default \code{NULL} (no export).
#'
#' @return A \code{ggplot} object. The plot is also drawn on the active
#'   graphics device.
#'
#' @details
#' Each batch is classified as either "In control" (T-squared <= UCL, shown in
#' cyan) or "Out of control" (T-squared > UCL, shown in burgundy), and the
#' legend is titled "Chart verdict" to say what the colours encode. That pair
#' of terms is the chart's verdict, and it is used identically in
#' \code{\link{plot_method_comparison}}, in \code{summary.afm_mcd_study} and
#' in the console output. It is deliberately distinct from the pair that names
#' what a batch \emph{is} - "Faulty batch" and "Fault-free batch" - which
#' appears only where ground truth was supplied by the caller.
#'
#' \strong{On the colour of the limit line.} It is burgundy here and in
#' \code{\link{plot_method_comparison}} because in both figures the line is a
#' limit: crossing it raises an alarm. In \code{\link{plot_afm_weights}} the
#' reference line is drawn in slate instead, because 1/K is not a limit and
#' crossing it means nothing. The two colours say which kind of line you are
#' looking at; it is not an oversight.
#'
#' This binary classification is determined exclusively by the UCL: the
#' function makes no use of contamination labels or ground truth, since in
#' production monitoring this information is not available. For Monte Carlo
#' studies that require ground-truth comparisons (true positives, false
#' alarms, missed detections), use the raw output of \code{monitor_*} and
#' compute confusion matrices directly rather than relying on the plot.
#'
#' That restriction belongs to this function, not to the package. When the
#' faulty batches genuinely are known - a simulation you generated, or a
#' labelled benchmark such as Tennessee Eastman -
#' \code{\link{plot_method_comparison}} will colour the points by truth if you
#' hand it that truth through its \code{faulty} argument, which is how the
#' published comparison figure is read: faulty batches sitting below a limit
#' that never fired. That is a validation device. It has no place in
#' production monitoring, where which batches are faulty is precisely the
#' question the chart is being asked.
#'
#' For paired comparisons of two methods on the same batches (classical
#' Hotelling against AFM-MCD), use \code{\link{plot_method_comparison}}, which
#' puts both charts on one figure with each batch in the same position in the
#' two panels. Calling this function twice also works and is what the vignette
#' still does when it introduces the two methods, but then the reader has to
#' match batches across two separate images by eye.
#'
#' @references
#' Montgomery, D. C. (2009). \emph{Introduction to Statistical Quality
#' Control}, 6th edition. John Wiley & Sons, Hoboken, NJ. ISBN 978-0-470-16992-6. Chapter 11.
#'
#' Wickham, H. (2016). \emph{ggplot2: Elegant Graphics for Data Analysis}
#' (2nd ed.). Springer.
#'
#' @importFrom ggplot2 ggplot aes geom_hline geom_line geom_point geom_text
#'   geom_segment scale_color_manual labs theme_minimal theme element_text
#'   element_blank element_line annotate scale_y_continuous scale_x_continuous
#'   margin guide_legend guides ggsave
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # End-to-end pipeline a quality engineer would run for a new lot:
#' # simulate historical and monitoring data, calibrate, compute the
#' # operational UCL, monitor Phase 2 and produce a labelled control chart.
#' data(afm_phase1)
#' data(afm_phase2)
#' vars <- paste0("Var", 1:4)
#' cal  <- calibrate_afm_mcd(afm_phase1, vars)
#' ucl  <- ucl_F_adjusted(cal, I = 20)
#' mon  <- monitor_afm_mcd(afm_phase2, cal, vars)
#'
#' plot_control_chart(
#'   mon,
#'   UCL          = ucl$UCL,
#'   method_label = "AFM-MCD Control Chart - Phase 2",
#'   alpha        = ucl$parameters$alpha
#' )
#'
#' # Save the chart at publication dimensions (PNG 300 dpi + PDF):
#' # plot_control_chart(mon, UCL = ucl$UCL, alpha = 0.001,
#' #                    save_path = file.path(tempdir(), "phase2_chart"))
plot_control_chart <- function(monitor_result,
                               UCL,
                               method_label = "Control chart",
                               alpha = NULL,
                               y_max = NULL,
                               show_values = FALSE,
                               save_path = NULL) {

  # --- Input validation ---
  if (!is.data.frame(monitor_result)) {
    stop("'monitor_result' must be a data frame.")
  }
  required_cols <- c("Batch", "T2")
  missing_cols <- setdiff(required_cols, colnames(monitor_result))
  if (length(missing_cols) > 0) {
    stop("'monitor_result' is missing required column(s): ",
         paste(missing_cols, collapse = ", "), ".")
  }
  if (nrow(monitor_result) < 1) {
    stop("'monitor_result' has zero rows; nothing to plot.")
  }
  if (!is.numeric(UCL) || length(UCL) != 1 || !is.finite(UCL) || UCL <= 0) {
    stop("'UCL' must be a single positive finite number.")
  }
  if (!is.character(method_label) || length(method_label) != 1) {
    stop("'method_label' must be a single character string.")
  }
  if (!is.null(y_max)) {
    if (!is.numeric(y_max) || length(y_max) != 1 || y_max <= 0) {
      stop("'y_max' must be a single positive number, or NULL.")
    }
  }
  if (!is.logical(show_values) || length(show_values) != 1) {
    stop("'show_values' must be TRUE or FALSE.")
  }
  if (!is.null(alpha)) {
    if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
      stop("'alpha' must be a single numeric in (0, 1), or NULL.")
    }
  }
  if (!is.null(save_path)) {
    if (!is.character(save_path) || length(save_path) != 1 ||
        nchar(save_path) == 0) {
      stop("'save_path' must be a single non-empty character string, or NULL.")
    }
  }

  # --- Build plot data with batch ordering preserved ---
  plot_df <- data.frame(
    Batch = factor(as.character(monitor_result$Batch),
                   levels = as.character(monitor_result$Batch)),
    BatchIdx = seq_len(nrow(monitor_result)),
    T2 = monitor_result$T2,
    stringsAsFactors = FALSE
  )

  # --- Classify each batch as In control or Out of control ---
  # Vocabulario del veredicto de la carta, comun a las tres figuras y a la
  # consola. "In control", no "under control", que es un calco: la literatura
  # de control estadistico dice in-control / out-of-control.
  plot_df$Status <- ifelse(plot_df$T2 > UCL, "Out of control", "In control")
  plot_df$Status <- factor(plot_df$Status,
                           levels = c("In control", "Out of control"))

  # --- Color palette ---
  color_map <- c("In control" = "#3FA9B6",     # cyan/turquoise
                 "Out of control" = "#A02D31")  # burgundy red

  # --- Compute y-axis upper bound ---
  # Solo hace falta hueco para las etiquetas de lote OOC, que van escalonadas
  # en dos alturas. Sin las etiquetas de valor, el 20% anterior era espacio
  # muerto.
  if (is.null(y_max)) {
    head_room <- if (isTRUE(show_values)) 1.20 else 1.12
    y_max <- max(head_room * UCL, head_room * max(plot_df$T2, na.rm = TRUE))
  }

  # --- Summary of OOC batches (for annotation) ---
  n_ooc     <- sum(plot_df$Status == "Out of control")
  n_total   <- nrow(plot_df)
  pct_ooc   <- if (n_total > 0) 100 * n_ooc / n_total else 0
  # "OOC" era jerga interna que ninguna leyenda definia.
  summary_label <- paste0(n_ooc, " out of control / ", n_total,
                          " (", format(round(pct_ooc, 1), nsmall = 1), "%)")

  ucl_label <- if (is.null(alpha)) {
    paste0("UCL = ", round(UCL, 2))
  } else {
    paste0("UCL = ", round(UCL, 2),
           "  (alpha = ", format(alpha, scientific = FALSE), ")")
  }

  # --- Segment colors: each segment colored "Out of control" if either
  #     endpoint is out of control, otherwise "In control" ---
  if (nrow(plot_df) >= 2) {
    seg_df <- data.frame(
      x_start = plot_df$BatchIdx[-nrow(plot_df)],
      x_end   = plot_df$BatchIdx[-1],
      y_start = plot_df$T2[-nrow(plot_df)],
      y_end   = plot_df$T2[-1],
      seg_color = ifelse(
        plot_df$Status[-nrow(plot_df)] == "Out of control" |
          plot_df$Status[-1]            == "Out of control",
        "Out of control", "In control"
      ),
      stringsAsFactors = FALSE
    )
    seg_df$seg_color <- factor(seg_df$seg_color,
                               levels = c("In control", "Out of control"))
  } else {
    seg_df <- NULL
  }

  # --- Build the ggplot object ---
  p <- ggplot2::ggplot()

  # Connecting line segments (colored by status)
  if (!is.null(seg_df)) {
    p <- p + ggplot2::geom_segment(
      data = seg_df,
      mapping = ggplot2::aes(x = .data$x_start, y = .data$y_start,
                             xend = .data$x_end, yend = .data$y_end,
                             color = .data$seg_color),
      linewidth = 0.6,
      show.legend = FALSE
    )
  }

  p <- p +

    # UCL reference line
    ggplot2::geom_hline(yintercept = UCL,
                        color = "#A02D31",
                        linetype = "dashed",
                        linewidth = 0.6) +

    # Points (colored by status)
    ggplot2::geom_point(
      data = plot_df,
      mapping = ggplot2::aes(x = .data$BatchIdx, y = .data$T2,
                             color = .data$Status),
      size = 2.8
    ) +

    # Color scale shared by points and segments
    ggplot2::scale_color_manual(values = color_map, drop = FALSE,
                                name = "Chart verdict") +

    # UCL label inline at the right margin, BELOW the line: the out-of-control
    # identifiers always sit above their points, so the strip just under the
    # limit at the right edge is the one place reliably free.
    ggplot2::annotate("text",
                      x = max(plot_df$BatchIdx) + 0.9,
                      y = UCL,
                      label = ucl_label,
                      hjust = 1.0, vjust = 1.5,
                      color = "#A02D31", size = 3.5,
                      fontface = "plain") +

    # OOC summary annotation, top-left inside the panel
    ggplot2::annotate("text",
                      x = min(plot_df$BatchIdx),
                      y = y_max * 0.97,
                      label = summary_label,
                      hjust = 0, vjust = 1,
                      color = if (n_ooc > 0) "#A02D31" else "#3FA9B6",
                      size = 3.6, fontface = "bold") +

    # Axes and title
    ggplot2::labs(
      title = method_label,
      x = "Batch",
      y = expression(T^2 ~ "statistic")
    ) +

    ggplot2::scale_y_continuous(limits = c(0, y_max),
                                expand = c(0, 0)) +
    # Margen a ambos lados: sin el, las etiquetas de lote de los extremos se
    # recortan contra el borde del panel.
    ggplot2::scale_x_continuous(
      breaks = plot_df$BatchIdx,
      labels = as.character(plot_df$Batch),
      expand = ggplot2::expansion(add = 0.9)
    ) +

    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 13,
                                               hjust = 0.5,
                                               margin = ggplot2::margin(b = 10)),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      axis.title       = ggplot2::element_text(face = "plain"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey90",
                                                 linewidth = 0.3),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_text(face = "bold", size = 10),
      legend.text      = ggplot2::element_text(size = 9),
      legend.margin    = ggplot2::margin(t = 4),
      plot.margin      = ggplot2::margin(t = 8, r = 12, b = 8, l = 8)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3)))

  # --- Optional: T-squared values above each point ---
  if (isTRUE(show_values)) {
    p <- p + ggplot2::geom_text(
      data = plot_df,
      mapping = ggplot2::aes(x = .data$BatchIdx, y = .data$T2,
                             label = format(round(.data$T2, 1), nsmall = 1)),
      vjust = -1.0, size = 2.8, color = "grey25",
      show.legend = FALSE
    )
  }

  # --- Bold red batch-ID labels above each OOC point ---
  # Escalonadas: dos lotes OOC contiguos en el eje x tienen etiquetas que se
  # tocan, y en una racha larga se funden. Se alterna la altura solo dentro
  # de la racha; al abrirse un hueco de mas de 2 posiciones se vuelve abajo.
  ooc_df <- plot_df[plot_df$Status == "Out of control", , drop = FALSE]
  if (nrow(ooc_df) > 0) {
    ooc_df <- ooc_df[order(ooc_df$BatchIdx), , drop = FALSE]
    level <- integer(nrow(ooc_df))
    level[1] <- 1L
    if (nrow(ooc_df) > 1) {
      for (i in 2:nrow(ooc_df)) {
        near <- (ooc_df$BatchIdx[i] - ooc_df$BatchIdx[i - 1]) <= 2
        level[i] <- if (near) 3L - level[i - 1] else 1L
      }
    }
    base_v <- if (isTRUE(show_values)) -2.8 else -1.2
    ooc_df$vjust <- ifelse(level == 1L, base_v, base_v - 1.3)

    p <- p + ggplot2::geom_text(
      data = ooc_df,
      mapping = ggplot2::aes(x = .data$BatchIdx, y = .data$T2,
                             label = as.character(.data$Batch),
                             vjust = .data$vjust),
      size = 3.0, color = "#A02D31", fontface = "bold",
      show.legend = FALSE
    )
  }

  # --- Optional export to PNG + PDF at publication dimensions ---
  if (!is.null(save_path)) {
    ggplot2::ggsave(paste0(save_path, ".png"), plot = p,
                    width = 7, height = 4.5, units = "in", dpi = 300)
    ggplot2::ggsave(paste0(save_path, ".pdf"), plot = p,
                    width = 7, height = 4.5, units = "in")
  }

  return(p)
}
