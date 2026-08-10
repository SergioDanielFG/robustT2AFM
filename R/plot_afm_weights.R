#' Plot the AFM Weights of the Phase 1 Calibration Batches
#'
#' One horizontal bar per Phase 1 batch, sorted by AFM weight, smallest at
#' the top. The dashed line is the average weight, 1/K.
#'
#' @param calibration A list returned by \code{\link{calibrate_afm_mcd}}.
#'   Only its \code{weights} element is plotted; the whole object is taken so
#'   the function can check it really comes from a calibration.
#' @param title Character. Plot title. Default
#'   \code{"AFM weights - Phase 1 calibration"}.
#' @param highlight_lowest Integer >= 0. Number of lowest-weight batches to
#'   draw in burgundy instead of the neutral colour, to point at specific
#'   batches in a report. Default 0 (every bar neutral). This is a display
#'   choice made by the caller, not a verdict computed by the package.
#' @param save_path Optional character. Path stem (without extension) to
#'   which the chart is exported as both PNG (300 dpi) and PDF. Default
#'   \code{NULL}, which writes nothing to disk.
#'
#' @details
#' A low weight means the batch is more spread out than the others. It does
#' not mean the batch is off-target, and it does not mean it is contaminated:
#' a batch shifted bodily away from the rest keeps its weight. The weighting
#' protects the reference covariance \code{Sw}, not the reference centre
#' \code{mu_r}.
#'
#' More than half the batches fall below 1/K even when Phase 1 is clean,
#' because the weights average to 1/K by construction and their distribution
#' is skewed to the right. That is why the bars are not coloured by that
#' line. Use \code{highlight_lowest} to point at specific batches.
#'
#' @return A ggplot object.
#'
#' @references
#' Abdi, H., Williams, L. J., & Valentin, D. (2013). Multiple factor
#' analysis: principal component analysis for multitable and multiblock data
#' sets. \emph{Wiley Interdisciplinary Reviews: Computational Statistics},
#' 5(2), 149-179. \doi{10.1002/wics.1246}
#'
#' Wickham, H. (2016). \emph{ggplot2: Elegant Graphics for Data Analysis}
#' (2nd ed.). Springer.
#'
#' @examples
#' data(afm_phase1)
#' cal <- calibrate_afm_mcd(afm_phase1, paste0("Var", 1:4))
#'
#' # Every bar neutral; the dashed line is the average weight 1/K.
#' plot_afm_weights(cal)
#'
#' # Point at the six most dispersed batches for a report. Four of those six
#' # carry outliers and two are clean, and two contaminated batches sit
#' # outside the six: the weight ranks dispersion, it does not classify.
#' plot_afm_weights(cal, highlight_lowest = 6)
#'
#' # Save at publication dimensions (PNG 300 dpi + PDF):
#' # plot_afm_weights(cal, save_path = file.path(tempdir(), "afm_weights"))
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_vline scale_fill_identity
#'   labs theme_minimal theme element_text element_blank element_line
#'   annotate scale_x_continuous margin ggsave
#' @importFrom rlang .data
#' @export
#'

plot_afm_weights <- function(calibration,
                             title = "AFM weights - Phase 1 calibration",
                             highlight_lowest = 0,
                             save_path = NULL) {

  # --- Input validation ---
  if (!is.list(calibration) || !"weights" %in% names(calibration)) {
    stop("'calibration' must be a list from calibrate_afm_mcd() containing ",
         "'weights'. You seem to have passed something else; call ",
         "calibrate_afm_mcd(phase1_data, variables) first and pass its result.")
  }
  weights <- calibration$weights
  if (!is.numeric(weights) || length(weights) < 2 || any(!is.finite(weights))) {
    stop("'calibration$weights' must be a numeric vector of at least 2 finite ",
         "weights. Found ", length(weights), " value(s).")
  }
  if (is.null(names(weights))) {
    stop("'calibration$weights' must be named with the batch identifiers. ",
         "Weights coming from calibrate_afm_mcd() always are.")
  }
  if (!is.character(title) || length(title) != 1) {
    stop("'title' must be a single character string.")
  }
  if (!is.numeric(highlight_lowest) || length(highlight_lowest) != 1 ||
      !is.finite(highlight_lowest) || highlight_lowest < 0 ||
      highlight_lowest != round(highlight_lowest)) {
    stop("'highlight_lowest' must be a single non-negative whole number ",
         "(0 = no batch highlighted).")
  }
  if (highlight_lowest > length(weights)) {
    stop("'highlight_lowest' (", highlight_lowest, ") exceeds the number of ",
         "calibration batches (", length(weights), ").")
  }
  if (!is.null(save_path)) {
    if (!is.character(save_path) || length(save_path) != 1 ||
        nchar(save_path) == 0) {
      stop("'save_path' must be a single non-empty character string, or NULL.")
    }
  }

  # --- Setup ---
  K <- length(weights)
  uniform_weight <- 1 / K

  # Orden ascendente por peso. En ggplot el primer nivel del factor se dibuja
  # abajo, asi que los niveles van de mayor a menor para que el peso mas bajo
  # quede arriba.
  ord <- order(weights, decreasing = FALSE)

  plot_df <- data.frame(
    Batch     = names(weights)[ord],
    Weight    = as.numeric(weights[ord]),
    Rank      = seq_len(K),
    Highlight = seq_len(K) <= highlight_lowest,
    stringsAsFactors = FALSE
  )
  plot_df$Batch <- factor(plot_df$Batch, levels = rev(plot_df$Batch))

  # --- Colours: neutral by default, burgundy only where the caller asked ---
  neutral_color   <- "#8A9AA5"   # slate grey, neutral by design
  highlight_color <- "#A02D31"   # burgundy, same accent as plot_control_chart
  reference_color <- "#2C3E50"   # dark slate for the 1/K line and its label
  plot_df$Fill <- ifelse(plot_df$Highlight, highlight_color, neutral_color)

  # --- Axis range: head-room on the right for the 1/K label ---
  x_max <- max(max(plot_df$Weight), uniform_weight) * 1.18

  ref_label <- paste0("1/K = ", format(round(uniform_weight, 4), nsmall = 4),
                      "  (K = ", K, ")")

  caption_label <- if (highlight_lowest > 0) {
    paste0("The ", highlight_lowest, " lowest-weight batch(es) are ",
           "highlighted at the caller's request. A low weight means high ",
           "internal dispersion, not contamination.")
  } else {
    "A low weight means high internal dispersion, not contamination."
  }

  # --- Build the ggplot object ---
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$Weight, y = .data$Batch, fill = .data$Fill)
  ) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::scale_fill_identity() +

    # Uniform weight reference line
    ggplot2::geom_vline(xintercept = uniform_weight,
                        color = reference_color,
                        linetype = "dashed",
                        linewidth = 0.6) +

    # Inline label for the reference line, at the top of the panel
    ggplot2::annotate("text",
                      x = uniform_weight,
                      y = K,
                      label = ref_label,
                      hjust = -0.06, vjust = 0.4,
                      color = reference_color, size = 3.3) +

    ggplot2::labs(
      title   = title,
      x       = "AFM weight",
      y       = "Batch",
      caption = caption_label
    ) +

    ggplot2::scale_x_continuous(limits = c(0, x_max), expand = c(0, 0)) +

    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(face = "bold", size = 13,
                                                 hjust = 0.5,
                                                 margin = ggplot2::margin(b = 10)),
      plot.caption       = ggplot2::element_text(size = 7.5, hjust = 0,
                                                 color = "grey35",
                                                 margin = ggplot2::margin(t = 8)),
      axis.text.y        = ggplot2::element_text(size = 7.5),
      axis.title         = ggplot2::element_text(face = "plain"),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "grey90",
                                                 linewidth = 0.3),
      plot.margin        = ggplot2::margin(t = 8, r = 12, b = 8, l = 8)
    )

  # --- Optional export to PNG + PDF ---
  # La altura crece con K: 30 barras no caben legibles en 4.5 pulgadas.
  if (!is.null(save_path)) {
    fig_height <- max(4.5, 0.16 * K + 1.5)
    ggplot2::ggsave(paste0(save_path, ".png"), plot = p,
                    width = 7, height = fig_height, units = "in", dpi = 300)
    ggplot2::ggsave(paste0(save_path, ".pdf"), plot = p,
                    width = 7, height = fig_height, units = "in")
  }

  return(p)
}
