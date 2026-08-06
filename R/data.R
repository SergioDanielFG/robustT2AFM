#' Phase 1 Calibration Batches (Simulated, Contaminated)
#'
#' \strong{These are simulated data, and two of their columns could not exist
#' in a real process.} \code{Status} and \code{ContaminationType} record which
#' batches are contaminated and how. That is knowable here only because the
#' data was generated on purpose; on a plant floor it is precisely the thing
#' the control chart is being asked to find out. Do not read this data set as
#' a template implying your own historical data should carry such a column: it
#' will not, and nothing in the method needs it. The two columns exist so that
#' the package can demonstrate validation against known truth, and the only
#' function that consumes them is the \code{faulty} argument of
#' \code{\link{plot_method_comparison}}.
#'
#' Thirty Phase 1 batches of 20 observations on 4 correlated process
#' variables, of which 6 batches carry 4 outlying observations each. This is
#' the Phase 1 of the base configuration in Frutos-Galarza et al. (2026), and
#' it is the calibration set used by every example in the package.
#'
#' @format A data frame with 600 rows and 8 columns:
#' \describe{
#'   \item{Batch}{Factor, batch identifier (\code{F1_B01} to \code{F1_B30}).}
#'   \item{Phase}{Factor, always "Phase 1" here.}
#'   \item{Status}{Factor, "Under Control" or "Out of Control". Ground truth,
#'         see the warning above.}
#'   \item{ContaminationType}{Factor, "Clean" or "Outliers". Ground truth,
#'         see the warning above.}
#'   \item{Var1, Var2, Var3, Var4}{Numeric process variables, equicorrelated
#'         with rho = 0.6 in control.}
#' }
#'
#' @details
#' The examples in this package write \code{data(afm_phase1)} before using the
#' object. With \code{LazyData} enabled that call is not strictly necessary,
#' since the object is available as soon as the package is attached. It is
#' kept because it tells the reader where the object comes from; without it
#' the name would appear out of nowhere.
#'
#' @source
#' Generated with:
#' \preformatted{
#' simulate_batch_process(
#'   K1 = 30, K2 = 20, I = 20, J = 4, rho = 0.6,
#'   outlier_batches_F1 = 6, outlier_rate = 0.20, outlier_shift = 4,
#'   prop_ooc_F2 = 0.5, shift_ooc = 1.0,
#'   seed = 20260425
#' )
#' }
#' then split by phase, with unused factor levels dropped. The full script is
#' \code{data-raw/make_datasets.R}, and \code{tests/testthat/test-data.R}
#' regenerates the data from that call and checks it against what is shipped.
#'
#' @seealso \code{\link{afm_phase2}} for the batches to monitor,
#'   \code{\link{simulate_batch_process}} to generate other scenarios.
#'
#' @examples
#' data(afm_phase1)
#' str(afm_phase1)
#' table(afm_phase1$ContaminationType)   # 480 clean rows, 120 with outliers
"afm_phase1"


#' Phase 2 Monitoring Batches (Simulated, Half Off-Target)
#'
#' \strong{These are simulated data, and two of their columns could not exist
#' in a real process.} \code{Status} and \code{ContaminationType} record which
#' batches are off-target. That is knowable here only because the data was
#' generated on purpose; on a plant floor it is precisely the thing the
#' control chart is being asked to find out. Do not read this data set as a
#' template implying your own monitoring data should carry such a column: it
#' will not, and nothing in the method needs it. The two columns exist so that
#' the package can demonstrate validation against known truth, and the only
#' function that consumes them is the \code{faulty} argument of
#' \code{\link{plot_method_comparison}}.
#'
#' Twenty Phase 2 batches of 20 observations on the same 4 variables, of which
#' 10 are shifted by 1 standard deviation. This is the Phase 2 of the base
#' configuration in Frutos-Galarza et al. (2026), and it is the monitoring set
#' used by every example in the package.
#'
#' @format A data frame with 400 rows and 8 columns:
#' \describe{
#'   \item{Batch}{Factor, batch identifier (\code{F2_B01} to \code{F2_B20}).}
#'   \item{Phase}{Factor, always "Phase 2" here.}
#'   \item{Status}{Factor, "Under Control" or "Out of Control". Ground truth,
#'         see the warning above.}
#'   \item{ContaminationType}{Factor, "Clean" or "OOC". Ground truth, see the
#'         warning above.}
#'   \item{Var1, Var2, Var3, Var4}{Numeric process variables.}
#' }
#'
#' @details
#' The examples in this package write \code{data(afm_phase2)} before using the
#' object. With \code{LazyData} enabled that call is not strictly necessary,
#' since the object is available as soon as the package is attached. It is
#' kept because it tells the reader where the object comes from; without it
#' the name would appear out of nowhere.
#'
#' To extract the ground truth for \code{plot_method_comparison}:
#' \preformatted{
#' faulty <- as.character(unique(
#'   afm_phase2$Batch[afm_phase2$Status == "Out of Control"]
#' ))
#' }
#'
#' @source
#' Same call as \code{\link{afm_phase1}}; see that help page and
#' \code{data-raw/make_datasets.R}.
#'
#' @seealso \code{\link{afm_phase1}} for the calibration batches.
#'
#' @examples
#' data(afm_phase2)
#' str(afm_phase2)
#' table(afm_phase2$Status)   # 10 batches in control, 10 off-target
"afm_phase2"
