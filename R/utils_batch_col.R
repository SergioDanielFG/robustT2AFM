# Shared check for the column that identifies the batch.
#
# Naming that column "Batch" is only the package default, not a requirement of
# the method, so the error has to name what was looked for, list what is
# actually there, and show how to point at the right column. Getting this
# wrong is the most likely first failure for a new user.
#
# Not exported: it is an internal helper for the five functions that read a
# user-supplied data frame.
check_batch_col <- function(data, batch_col, arg_name, call_hint) {
  if (!is.character(batch_col) || length(batch_col) != 1 ||
      is.na(batch_col) || nchar(batch_col) == 0) {
    stop("'batch_col' must be a single non-empty column name.", call. = FALSE)
  }
  if (!batch_col %in% colnames(data)) {
    stop("Column '", batch_col, "' not found in '", arg_name,
         "'. The available columns are: ",
         paste(colnames(data), collapse = ", "),
         ". If your batch identifier has a different name, pass it: ",
         call_hint, ".", call. = FALSE)
  }
  invisible(TRUE)
}
