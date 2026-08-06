# Shared input checks. Kept in one place so the same mistake produces the same
# message wherever the user hits it.
#
# The rule they follow: an error says what to do, not only what failed. For a
# "not found" mistake that means naming what was looked for AND listing what is
# actually there, because the user cannot fix a typo they cannot see.


# Shared check for the process variables.
#
# The batch-column check below listed the available columns from the start;
# this one used to say only "Variables not found in data: Nope", three lines
# away in the same function. Same mistake, half the help.
check_variables <- function(data, variables, arg_name, call_hint) {
  if (!is.character(variables) || length(variables) < 1) {
    stop("'variables' must be a non-empty character vector with the names of ",
         "the process variables to monitor.", call. = FALSE)
  }
  missing_vars <- setdiff(variables, colnames(data))
  if (length(missing_vars) > 0) {
    numeric_cols <- colnames(data)[vapply(data, is.numeric, logical(1))]
    stop("Variable(s) not found in '", arg_name, "': ",
         paste(missing_vars, collapse = ", "),
         ". The numeric columns available are: ",
         if (length(numeric_cols) > 0) {
           paste(numeric_cols, collapse = ", ")
         } else {
           "none"
         },
         ". Pass the ones you want to monitor: ", call_hint, ".",
         call. = FALSE)
  }
  invisible(TRUE)
}


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
