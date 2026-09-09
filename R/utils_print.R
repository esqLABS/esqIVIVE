# Internal helper for verbose = TRUE output across exported functions.
# Not exported: prints the inputs and result of a function call in a
# consistent format, gated behind each function's verbose argument.
.print_ivive_result <- function(fn_name, inputs, result) {
  cat("---", fn_name, "---\n")
  cat("Inputs:\n")
  for (nm in names(inputs)) {
    cat(sprintf("  %s = %s\n", nm, paste(inputs[[nm]], collapse = ", ")))
  }
  cat("Result:\n")
  if (is.list(result) || (is.vector(result) && !is.null(names(result)))) {
    for (nm in names(result)) {
      cat(sprintf("  %s = %s\n", nm, paste(result[[nm]], collapse = ", ")))
    }
  } else {
    cat(sprintf("  %s\n", paste(result, collapse = ", ")))
  }
  invisible(NULL)
}
