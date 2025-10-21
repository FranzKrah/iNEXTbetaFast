#' Helper functions
#' @param fun Function to be loaded from a package function that is not exported from that package (not in NAMESPACE)
#' @param pkg R package from which to load the function from
#' @details Workaround proposed here \link{https://stackoverflow.com/questions/63023526/unexported-object-imported-by-a-call-tsfeaturesscalets}
#' @importFrom utils getFromNamespace
#' @keywords internal

safe_get <- function(fun, pkg) {
  if (fun %in% getNamespaceExports(pkg)) {
    getExportedValue(pkg, fun)
  } else if (exists(fun, envir = asNamespace(pkg), inherits = FALSE)) {
    getFromNamespace(fun, pkg)
  } else {
    stop(sprintf("Function '%s' not found in package '%s'", fun, pkg))
  }
}

