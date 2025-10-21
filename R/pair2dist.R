#' Helper function to transform long matrix to distance matrix
#' @param dat Long matrix with three columns (site1, site2 and values)
#' @param colname Name of the column that contains distance value
#' @return Object of class \code{dist}
#' @export

pair2dist <- function(dat, colname) {
  dat <- dat[order(dat$site1, dat$site2), ]
  sites <- sort(unique(c(dat$site1, dat$site2)))
  vals <- dat[[colname]]
  attr(vals, "Size") <- length(sites)
  attr(vals, "Labels") <- sites
  attr(vals, "Diag") <- FALSE
  attr(vals, "Upper") <- FALSE
  class(vals) <- "dist"
  vals
}
