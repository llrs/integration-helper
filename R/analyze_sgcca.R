#' Analyze a sgcca object
#'
#' Calculates the correlation between the first components to reach the
#' canonical correlation as well.
#' @param sgcca SGCCA object from the RGCCA package
#' @return A vector with the correlation between components, AVE (both inner
#' and outer), the canonical correlation, the weight in the design matrix, and
#' the number of interactions that exists.
#' @export
analyze <- function(sgcca) {
  ind <- index(sgcca)

  cY <- canonical_correlation(sgcca)
  cc <- helper_cc(sgcca, cY)

  # Values of the correlation with the dimensions
  var <- cY[upper.tri(cY)]
  names(var) <- paste0("vs", ind)

  # Values of the design matrix
  vars <- sgcca$C[upper.tri(sgcca$C)]
  names(vars) <- paste0("var", ind)

  # weights used
  weight <- sum(vars != 0)
  names(weight) <- "weights"

  # Output
  c(var, unlist(sgcca$AVE[c("AVE_inner", "AVE_outer")]), cc1 = cc,
    vars, weight)
}

canonical_correlation <- function(sgcca) {
  # Correlation between Y
  Y <- simplify2array(sgcca$Y, higher = FALSE)
  cor(Y)
}

helper_cc <- function(sgcca, cY) {
  d <- cY * sgcca$C
  # Canonical correlation
  sum(d[upper.tri(d)])
}

index <- function(x) {
  apply(which(upper.tri(x$C), arr.ind = TRUE), 1,
        paste0, collapse = "")
}
