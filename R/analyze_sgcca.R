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

  cY <- dimensions_correlation(sgcca)
  cc <- helper_cc(sgcca, cY)

  # Values of the correlation between the dimensions
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

dimensions_correlation <- function(sgcca) {
  # Correlation between Y
  Y <- simplify2array(sgcca$Y, higher = FALSE)
  cor(Y)
}

helper_cc <- function(sgcca, cY) {
  d <- cY * sgcca$C
  switch(sgcca$scheme,
         centroid = sum(abs(d[upper.tri(d)])),
         horst = sum(d[upper.tri(d)]),
         factorial = sum(d[upper.tri(d)]^2))
}

index <- function(x) {
  apply(which(upper.tri(x$C), arr.ind = TRUE), 1,
        paste0, collapse = "")
}

#' Method to simplify AVE
#'
#' This simplifies the AVE_X to make it easier to understand
#' @param x rgcca or sgcca object
#' @return The same object with AVE_X simplified
#' @export
aves <- function(x){
  x$AVE$AVE_X <- simplify2array(x$AVE$AVE_X)
  x
}
