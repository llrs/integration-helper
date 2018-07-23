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
  ind <- apply(which(upper.tri(result.sgcca$C), arr.ind = TRUE), 1,
               paste0, collapse = "")

  # Correlation between Y
  Y <- simplify2array(result.sgcca$Y, higher = FALSE)
  cY <- cor(Y)
  d <- cY * x
  var <- cY[upper.tri(cY)]
  names(var) <- paste0("vs", ind)
  # Canonical correlation
  cc <- sum(d[upper.tri(d)])
  # Values of the design matrix
  vars <- result.sgcca$C[upper.tri(result.sgcca$C)]
  names(vars) <- paste0("var", ind)
  # weights used
  weight <- sum(vars != 0)
  names(weight) <- "weights"

  # Output
  c(var, unlist(result.sgcca$AVE[c("AVE_inner", "AVE_outer")]), cc1 = cc,
    vars, weight)
}