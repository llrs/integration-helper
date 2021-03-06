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


#' Improve the information on sgcca classes
#'
#' Add names to data, simplify AVE output
#' @param sgcca An object of class `sgcca`.
#' @param namesA The names of the original data
#' @return An object of class sgcca
#' @export
improve.sgcca <- function(sgcca, namesA) {
  if (is.null(namesA)) {
    stop("namesA shouldn't be NULL\n",
         "Consider adding names to A.")
  }
  names(sgcca$Y) <- namesA
  names(sgcca$a) <- namesA
  names(sgcca$astar) <- namesA
  names(sgcca$AVE$AVE_X) <- namesA
  colnames(sgcca$C) <- namesA
  rownames(sgcca$C) <- namesA
  aves(sgcca)
}

#' Extract in a tidy way the information about the variables
#'
#' @param sgcca A `sgcca` object with named rows.
#' @return A data.frame of variables, component and origin.
#' @export
#' @seealso [plot_variables()]
variables <- function(sgcca){
  if (!is(sgcca, "sgcca")) {
    stop("Not of sgcca class")
  }

  a_rownames <- lapply(sgcca$a, rownames)
  a_len <- vapply(sgcca$a, nrow, numeric(1L))
  if (any(lengths(a_rownames) != a_len)) {
    stop("Some rows are not named.")
  }
  variables <- data.frame(
    comp1 = unlist(lapply(sgcca$a, function(x) {x[, 1]}),
                   use.names = FALSE, recursive = FALSE),
    comp2 = unlist(lapply(sgcca$a, function(x) {x[, 2]}),
                   use.names = FALSE, recursive = FALSE),
    Origin = rep(names(sgcca$a), lengths(a_rownames))
  )
  variables$var <- unlist(a_rownames, use.names = FALSE, recursive = FALSE)
  variables
}
