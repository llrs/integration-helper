#' Intracorrelation
#'
#' Calculates the intercorrelation between ids
#' (the correlation between different ids)
#' @param A Matrix
#' @param ids Vector of ids of the sample
#' @return the pearson coefficient
intercorrelation <- function(A, ids) {
  stopifnot(is.matrix(A))
  stopifnot(length(ids) != ncol(A))

  x <- cor(A)
  sameIDs <- outer(ids, ids, "==")
  unlist(x[upper.tri(x) & !sameIDs], use.names = FALSE)
}

#' Intracorrelation
#'
#' Calculates the intracorrelation between ids
#' @param A Matrix
#' @param ids Vector of ids of the sample
#' @return the pearson coefficient
intracorrelation <- function(A, ids) {
  stopifnot(is.matrix(A))
  stopifnot(length(ids) != ncol(A))

  x <- cor(A)
  sameIDs <- outer(ids, ids, "==")
  unlist(x[upper.tri(x) & sameIDs], use.names = FALSE)
}
