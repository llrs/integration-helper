#' Indexes without one sample
#'
#' Performs the leave-out-one indexing
#' @param size The number of samples
#' @return A list of indices where one sample is excluded
#' @export
#' @examples
#' looIndex(15)
looIndex <- function(size){

  l <- seq_len(size)
  vl <- vector("list", length = size)
  for (i in l) {
    vl[[i]] <- l[-i]
  }
  vl
}

#' Number of samples
#'
#' @param x A list or a matrix where rows are the samples
#' @return The number of samples
#' @export
size <- function(x){
  if (is(x, "list")) {
    stopifnot(length(unique(vapply(x, nrow, numeric(1L)))) == 1L)
    nrow(x[[1]])
  } else {
    nrow(x)
  }
}

# Given an index keep those of the list of matrices or data.frames
#' Subset a list
#'
#' Given a list and an index subset each element of the list.
#' @param A A list of an array with samples in rows.
#' @param index The number of samples to keep
#' @return A list with less samples
#' @export
subsetData <- function(A, index) {
  lapply(A, function(x){
    y <- x[index, , drop = FALSE] # subset
    y[, apply(y, 2, sd) != 0] # Remove variables that are constant.
    })
}



#  K-fold ####
#' Indexes for K-folds
#'
#' Calcules the indices of the k-fold bootstrapping for other functions
#' @param k Number of k-folds (two divides the dataset in 2)
#' @param n Number of k-folds to do
#' @param size The number of samples required
#' @param ... Other arguments passed to sample the
#' @return An list of length \code{n} with two elements the training and the
#' testing elements with the index of the samples
#' @export
#' @seealso \code{\link{looIndex}}
kFolding <- function(k, size, n, ...) {

  stopifnot(k < size)
  if (n < k) {
    warning("Less k-folds samples than possible")
  }
  if (n < size) {
    warning("More k-folds than samples")
  }
  keep <- ceiling(size/k)
  if (n > choose(size, keep)) {
    warning("More k-folds than possibilities")
  }

  x <- seq_len(size)
  out <- vector("list", length = n)

  for (i in seq_len(n)) {
    samples <- sample(x, keep, ...)
    # to avoid random subsetting  it could use something of
    # seq(5, length.out = 15) To select the testing group
    out[[i]] <- list("testing" = NULL, "training" = NULL)
    out[[i]]$testing <- samples
    out[[i]]$training <- x[!x %in% samples]
  }
  out
}
