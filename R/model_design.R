#'
#' Check a design for different weights
#' @param weights A numeric value with the number of weights to be used
#' @param size A numveric value with the numer of datasets on the design.
#' @return A list of matrices with the designs with different weights
#' @export
#' @author Flodel \url{https://codereview.stackexchange.com/a/203517/36067}
#' @examples
#' out <- weight_design(4, 4)
#' head(out)
weight_design <- function(weights = 4, size){

  p <- size * (size - 1) / 2    # 6
  w <- seq(from = 0, to = 1, length.out = weights)

  # all possible combinations by doing:
  W <- as.matrix(expand.grid(rep(list(w), p)))

  X <- matrix(1:(size*size), size, size) # pattern matrix of indices
  A <- matrix(0, nrow(W), size * size)

  # Replace the positions by the weights
  A[,    X[lower.tri(X)]] <- W
  A[, t(X)[lower.tri(X)]] <- W

  # A 3D array with the weights
  dim(A) <- c(nrow(W), size, size)

  # Convert to a list
  lapply(seq(nrow(A)), function(i) A[i, , ])
}

#' Validate designs
#'
#' @param designs A list of design matrices as obtained from
#' \code{\link{weight_design}}
#' @return A logical vector that validates if the designs are correct or not.
#' @export
#' @seealso \code{\link{symm}}, \code{\link{subSymm}}
#' @examples
#' designs <- weight_design(4, 4)
#' keep <- check_design(designs)
#' summary(keep)
check_design <- function(designs) {
  # Validation of the input
  stopifnot(is(designs, "list"))
  nCols <- vapply(designs, ncol, numeric(1L))
  nRows <- vapply(designs, nrow, numeric(1L))
  stopifnot(length(unique(nCols)) == 1)
  stopifnot(length(unique(nRows)) == 1)
  stopifnot(unique(nRows) == unique(nCols))
  # Actual function
  vapply(designs, valid,  logical(1L))
}

# Check that each data block is connected to others
valid <- function(x){
  all(rowSums(x != 0) != 0)
}


# Check that each data block is to all
# See this answer:
# https://math.stackexchange.com/a/551947
# There should work
#' Check that the network is fully connected
#'
#' Given the design matrix, checks that all the blocks are connected between them
#' @param x Design matrix, a symmetric matrix with
#' @return A logical value if it is fully connected or not.
#' @references \url{https://math.stackexchange.com/a/551947}
#' @export
correct <- function(x){
  A <- x != 0 # Adjacency
  # Repeat the adjaceny as much as it is needed.
  l <- lapply(seq_len(ncol(A) - 1), function(y){A})
  # Calculate the power (there are more efficient ways but for small matrices it should work)
  red <- Reduce(`%*%`, l, init = A, accumulate = TRUE)
  # Add them up (S)
  final <- Reduce(`+`, red)
  all(final != 0)
}

#' Prepare data
#'
#' Prepares the factors into their vectors.
#' @param data A data.frame with the information about the samples
#' @param columns The name of the columns to be used to build the matrix
#' @param intercept A logical value if you want one column with all 1 or not.
#' @return A matrix with each factor is decomposed in as much columns as
#' factors has minus 1 and with the numeric values as they were.
#' @export
model_RGCCA <- function(data, columns, intercept = FALSE){

  m <- data[, columns, drop = FALSE]
  num <- vapply(m, is.numeric, logical(1L))
  if (any(!num)) {
    if (sum(!num) > 1) {
      o <- sapply(m[, !num, drop = FALSE], function(x){
        levels <- unique(x)
        levels <- levels[!is.na(levels)]
        o <- vapply(levels, function(level) {
          as.numeric(x %in% level)
        }, numeric(nrow(data)))
        o[, -1, drop = FALSE]
      })
      o <- do.call(cbind, o)
    } else {
      levels <- unique(m[, !num])
      levels <- levels[!is.na(levels)]
      o <- vapply(levels, function(level) {
        as.numeric(m[, !num] %in% level)
      }, numeric(nrow(data)))
      o <- o[, -1, drop = FALSE]

    }
  }

  if (any(!num) & any(num)) {
    out <- cbind(o, m[, num, drop = FALSE])
  } else if (any(!num)) {
    out <- o
  } else {
    out <- m[, num, drop = FALSE]
  }

  colnames(out)[colnames(out) == ""] <- seq_len(sum(colnames(out) == ""))

  if (intercept) {
    cbind(1, out)
  } else {
    out
  }
}

#' Create symmetric matrix
#'
#' @param m Square matrix.
#' @param data Numeric values of the upper triangular side of the matrix
#' @note After the upper case there can be other values that are assumed to be
#' in the diagonal.
#' @return A square symmetric matrix.
#' @seealso \code{\link{subSymm}}, \code{\link{correct}}, \code{\link{check_design}}
#' @export
symm <- function(m, data) {

  if (!is(data, "numeric")) {
    data <- unlist(data)
  }
  weights <- ncol(m)^2 - ncol(m)
  upper <- data[seq_len(weights)]
  if (length(data) > weights) {
    Diag <- rep(0, ncol(m))
    Diag[seq_len(length(data) - weights)] <- data[seq_along(data) > weights]
    diag(m) <- Diag
  }
  m[upper.tri(m)] <- upper
  as.matrix(Matrix::forceSymmetric(m, "U"))
}
