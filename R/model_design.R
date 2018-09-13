#'
#' Check a design for different weights
#' @param weights A numeric value with the number of weights to be used
#' @param size A numveric value with the numer of datasets on the design.
#' @return A list of matrices with the designs with different weights
#' @export
#' @examples
#' weight_design(4, 4)
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
  if (any(!num)){
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
