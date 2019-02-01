# Having a loo list flatten it to a table (only the first dimension)
simplify_weights <- function(x, type) {
  colnames(x) <- paste0("Weights.", type)
  out <- cbind.data.frame(rownames(x), x)
  colnames(out)[1] <- type
  rownames(out) <- NULL
  out
}

name_weights <- function(x, pos) {
  if (pos == 1) {
    type <- "Gene"
  } else if (pos == 2) {
    type <- "Microorganism"
  } else if (pos == 4) {
  } else if (pos == 5) {
  } else {
    type <- pos
    pos <- 3
  }
  simplify_weights(x$a[[pos]], type)
}

remove_zeros <- function(x){
  # Assumes that the first column is the rownames
  k <- rowSums(x[, seq(from = 2, to = ncol(x))])
  x <- x[k != 0, ]
  rownames(x) <- x[, 1]
  as.matrix(x[, -1])
}

#' Correlation weights
#'
#' Select those variables with some weights, filter out those that appear in
#' less than half the list, and make the correlation (using spearman)
#' @param x A list with different sgcca objects (usually from the leave-one-out
#' procedure)
#'
#' @return A matrix with correlations of genes and OTUs
#' @export
weights_correlation <- function(x) {

  sa <- seq_along(x)
  genes_tables <- lapply(x, name_weights, pos = 1)
  genes_df <- Reduce(function(x, y){merge(x, y, by = "Gene")}, genes_tables)
  colnames(genes_df) <- c("Gene", paste0("Weights.", sa))
  genes_df <- remove_zeros(genes_df)


  micro_tables <- lapply(x, name_weights, pos = 2)
  micro_df <- Reduce(function(x, y){merge(x, y, by = "Microorganism")},
                     micro_tables)

  colnames(micro_df) <- c("Microorganism", paste0("Weights.", sa))
  micro_df <- remove_zeros(micro_df)
  # Remove some more
  k <- rowSums(micro_df != 0)
  micro_df <- micro_df[k > max(sa)/2, ]
  k <- rowSums(genes_df != 0)
  genes_df <- genes_df[k > max(sa)/2, ]
  cor(t(micro_df), t(genes_df), method = "spearman")
}



#' Independence between genes and OTUs
#'
#' @param x A list containing sgcc objects
#' @param absolute A logical value to check if the weight should be taken in
#' absolute value or not
#'
#' @return A matrix of p-values of the fisher test to evaluate the independence
#' of the pair.
#' @note We want related pairs, so the higher the p-value, the better.
#' @export
weights_bayes <- function(x, absolute = FALSE) {

  sa <- seq_along(x)
  genes_tables <- lapply(x, name_weights, pos = 1)
  genes_df <- Reduce(function(x, y){merge(x, y, by = "Gene")}, genes_tables)
  colnames(genes_df) <- c("Gene", paste0("Weights.", sa))
  genes_df <- sign(remove_zeros(genes_df))


  micro_tables <- lapply(x, name_weights, pos = 2)
  micro_df <- Reduce(function(x, y){merge(x, y, by = "Microorganism")},
                     micro_tables)

  colnames(micro_df) <- c("Microorganism", paste0("Weights.", sa))
  micro_df <- sign(remove_zeros(micro_df))


  if (absolute) {
    micro_df <- abs(micro_df)
    genes_df <- abs(genes_df)
  }
  # Preallocate
  out <- matrix(0, ncol = nrow(genes_df), nrow = nrow(micro_df),
                dimnames = list(rownames(micro_df), rownames(genes_df)))
  for (i in seq_len(ncol(out))) {
    for (j in seq_len(nrow(out))) {
      out[j, i] <- helper_fisher(micro_df[j, ], genes_df[i, ],
                                 absolute = absolute)
    }
  }
  out
}


helper_fisher <- function(x, y, absolute) {
  if (!absolute) {
    m <- matrix(0, ncol = 3, nrow = 3, dimnames = list(c("1", "0", "-1"),
                                                       c("1", "0", "-1")))
  } else {
    m <- matrix(0, ncol = 2, nrow = 2, dimnames = list(c("1", "0"),
                                                       c("1", "0")))

  }
  tm <- table(x, y)
  dmTm <- dimnames(tm)
  m[dmTm$x, dmTm$y] <- tm
  fisher.test(m)$p.value
}
