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
  x[, -1]
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
