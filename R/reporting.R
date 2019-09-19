#' Select the weights and adds information for human RNAseq
#'
#' @param x a RGCCA outputs
#' @return a data.frame
#' @importFrom AnnotationDbi select
#' @export
weights <- function(x) {
  loadings <- x$a$RNAseq
  colnames(loadings) <- paste0("Weights_", seq_len(ncol(loadings)))
  ensemblID <- rownames(loadings)
  ensemblID <- trimVer(ensemblID)
  rownames(loadings) <- trimVer(rownames(loadings))
  symbolID <- select(org.Hs.eg.db,
    keys = ensemblID, keytype = "ENSEMBL", columns = c("SYMBOL", "GENENAME")
  )
  a <- match(symbolID$ENSEMBL, rownames(loadings))
  rownames(loadings) <- seq_len(nrow(loadings))
  loadings <- loadings[a, ]
  rownames(loadings) <- seq_len(nrow(loadings))
  out <- cbind(loadings, symbolID)

  keep <- apply(out[, seq_len(ncol(loadings))], 1, function(x) {
    any(x != 0)
  })
  out <- out[keep, ]
  rownames(out) <- seq_len(nrow(out))
  out
}

#' Select the weights and adds information for taxa RNAseq
#'
#' @param x a RGCCA outputs
#' @param taxa The taxonomic levels
#' @return a data.frame
#' @export
weights_otus <- function(x, taxa) {
  loadings <- x$a$`16S`
  colnames(loadings) <- paste0("Weights_", seq_len(ncol(loadings)))
  taxa <- taxa[rownames(loadings), ]
  out <- cbind(loadings, taxa)
  keep <- apply(out[, seq_len(ncol(loadings))], 1, function(x) {
    any(x != 0)
  })
  out <- out[keep, ]
  rownames(out) <- seq_len(nrow(out))
  out
}
