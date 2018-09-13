#' Function to read the epithelium signature
#'
#' Return the Epithelium in Entrez ids
#' @param x A character vector with SYMBOL names
#' @return A character vector with the Entrez ids
#' @export
epitheliumE <- function(x){
  epitheliumE <- mapIds(
    org.Hs.eg.db,
    keys = as.character(x),
    keytype = "SYMBOL", column = "ENTREZID"
  )

  epitheliumE <- unlist(epitheliumE, use.names = TRUE)
  epitheliumE[!is.na(epitheliumE)]
}