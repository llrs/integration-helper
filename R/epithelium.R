#' Function to read the epithelium signature
#'
#' Return the Epithelium in Entrez ids
#' @param x A character vector with SYMBOL names
#' @return A character vector with the Entrez ids
#' @export
epitheliumE <- function(x){
  if (!requireNamespace("AnnotationDbi")) {
    stop("Install AnnotationDbi from Bioconductor", call. = FALSE)
  }
  if (!requireNamespace("org.Hs.eg.db")) {
    stop("Install org.Hs.eg.db from Bioconductor", call. = FALSE)
  }
  epitheliumE <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = as.character(x),
    keytype = "SYMBOL", column = "ENTREZID"
  )

  epitheliumE <- unlist(epitheliumE, use.names = TRUE)
  epitheliumE[!is.na(epitheliumE)]
}
