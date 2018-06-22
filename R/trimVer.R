
#' Trim version number of genes
#'
#' Remove the trailing version of gene names in the ENSEMBEL versioning
#' scheme.
#' @return A vector of genes names
#' @param genes A character vector.
#' @examples
#' trimVer("ENSG00000049192.14")
#' trimVer(c("ENSG00000049192.14", "ENSG00000065320.8"))
#' @export
trimVer <- function(genes){
  gsub("\\..+", "", genes)
}
