#' Compares two objects of class sgcca
#'
#' Returns the number of genes in common, it assumes the first componenet is
#' always the genes.
#' @param object1,object2 sgcca objects
#' @return Genes In common
#' @export
compare <- function(object1, object2) {

  if (!is(object1, "sgcca") || !is(object2, "sgcca")){
    stop("they should be an sgcca objects from the RGCCA package")
  }

  comp11 <- object1$a[[1]][, 1]
  comp21 <- object2$a[[1]][, 1]

  genes1 <- names(comp11)[diff0(comp11)]
  genes2 <- names(comp21)[diff0(comp21)]

  message("object1: ", length(genes1))
  message("object2: ", length(genes2))
  inter <- intersect(genes1, genes2)
  message("intersection: ", length(inter))

  inter
}



diff0 <- function(x){
  x != 0
}