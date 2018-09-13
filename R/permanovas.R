#' Permanova
#'
#' @param otus Matrix of the variables to analyze
#' @param meta Data frame with the variables to compare
#' @return An anova.cca object
permanova_otus <- function(otus, meta) {
  adonis(as.matrix(otus) ~ ., data = meta, method = "jaccard")
}

#' Permanova
#'
#' @param expr Matrix of the expression
#' @param meta Data frame with the variables to compare
#' @return An anova.cca object
permanova_expr <- function(expr, meta) {
  diss <- 1 - cor(expr)
  adonis(diss ~ ., data = meta, method = "euclidian")
}

#
# b <- adonis(as.data.frame(t(otus_table_i)) ~ AGE_SAMPLE + diagTime + AgeDiag, data = (meta_r[, nam]), method = "jaccard")
#
# > b <- adonis(t(as.matrix(otus_table_i)) ~ AGE_SAMPLE + diagTime + AgeDiag, data = as.data.frame(t(meta_r[, nam])), method = "jaccard")
# Error in eval(predvars, data, env) : object 'AGE_SAMPLE' not found
