.onAttach <- function(libname, pkgname) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::theme_set(ggplot2::theme_bw())
  }
  invisible()
}
