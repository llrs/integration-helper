.onAttach <- function(libname, pkgname) {
  ggplot2::theme_set(ggplot2::theme_bw())
  invisible()
}