#' Significative correlation
#'
#' Calculates the correlation at which they become significant
#' @param n The number of data points where the correlation has been done
#' @return the correlation value at which it becomes below 0.05
#' @export
cor_sign <- function(n) {
  stopifnot(n > 2)
  stopifnot(is.numeric(n))

  df <- n - 2
  r <- 0
  pval <- 1
  # Compare it with https://stackoverflow.com/a/45412227/2886003
  while (pval >= 0.05) {
    r <- r + 0.0005
    t <- sqrt(df) * r / sqrt(1 - r^2)
    pval <- 2 * min(pt(t, df), pt(t, df, lower.tail = FALSE)) ## th
  }
  r
}

#' Calculates the p.value
#'
#' @param r Correlation coefficient
#' @param n Number of samples
#' @return the p-value
#' @export
pvalue <- function(r, n) {
  df <- n - 2
  if (is.matrix(r)) {
    apply(r, 1:2, function(x) {
      t <- sqrt(df) * x / sqrt(1 - x^2)
      2 * min(pt(t, df), pt(t, df, lower.tail = FALSE))
    })
  } else {
    t <- sqrt(df) * r / sqrt(1 - r^2)
    2 * min(pt(t, df), pt(t, df, lower.tail = FALSE)) ## th
  }
}
