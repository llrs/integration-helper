test_that("clean_unvariable works", {
  x <- list(x = matrix(runif(10), ncol = 2, nrow = 5),
            y = matrix(c(0, 0, 0, 0, 0, 1, 2, 3, 4, 5), ncol = 2, nrow = 5))
  xb <- clean_unvariable(x)
  expect_equal(names(xb), names(x))
  expect_equal(dim(xb$y), c(5, 1))
})
