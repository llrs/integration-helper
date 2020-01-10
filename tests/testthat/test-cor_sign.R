test_that("cor_sign works", {
  expect_equal(cor_sign(3), 0.996999999999946)
})

test_that("pvalue works", {
  expect_equal(pvalue(0.5, 5), 0.391002218955771)

  o <- pvalue(matrix(runif(25), ncol = 5, nrow = 5), 5)
  expect_equal(dim(o), c(5, 5))
  expect_true(is(o, "matrix"))
})
