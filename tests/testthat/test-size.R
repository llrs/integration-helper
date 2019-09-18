test_that("size works", {
  expect_equal(size(vector(length = 3L)), 3L)
  expect_equal(size(matrix()), 1L)
  expect_equal(size(list(a = matrix())), 1L)
})
