context("test-check_design")

test_that("multiplication works", {
  designs <- weight_design(4, 4)
  keep <- check_design(designs)
  expect_true(keep[145])
  expect_false(keep[140])
})
