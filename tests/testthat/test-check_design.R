context("test-check_design")

test_that("weight_design works", {
  designs <- weight_design(4, 4)
  keep <- check_design(designs)
  expect_true(keep[145])
  expect_false(keep[140])
})

test_that("weight_design works with diff0", {
  d <- weight_design(4, 4)[[60]]
  pos <- which(lower.tri(d) & d != 0)
  w <- weight_design(4, 4, pos)
  keep <- sapply(w, function(x){
    all(x[1, 2:4] != 0)
  })
  expect_true(all(keep))
})
