context("test-allcomb")

test_that("allComb works", {
  disease <- factor(sample(c("disease", "control"), 15, TRUE))
  ileum <- factor(sample(c("ileum", "colon"), 15, TRUE))
  df <- cbind.data.frame(disease, ileum)
  out <- allComb(df, c("disease", "ileum"))

  expect_true(any(out))
  expect_true(all(apply(out, 1, sum) == 1))
})
