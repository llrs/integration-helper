test_that("multiplication works", {
  presence <- structure(c(13, 9, 8, 0, 1, 10, 8, 13, 0, 2),
                        .Dim = c(5L, 2L),
                        .Dimnames = list(NULL, c("0", "14")))

  absence <- structure(c(29, 33, 34, 42, 41, 31, 33, 28, 41, 39),
                       .Dim = c(5L, 2L), .Dimnames = list(NULL, c("0", "14")))

  expect_equal(prevalence(presence, absence),
               c(0.625176302023975, 1, 0.214233514422752, 1, 0.615796519410977
  ))


})
