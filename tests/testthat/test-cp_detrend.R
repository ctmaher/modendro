test_that("Throws error if not fed a rwl, data.frame, or matrix", {
  expect_error(cp_detrend(c(1:10)))
  expect_error(cp_detrend(1))
})

test_that("Throws error if fed data with no rownames (years) or colnames (series IDs)", {
  mat.test <- matrix(nrow = 10, ncol = 10)
  expect_error(cp_detrend(mat.test))
  rownames(mat.test) <- 1:10
  expect_error(cp_detrend(mat.test))
  rownames(mat.test) <- NULL
  colnames(mat.test) <- 1:10
  expect_error(cp_detrend(mat.test))
})

test_that("Throws error if any series are all NAs", {
  mat.test <- matrix(nrow = 10, ncol = 10)
  rownames(mat.test) <- 1:10
  colnames(mat.test) <- 1:10
  mat.test[, 1:5] <- runif(10)
  expect_error(cp_detrend(mat.test))
})

test_that("Names (series IDs) are equal in the input rwl and output rwl", {
  mat.test <- matrix(nrow = 50, ncol = 10)
  rownames(mat.test) <- 1:50
  colnames(mat.test) <- 1:10
  mat.test <- apply(mat.test, MARGIN = 2, FUN = \(x) runif(length(x), 0.1, 2))
  out.test <- cp_detrend(mat.test)[[1]]
  expect_equal(colnames(mat.test), colnames(out.test))
})

####
test_that("Essential function works - and returns a list", {
  # Make some fake tree ring data
  mat.test <- matrix(nrow = 50, ncol = 10)
  rownames(mat.test) <- 1:50
  colnames(mat.test) <- 1:10
  mat.test <- apply(mat.test, MARGIN = 2, FUN = \(x) runif(length(x), 0.1, 2))
  expect_vector(cp_detrend(mat.test), ptype = list())
})
