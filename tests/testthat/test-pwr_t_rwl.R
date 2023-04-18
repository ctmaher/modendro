test_that("Throws error if not fed a rwl, data.frame, or matrix", {
  expect_error(pwr_t_rwl(c(1:10)))
  expect_error(pwr_t_rwl(1))
})

test_that("Throws error if fed data with no rownames (years) or colnames (series IDs)", {
  mat.test <- matrix(nrow = 10, ncol = 10)
  expect_error(pwr_t_rwl(mat.test))
  rownames(mat.test) <- 1:10
  expect_error(pwr_t_rwl(mat.test))
  rownames(mat.test) <- NULL
  colnames(mat.test) <- 1:10
  expect_error(pwr_t_rwl(mat.test))
})

test_that("Throws error if any series are all NAs", {
  mat.test <- matrix(nrow = 10, ncol = 10)
  rownames(mat.test) <- 1:10
  colnames(mat.test) <- 1:10
  mat.test[, 1:5] <- runif(10)
  expect_error(pwr_t_rwl(mat.test))
})

test_that("Names are equal in the rwl and messages outputs", {
  mat.test <- matrix(nrow = 50, ncol = 10)
  rownames(mat.test) <- 1:50
  colnames(mat.test) <- 1:10
  mat.test <- apply(mat.test, MARGIN = 2, FUN = \(x) runif(length(x), 0.1, 2))
  out.test <- pwr_t_rwl(mat.test)
  expect_equal(colnames(out.test[[1]]), names(out.test[[2]]))
})
