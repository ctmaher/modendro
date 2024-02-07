####
test_that("Throws error if rwl is not a rwl, data.frame, or matrix", {
  expect_error(yrs_to_pith(c(1:10)))
  expect_error(yrs_to_pith(1))
})
####
test_that("Throws error if d2pith is not data.frame or matrix", {
  expect_error(yrs_to_pith(rwl = matrix(nrow = 10, ncol = 10), d2pith = c(1,10,1)))
  expect_error(yrs_to_pith(rwl = matrix(nrow = 10, ncol = 10), d2pith = 1))
})
####
test_that("Throws error if n.rings is not a length = 1 numeric", {
  expect_error(yrs_to_pith(rwl = matrix(nrow = 10, ncol = 10),
                           d2pith = data.frame(series = colnames(mat.test),
                                               d2pith = runif(ncol(mat.test), 0.2, 5)),
                           n.rings = "5"))
  expect_error(yrs_to_pith(rwl = matrix(nrow = 10, ncol = 10),
                           d2pith = data.frame(series = colnames(mat.test),
                                               d2pith = runif(ncol(mat.test), 0.2, 5)),
                           n.rings = c(5,10)))
})
####
test_that("Throws error if fed an rwl with no rownames (years) or colnames (series IDs)", {
  mat.test <- matrix(nrow = 10, ncol = 10)
  expect_error(yrs_to_pith(mat.test))

  rownames(mat.test) <- 1:10
  expect_error(yrs_to_pith(mat.test))

  rownames(mat.test) <- NULL
  colnames(mat.test) <- paste0("series", 1:10)
  expect_error(yrs_to_pith(mat.test))
})
####
test_that("Throws error if any rwl series are all NAs", {
  mat.test <- matrix(nrow = 10, ncol = 10)
  rownames(mat.test) <- 1:10
  colnames(mat.test) <- paste0("series", 1:10)
  mat.test[, 1:5] <- runif(10)
  expect_error(yrs_to_pith(mat.test))
})
####
test_that("Names are equal in the rwl and messages outputs", {
  mat.test <- matrix(nrow = 50, ncol = 10)
  mat.test <- apply(mat.test, MARGIN = 2, FUN = \(x) runif(length(x), 0.1, 2)) |> as.data.frame()
  rownames(mat.test) <- 1:50
  colnames(mat.test) <- paste0("series", 1:10)

  d2pith.test <- data.frame(series = colnames(mat.test), d2pith = runif(ncol(mat.test), 0.2, 5))
  out.test <- yrs_to_pith(mat.test, d2pith = d2pith.test)
  expect_equal(colnames(out.test[[1]]), names(out.test[[2]]))
})
