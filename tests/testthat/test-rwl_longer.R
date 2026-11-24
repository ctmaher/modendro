# Need to generate some fake tree-ring data for testing
rwl <- matrix(nrow = 50, ncol = 25)
rwl <- apply(rwl, MARGIN = 2, FUN = \(x) runif(length(x), 0.1, 2)) |>
  as.data.frame()
colnames(rwl) <- paste0("t", 1:25)
rownames(rwl) <- 1951:2000
class(rwl) <- c("rwl","data.frame")
rwl1 <- rwl

rwl.long <- rwl_longer(rwl)

test_that("Throws error if not fed a rwl", {
  expect_error(rwl_longer(c(1:10)))
  expect_error(rwl_longer(1))
})


test_that("Throws error if not fed a rwl-like object", {
  expect_error(rwl_longer(rwl.long))
})

test_that("Throws error if rownames are not consecutive integers", {
  rownames(rwl) <- (1951:2000) + 0.5
  expect_error(rwl_longer(rwl))
  rownames(rwl) <- paste0(1951:2000, "a")
  expect_error(rwl_longer(rwl))
  rownames(rwl) <- c(1951:1970, 1981:2010)
  expect_error(rwl_longer(rwl))
})

test_that("Throws error if any series is all NA", {
  rwl[, "t1"] <- NA
  expect_error(rwl_longer(rwl)) # triggers not numeric error
  rwl[, "t1"] <- NA_real_
  expect_error(rwl_longer(rwl)) # triggers all NAs error
})


