####
test_that("Throws error if not fed a chrono, data.frame, or matrix", {
  expect_error(n_mon_corr(chrono = c(1:10), clim = c(1:10)))
  expect_error(n_mon_corr(chrono = 1, clim = 1))
})
####
test_that("Throws error if not fed a valid clim var name", {
  clim <- matrix(nrow = 10, ncol = 3)
  colnames(clim) <- c("year", "month", "clim.var")
  expect_error(n_mon_corr(chrono = matrix(nrow = 10, ncol = 2), clim = clim,
                          var = "wack-a-doodle"))
})
####
test_that("Throws error if not fed a valid chrono.col name", {
  clim <- matrix(nrow = 10, ncol = 3)
  colnames(clim) <- c("year", "month", "clim.var")
  chrono <- matrix(nrow = 10, ncol = 2)
  colnames(chrono) <- c("samp.depth", "std")
  expect_error(n_mon_corr(chrono = chrono, clim = clim,
                          var = "clim.var", chrono.col = "wack-a-doodle"))
})
####
test_that("Throws error if not fed a valid aggregation function", {
  expect_error(n_mon_corr(chrono = matrix(nrow = 10, ncol = 2), clim = matrix(nrow = 10, ncol = 3),
                          agg.fun = "wack-a-doodle"))
})
####
test_that("Throws error if not fed a valid correlation method", {
  clim <- matrix(nrow = 10, ncol = 3)
  colnames(clim) <- c("year", "month", "clim.var")
  chrono <- matrix(nrow = 10, ncol = 2)
  colnames(chrono) <- c("samp.depth", "std")
  expect_error(n_mon_corr(chrono = chrono, clim = clim,
                          var = "clim.var", chrono.col = "std",
                          corr.method = "catty-wompus"))
})
####
test_that("Essential function works - and returns a vector", {
  # Make some fake tree ring data
  mat.test <- matrix(nrow = 50, ncol = 10)
  rownames(mat.test) <- 1:50
  colnames(mat.test) <- 1:10
  mat.test <- apply(mat.test, MARGIN = 2, FUN = \(x) runif(length(x), 0.1, 2))
  chrono <- dplR::chron(mat.test)
  # Make some fake climate data
  mo <- 1:12
  x <- seq(-4, 4, length.out = 12)
  gauss.curv <- \(x) {(10/sqrt(2*pi*1.6))*exp(-((x^2)/(2*1.6^2)))}
  clim.mo <- data.frame(month = mo, clim.var = gauss.curv(x))
  clim.list <- vector("list", length = 50)
  for (y in seq_along(clim.list)) {
    clim.y <- clim.mo
    clim.y[,"clim.var"] <- clim.y[,"clim.var"] + runif(1, min = 0, max = 4)
    clim.y$year <- y
    clim.list[[y]] <- clim.y
  }
  clim <- do.call("rbind", clim.list)
  expect_vector(n_mon_corr(chrono = chrono, clim = clim, var = "clim.var",
             chrono.col = "std", chrono.name = "Synthetic"))
})

