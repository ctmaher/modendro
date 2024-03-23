####
test_that("Throws error if not fed a chrono, data.frame, or matrix", {
  expect_error(n_mon_corr(rw = c(1:10), clim = c(1:10), rel.per.begin = 10, hemisphere = "N"))
  expect_error(n_mon_corr(rw = 1, clim = 1, rel.per.begin = 10, hemisphere = "N"))
})
####
test_that("Throws error if not fed a valid clim var name", {
  clim <- matrix(nrow = 10, ncol = 3)
  colnames(clim) <- c("year", "month", "clim.var")
  expect_error(n_mon_corr(rw = matrix(nrow = 10, ncol = 2), clim = clim,
                          rel.per.begin = 10, hemisphere = "N",
                          clim.var = "wack-a-doodle"))
})
####
test_that("Throws error if not fed a valid chrono.col name", {
  clim <- matrix(nrow = 10, ncol = 3)
  colnames(clim) <- c("year", "month", "clim.var")
  chrono <- matrix(nrow = 10, ncol = 2)
  colnames(chrono) <- c("samp.depth", "std")
  expect_error(n_mon_corr(rw = chrono, clim = clim, rel.per.begin = 10, hemisphere = "N",
                          clim.var = "clim.var", rw.col = "wack-a-doodle"))
})
####
test_that("Throws error if not fed a valid aggregation function", {
  expect_error(n_mon_corr(rw = matrix(nrow = 10, ncol = 2), clim = matrix(nrow = 10, ncol = 3),
                          rel.per.begin = 10, hemisphere = "N",
                          agg.fun = "wack-a-doodle"))
})
####
test_that("Throws error if not fed a valid correlation method", {
  clim <- matrix(nrow = 10, ncol = 3)
  colnames(clim) <- c("year", "month", "clim.var")
  chrono <- matrix(nrow = 10, ncol = 2)
  colnames(chrono) <- c("samp.depth", "std")
  expect_error(n_mon_corr(rw = chrono, clim = clim,
                          rel.per.begin = 10, hemisphere = "N",
                          clim.var = "clim.var", rw.col = "std",
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
  expect_vector(n_mon_corr(rw = chrono, clim = clim,
                           rel.per.begin = 3, hemisphere = "S", clim.var = "clim.var",
                           rw.col = "std", rw.name = "Synthetic", silent = TRUE))
})
####
test_that("The order of months in the input climate data has no effect on the accuracy of the
          output", {
test.list <- vector("logical", length = 50)
for (i in 1:length(test.list)) {

  mat.test <- matrix(nrow = 50, ncol = 10)
  rownames(mat.test) <- 1:50
  colnames(mat.test) <- 1:10
  mat.test <- apply(mat.test, MARGIN = 2, FUN = \(x) runif(length(x), 0.1, 2))
  chrono <- dplR::chron(mat.test)
  chrono$year <- rownames(chrono) |> as.numeric()
  # Make some fake climate data
  mo <- 1:12
  x <- seq(-4, 4, length.out = 12)
  gauss.curv <- \(x) {(10/sqrt(2*pi*1.6))*exp(-((x^2)/(2*1.6^2)))}
  clim.mo <- data.frame(month = mo, clim.var = gauss.curv(x))

  clim.list <- vector("list", length = 12)
  for (y in seq_along(clim.list)) {
    clim.y <- data.frame(year = 1:50, month = y)
    clim.y[,"clim.var"] <- runif(50, min = 0.1, max = 2) + clim.mo[clim.mo$month %in% y, "clim.var"]
    clim.list[[y]] <- clim.y
  }
  clim <- do.call("rbind", clim.list)
  clim <- clim[order(clim$year, clim$month),]

  # Let's create a marker month, such that it is identical to the chronology
  marker.month <- sample(1:12, 1)
  # Actually have to account for the grow year concept here, used in n_mon_corr
  # I.e., if rel.per.begin = 10, then months 10-12 of the previous calendar year
  # are part of the grow year, but that will confuse the signal here.
  # Since I use rel.per.begin = 10 below, let's plan for this
  if(marker.month %in% c(10, 11, 12)){
    clim[order(clim$year), "clim.var"][clim$month %in% marker.month] <-
      c(chrono[order(chrono$year),"std"][2:50],NA)
  } else {
    clim[order(clim$year), "clim.var"][clim$month %in% marker.month] <-
      chrono[order(chrono$year),"std"]
  }

  # ggplot(clim[clim$month %in% marker.month,], aes(year, clim.var)) +
  #   geom_line() +
  #   facet_wrap(~month)
  #
  # ggplot(chrono, aes(year, std)) +
  #   geom_line()

  test10 <- n_mon_corr(rw = chrono, clim = clim, clim.var = "clim.var",
                       rel.per.begin = 10, hemisphere = "N", rw.name = "Synthetic",
                       corr.method = "pearson", silent = TRUE)
  cor.res <- test10[["Correlation results"]]
  #
  test.list[[i]] <- cor.res$months[1] %in% marker.month
  # The top row is the highest correlation (output is sorted),
  # thus marker month should be the 1st row
}
expect_true(all(test.list))
})

