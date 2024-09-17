#### Some set up vectors
# Make some fake tree ring data
rwl <- matrix(nrow = 50, ncol = 10)
rwl <- apply(rwl, MARGIN = 2, FUN = \(x) runif(length(x), 0.1, 2)) |>
  as.data.frame()
colnames(rwl) <- paste0("t", 1:10)
rownames(rwl) <- 1951:2000
# Make some fake climate data
mo <- 1:12
x <- seq(-4, 4, length.out = 12)
gauss.curv <- \(x) {(10/sqrt(2*pi*1.6))*exp(-((x^2)/(2*1.6^2)))}
clim.mo <- data.frame(month = mo, clim.var = gauss.curv(x))

clim.list <- vector("list", length = 12)
for (y in seq_along(clim.list)) {
  clim.y <- data.frame(year = 1951:2000, month = y)
  clim.y[,"clim.var"] <- runif(50, min = 0.1, max = 2) + clim.mo[clim.mo$month %in% y, "clim.var"]
  clim.list[[y]] <- clim.y
}
clim <- do.call("rbind", clim.list)

clim <- clim[order(clim$year, clim$month),]
# These will allow testing of individual elements


####
test_that("Throws error if rwl is not an rwl or data.frame", {
  expect_error(
    n_mon_corr(
      rwl = 1:10, #
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "mean",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = "what what", #
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "mean",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
})

####
test_that("Throws error if clim is not a data.frame or matrix", {
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = 1:10, #
      clim.var = "clim.var",
      agg.fun = "mean",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = "what what", #
      clim.var = "clim.var",
      agg.fun = "mean",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
})
####
test_that("Throws errors if clim.var isn't a character or doesn't match a column in clim", {
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = runif(50), #
      agg.fun = "mean",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "bad.name", #
      agg.fun = "mean",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
})

####
test_that("Throws error if not fed a valid agg.fun", {
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = 1, #
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "average", #
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
})

####
test_that("Throws error if not fed a valid max.lag", {
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "sum",
      max.lag = "wayback machine", #
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "sum",
      max.lag = 0, #
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
})

####
test_that("Throws error if not fed a valid hemisphere", {
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "sum",
      max.lag = 1,
      hemisphere = 19, #
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "sum",
      max.lag = 1,
      hemisphere = "the earth is flat", #
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
})

####
test_that("Throws error if not fed a valid prewhiten arg", {
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "sum",
      max.lag = 1,
      hemisphere = "S",
      prewhiten = 27, #
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "sum",
      max.lag = 1,
      hemisphere = "S",
      prewhiten = "Nope", #
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
})

####
test_that("Throws error if not fed a valid correlation method", {
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "sum",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = 5, #
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "sum",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "Daryl's", #
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
})

####
test_that("Throws error if not fed a valid gro.period.end", {
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "sum",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 13, #
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "sum",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = c(1:2), #
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "sum",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = "March", #
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "sum",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 4.5, #
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
})

####
test_that("Throws error if not fed a valid make.plot arg", {
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "sum",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = 5,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "sum",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = "yessiree",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
})

####
test_that("Essential function works - and returns a vector under all basic variations of args", {
  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "sum",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )

  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "mean", #
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )

  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "mean", #
      max.lag = 3, #
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )

  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "mean", #
      max.lag = 3, #
      hemisphere = "S", #
      prewhiten = TRUE,
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )

  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "mean", #
      max.lag = 3, #
      hemisphere = "S", #
      prewhiten = FALSE, #
      corr.method = "spearman",
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )

  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "mean", #
      max.lag = 3, #
      hemisphere = "S", #
      prewhiten = FALSE, #
      corr.method = "pearson", #
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )

  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "mean", #
      max.lag = 3, #
      hemisphere = "S", #
      prewhiten = FALSE, #
      corr.method = "kendall", #
      gro.period.end = 9,
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )

  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "mean", #
      max.lag = 3, #
      hemisphere = "S", #
      prewhiten = FALSE, #
      corr.method = "kendall", #
      gro.period.end = 4, #
      make.plot = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )

  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      agg.fun = "mean", #
      max.lag = 3, #
      hemisphere = "S", #
      prewhiten = FALSE, #
      corr.method = "kendall", #
      gro.period.end = 4, #
      make.plot = FALSE, #
      group.IDs.df = NULL,
      group.var = NULL
    )
  )

})

####### Left off here 17 Sep

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

