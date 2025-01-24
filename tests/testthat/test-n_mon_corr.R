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
      common.years = 1951:2000,
      agg.fun = "mean",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = "my trees", #
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "mean",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
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
      common.years = 1951:2000,
      agg.fun = "mean",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = "what what", #
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "mean",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
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
      common.years = 1951:2000,
      agg.fun = "mean",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "bad.name", #
      common.years = 1951:2000,
      agg.fun = "mean",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
})


####
test_that("Throws errors if common.years isn't a numeric vector of sufficent length or if
          common years are outside of data coverage", {
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = "the golden years", #
      agg.fun = "mean",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:1953, #
      agg.fun = "mean",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 79:130, #
      agg.fun = "mean",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
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
      common.years = 1951:2000,
      agg.fun = 1, #
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "average", #
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
})

####
test_that("Throws error if not fed a valid max.win", {
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = 1,
      max.win = "6", #
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "average",
      max.win = 50, #
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
})


####
test_that("Throws error if not fed a valid win.align", {
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = 1,
      max.win = "6",
      win.align = "center", #
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "average",
      max.win = 50,
      win.align = 2, #
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
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
      common.years = 1951:2000,
      agg.fun = "sum",
      max.win = 6,
      win.align = "left",
      max.lag = "wayback machine", #
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "sum",
      max.win = 6,
      win.align = "left",
      max.lag = 0, #
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
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
      common.years = 1951:2000,
      agg.fun = "sum",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = 19, #
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "sum",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "the earth is flat", #
      prewhiten = TRUE,
      corr.method = "spearman",
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
      common.years = 1951:2000,
      agg.fun = "sum",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "S",
      prewhiten = 27, #
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "sum",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "S",
      prewhiten = "Nope", #
      corr.method = "spearman",
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
      common.years = 1951:2000,
      agg.fun = "sum",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = 5, #
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "sum",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "Daryl", #
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
})


####
test_that("Throws error if not fed a valid group.IDs.df arg", {
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "sum",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = 2,
      group.var = NULL
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "sum",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = "guess what?",
      group.var = NULL
    )
  )
})

####
test_that("Throws error if not fed a valid group.var arg", {
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "sum",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = 7
    )
  )
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "sum",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = "chicken gut"
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
      common.years = 1951:2000,
      agg.fun = "sum",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )

  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "mean", #
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )

  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1961:2000, #
      agg.fun = "mean",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )


  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1961:2000,
      agg.fun = "mean",
      max.win = 12, #
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )

  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1961:2000,
      agg.fun = "mean",
      max.win = 6,
      win.align = "right", #
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )

  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "mean",
      max.win = 6,
      win.align = "left",
      max.lag = 3, #
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )

  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "mean",
      max.win = 6,
      win.align = "left",
      max.lag = 3,
      hemisphere = "S", #
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )

  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "mean",
      max.win = 6,
      win.align = "left",
      max.lag = 3,
      hemisphere = "S",
      prewhiten = FALSE, #
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )

  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "mean",
      max.win = 6,
      win.align = "left",
      max.lag = 3,
      hemisphere = "S",
      prewhiten = FALSE,
      corr.method = "pearson", #
      group.IDs.df = NULL,
      group.var = NULL
    )
  )

  expect_vector(
    n_mon_corr(
      rwl = rwl,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "mean",
      max.win = 6,
      win.align = "left",
      max.lag = 3,
      hemisphere = "S",
      prewhiten = FALSE,
      corr.method = "kendall", #
      group.IDs.df = NULL,
      group.var = NULL
    )
  )


})

####
test_that("The order of years & months in the input climate data has no effect on the accuracy
of the output (ungrouped data)", {
test.list <- vector("logical", length = 50)
for (i in 1:length(test.list)) {

  mat.test <- matrix(nrow = 50, ncol = 10)
  rownames(mat.test) <- 1:50
  colnames(mat.test) <- paste0("tr",1:10)
  mat.test <- apply(mat.test, MARGIN = 2, FUN = \(x) runif(length(x), 0.1, 2)) |>
    as.data.frame()

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

  # Let's create a marker month, such that it is identical to some of the tree-ring data
  marker.month <- sample(1:12, 1)

  # Insert the tree-ring data from one of the series
  this.tr.data <- mat.test[, sample(colnames(mat.test), 1)]
  clim$clim.var[clim$month %in% marker.month] <-
    this.tr.data

  # Visual check
  # ggplot(clim[clim$month %in% marker.month,], aes(year, clim.var)) +
  #   geom_line() +
  #   geom_line(data = data.frame(year = unique(clim$year),
  #                               clim.var = this.tr.data),
  #             color = "red")

  # Now scramble the order of the years & months
  clim <- clim[sample(1:nrow(clim), nrow(clim)),]

  # Run the function
  test10 <- n_mon_corr(
    rwl = mat.test,
    clim = clim,
    clim.var = "clim.var",
    common.years = 1:50,
    agg.fun = "mean",
    max.win = 6,
    win.align = "left",
    max.lag = 1,
    hemisphere = "N",
    prewhiten = FALSE,
    corr.method = "spearman",
    group.IDs.df = NULL,
    group.var = NULL
  )

  cor.res <- test10[["Correlation results"]]

  # The top (or first) row is the highest correlation (output is sorted),
  # thus marker month should be the 1st row
  test.list[[i]] <- cor.res$month[cor.res$win.len %in% 1 &
                                    cor.res$lag %in% "0"][1] %in% marker.month
}
expect_true(all(test.list))
})


####
test_that("The order of years & months in the input climate data has no effect on the accuracy
of the output (grouped data)", {
  test.list <- vector("logical", length = 10)
  for (i in 1:length(test.list)) {

    mat.test <- matrix(nrow = 50, ncol = 10)
    rownames(mat.test) <- 1:50
    colnames(mat.test) <- paste0("tr",1:10)
    mat.test <- apply(mat.test, MARGIN = 2, FUN = \(x) runif(length(x), 0.1, 2)) |> as.data.frame()

    # Make some fake climate data
    mo <- 1:12
    x <- seq(-4, 4, length.out = 12)
    gauss.curv <- \(x) {(10/sqrt(2*pi*1.6))*exp(-((x^2)/(2*1.6^2)))}
    clim.mo <- data.frame(month = mo, clim.var = gauss.curv(x))

    clim.list <- vector("list", length = 12)
    for (y in seq_along(clim.list)) {
      clim.y <- data.frame(year = 1:50, month = y)
      clim.y[,"clim.var"] <- runif(50, min = 0.1, max = 2) +
        clim.mo[clim.mo$month %in% y, "clim.var"]
      clim.list[[y]] <- clim.y
    }
    clim <- do.call("rbind", clim.list)
    clim <- clim[order(clim$year, clim$month),]

    # Make groups by replicating climate data and giving different labels
    clim <- rbind(clim, clim)
    clim$group <- c(rep("g1", 600), rep("g2", 600))

    # We'll need a group.IDs.df too
    group.IDs.df <- data.frame(series = colnames(mat.test),
                               group = c(rep("g1", 5), rep("g2", 5)))

    # Let's create a marker month, such that it is identical to some of the tree-ring data
    marker.month <- sample(1:12, 1)

    # Insert the tree-ring data from one of the g1 series into the g1 climate data
    this.tr.data <- mat.test[, sample(colnames(mat.test)[1:5], 1)]
    clim$clim.var[clim$month %in% marker.month &
                                               clim$group %in% "g1"] <-
      this.tr.data

    # Visual check
    # ggplot(clim[clim$month %in% marker.month,], aes(year, clim.var)) +
    #   geom_line() +
    #   geom_line(data = data.frame(year = unique(clim$year),
    #                               clim.var = this.tr.data),
    #             color = "red") +
    #   facet_wrap( ~ group, nrow = 2)

    # Now scramble the order of the years & months
    clim <- clim[sample(1:nrow(clim), nrow(clim)),]

    # Run the function
    test10 <- n_mon_corr(
      rwl = mat.test,
      clim = clim,
      clim.var = "clim.var",
      common.years = 1:50,
      agg.fun = "mean",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = FALSE,
      corr.method = "spearman",
      group.IDs.df = group.IDs.df,
      group.var = "group"
    )

    cor.res <- test10[["Correlation results"]]
    # The top (or first) row is the highest correlation (output is sorted),
    # thus marker month should be the 1st row
    # We also know that the signal should be for win.len = 1 & in group = "g1"
    test.list[[i]] <- cor.res$month[cor.res$win.len == 1 &
                                      cor.res$lag %in% "0"][1] %in% marker.month &
      cor.res$series[cor.res$win.len == 1 &
                       cor.res$lag %in% "0"][1] %in%
      group.IDs.df[group.IDs.df$group %in% "g1", "series"]
  }
  expect_true(all(test.list))
})


#### Missing data in clim
test_that("Throws an error when clim.var has internal NAs",{
  clim.miss <- clim
  clim.miss[clim.miss$year %in% sample(1951:2000, 1) &
                           clim.miss$month %in% sample(1:12, 1), "clim.var"] <- NA
expect_error(
  n_mon_corr(
    rwl = rwl,
    clim = clim.miss,
    clim.var = "clim.var",
    common.years = 1951:2000,
    agg.fun = "sum",
    max.win = 6,
    win.align = "left",
    max.lag = 1,
    hemisphere = "N",
    prewhiten = TRUE,
    corr.method = "spearman",
    group.IDs.df = NULL,
    group.var = NULL
  )
)
})

test_that("Throws an error when clim.var is missing observations",{
  clim.miss <- clim[!(clim$year %in% sample(1951:2000, 1) &
                        clim$month %in% sample(1:12, 1)),]
  expect_error(
    n_mon_corr(
      rwl = rwl,
      clim = clim.miss,
      clim.var = "clim.var",
      common.years = 1951:2000,
      agg.fun = "sum",
      max.win = 6,
      win.align = "left",
      max.lag = 1,
      hemisphere = "N",
      prewhiten = TRUE,
      corr.method = "spearman",
      group.IDs.df = NULL,
      group.var = NULL
    )
  )
})


#### Check that correlation coefficients are correct with an independent check

test_that("Correlation coefficients are equal with an independent check", {
  # Get the moving averages of the climate data
  clim.ma <- moving_win_multi(clim,
                              clim.var = "clim.var",
                              win.lens = 1:2,
                              win.align = "left",
                              agg.fun = "mean")

  # Set common.years here
  common.years <- 1960:2000

  rwl.comm <- rwl[rownames(rwl) %in% common.years,]
  # Disassemble the rwl and run correlations with each climate series
  suppressWarnings(
  ind.test <- mapply(FUN = \(tr, tr.names) {
    lapply(split(clim.ma[clim.ma$win.len == 2,],
                 f = clim.ma[clim.ma$win.len == 2, "month"]), FUN = \(mon) {
                   mon <- mon[order(mon[,"year"], decreasing = FALSE),]
                   mon <- mon[mon[,"year"] %in% common.years,]
                   cor.res <- corTESTsrd(x = mon[, "clim.var"],
                              y = tr,
                              method = "spearman",
                              iid = FALSE,
                              alternative = "two.sided")

                   data.frame(series = tr.names,
                              month = unique(mon[,"month"]),
                              ind.coef = cor.res[["rho"]])
                 }) |> do.call(what = "rbind")

  }, tr = as.list(rwl.comm), tr.names = colnames(rwl.comm),
  SIMPLIFY = FALSE) |> do.call(what = "rbind")
  )

  fun.test <- n_mon_corr(rwl = rwl,
                         clim = clim,
                         clim.var = "clim.var",
                         common.years = common.years,
                         agg.fun = "mean",
                         max.win = 2,
                         win.align = "left",
                         max.lag = 1,
                         hemisphere = "N",
                         prewhiten = FALSE,
                         corr.method = "spearman")
  fun.test.res <- fun.test[["Correlation results"]]
  # Subset to match the 2-month moving window and current year
  fun.test.res <- fun.test.res[fun.test.res$win.len == 2 & fun.test.res$lag %in% "0",]
  test.merge <- merge(fun.test.res, ind.test, by = c("series","month"))

  expect_equal(test.merge$coef, test.merge$ind.coef)

  })

