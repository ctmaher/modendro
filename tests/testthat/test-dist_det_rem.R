test_that("Throws error if not fed a rwl, data.frame, or matrix", {
  expect_error(dist_det_rem(c(1:10)))
  expect_error(dist_det_rem(1))
})

test_that("Throws error if fed data with no rownames (years) or colnames (series IDs)", {
  mat.test <- matrix(nrow = 10, ncol = 10)
  expect_error(dist_det_rem(mat.test))
  rownames(mat.test) <- 1:10
  expect_error(dist_det_rem(mat.test))
  rownames(mat.test) <- NULL
  colnames(mat.test) <- 1:10
  expect_error(dist_det_rem(mat.test))
})

test_that("Throws error if any series are all NAs", {
  mat.test <- matrix(nrow = 10, ncol = 10)
  rownames(mat.test) <- 1:10
  colnames(mat.test) <- 1:10
  mat.test[, 1:5] <- runif(10)
  expect_error(dist_det_rem(mat.test))
})

test_that("Names (series IDs) are equal in the input rwl and output rwls", {
  mat.test <- matrix(nrow = 50, ncol = 10)
  mat.test <- apply(mat.test, MARGIN = 2, FUN = \(x) runif(length(x), 0.1, 2), simplify = TRUE) |>
    as.data.frame()
  rownames(mat.test) <- 1:50
  colnames(mat.test) <- paste0("series", 1:10)
  out.test <- dist_det_rem(mat.test)
  expect_equal(colnames(mat.test), names(out.test[["Corrected RWI"]]))
})

####
test_that("Essential function works - and returns a list", {
  # Make some fake tree ring data
  mat.test <- matrix(nrow = 50, ncol = 10)
  mat.test <- apply(mat.test, MARGIN = 2, FUN = \(x) runif(length(x), 0.1, 2), simplify = TRUE) |>
    as.data.frame()
  rownames(mat.test) <- 1:50
  colnames(mat.test) <- paste0("series", 1:10)
  expect_vector(dist_det_rem(mat.test), ptype = list())
})

#### How about a test that adds specific synthetic disturbances (constructed out of Hugershoff curves) for
# specific start years, and we test to see that those are detected and removed.
# This is a somewhat arbitrary test, as we've set the n trials to 100 and specified just
# one type of underlying curve shape.
test_that("Function works - and detects & removes a specified disturbance most of the time within ±1 year", {
  n.trials <- 100
  dist.year <- vector("numeric", length = n.trials)
  for (i in 1:n.trials) {
  # Make some fake tree ring data
  mat.test <- matrix(nrow = 50, ncol = 10)
  mat.test <- apply(mat.test, MARGIN = 2, FUN = \(x) runif(length(x), -1, 1), simplify = TRUE) |>
    as.data.frame()
  rownames(mat.test) <- 1:50
  colnames(mat.test) <- paste0("series", 1:10)
  # The Hugershoff curve
  hug.fun <- function(x, y) {
    d <- mean(y, na.rm = TRUE)
    dist <- sapply(x, FUN = \(x1){
      (1.5*(x1 - 0)^1)*exp(-0.15*(x1 - 0))
    })

  }

  mat.test[,"series1"] <- (c(rep(0, 40), hug.fun(x = 1:10, y = mat.test[,"series1"])) + mat.test[,"series1"])

  # Turn mat.test into a list
  mat.test.list <-   asplit(mat.test, MARGIN = 2)

  out.test <- dist_det_rem(rwi = mat.test.list, min.win = 6)

  # plot(x = 1:50, y = mat.test[, "series1"], type = "l", col = "grey50")
  # lines(out.test[["Corrected RWI"]][["series1"]], col = "black")
  # lines(x = out.test[["Disturbance curves"]][["series1"]]$year, y = out.test[["Disturbance curves"]][["series1"]]$curve, col = "red")
  # lines(x = as.numeric(names(out.test[["Corrected RWI"]][["series1"]])),
  #       y = c(rep(0, 40), hug.fun(x = 1:10, y = mat.test[,"series1"])), col = "blue")

  if (is.character(out.test$`Disturbance curves`$series1)) {
    dist.year[i] <- NA
  } else {
  dist.year[i] <- out.test$`Disturbance curves`$series1[1, "year"]
  }
  }

  #hist(dist.year)
  #summary(dist.year) # 8 % of runs didn't detect a disturbance
  expect_true(median(dist.year, na.rm = TRUE) < 41 & median(dist.year, na.rm = TRUE) > 39)
  })
