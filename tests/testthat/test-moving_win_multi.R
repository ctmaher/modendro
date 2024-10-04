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


####
test_that("The moving window means match outputs from 2 common packages", {

  func.val <- moving_win_multi(df = clim,
                               clim.var = "clim.var",
                               win.lens = 2:12,
                               agg.fun = "mean")

  sapply(2:12, FUN = \(this.win) {

  forecast.val <- forecast::ma(x = clim$clim.var,
                               order = this.win,
                               centre = FALSE)

  DescTools.val <- DescTools::MoveAvg(x = clim$clim.var,
                               order = this.win,
                               align = "left")

  expect_equal(na.omit(func.val[func.val$win.len == this.win, "clim.var"]) |> as.numeric(),
               na.omit(as.numeric(forecast.val)) |> as.numeric())

  expect_equal(na.omit(func.val[func.val$win.len == this.win, "clim.var"]) |> as.numeric(),
               na.omit(DescTools.val) |> as.numeric())
})
})


####
test_that("The moving window means are labeled correctly in the output", {

  func.val <- moving_win_multi(df = clim,
                               clim.var = "clim.var",
                               win.lens = 2:12,
                               agg.fun = "mean")

  lapply(split(func.val, f = func.val$win.len), FUN = \(win) {
    this.win.len <- unique(as.numeric(win[,"win.len"]))
    expect_equal(this.win.len, nrow(win[is.na(win[,"clim.var"]),])+1)
  })
})
