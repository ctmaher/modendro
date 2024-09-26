#' Calculate multiple-length moving windows of monthly climate data
#'
#' @description
#' This function is a workhorse function behind \code{\link{n_mon_corr}}. It very efficiently
#' calculates the moving window sums or means by taking lagged differences from the cumsum() of
#' the properly ordered climate data.
#'
#' @param df a data.frame or matrix of climate data, with at least 3 columns: year, month (numeric),
#' and a climate variable.
#' @param clim.var a character vector of length = 1 that specifies the column name in `df` that
#' refers to the climate variable you wish to calculate moving windows for.
#' @param win.lens a numeric sequence specifying which moving windows (in months) to calculate. In
#' \code{\link{n_mon_corr}}, this will be 2:12.
#' @param agg.fun character vector of length = 1 that specifies whether to take the mean (e.g., for
#' temperature) or sum (e.g., for precipitation) of the moving windows. Options are
#' `c("mean","sum")`.
#'
#' @return A data.frame with the same number of rows as `df` - essentially it is `df` with the
#' moving windows attached with \code{\link[base]{cbind}}. For each of the moving windows, the
#' numeric month variable indicates the start month of the window. The length of the window is
#' indicated in the column name (e.g., "mean2", "mean6", "mean12", etc.).
#'
#' @export
#'
#' @examples
#' # Make some fake climate data for testing and building the function
#' mo <- 1:12
#' x <- seq(-4, 4, length.out = 12)
#' # Give the clim.var a seasonal fluctuation, kind of like temperature outside of
#' # tropical latitudes
#' gauss.curv <- \(x) {(10/sqrt(2*pi*1.6))*exp(-((x^2)/(2*1.6^2)))}
#' clim.mo <- data.frame(month = mo, clim.var = gauss.curv(x))
#'
#' clim.list <- vector("list", length = 12)
#' for (y in seq_along(clim.list)) {
#'   clim.y <- data.frame(year = 1951:2000, month = y)
#'   clim.y[,"clim.var"] <- runif(50, min = 0.1, max = 2) +
#'   clim.mo[clim.mo$month %in% y, "clim.var"]
#'   clim.list[[y]] <- clim.y
#' }
#'
#' fake.clim <- do.call("rbind", clim.list)
#'
#' test.output <- moving_win_multi(df = fake.clim,
#' clim.var = "clim.var",
#' win.lens = 2:12,
#' agg.fun = "mean")
#'
#' # Make a plot of 2-, 6-, and 12-month moving window means.
#' library(ggplot2)
#'
#' # Add a date column to make plotting easier
#' test.output$yr.mo <- as.Date(paste(test.output$year,
#' ifelse(nchar(test.output$month) == 1,
#'        paste0("0", test.output$month),
#'        test.output$month), "01",
#' sep = "-"), format = "%Y-%m-%d")
#'
#' # The continuous moving windows - note that the month on the x-axis refers to the start month
#' # of each moving window.
#' # Trim down to a 10-year period so that the results are easier to see.
#'
#' ggplot(clim1[clim1$year %in% 1991:2000,]) +
#' geom_path(aes(yr.mo, clim.var, color = win.len), na.rm = TRUE)


moving_win_multi <- function(df, clim.var, win.lens, agg.fun = NULL) {
  ## Initial error catching and interactive prompts
  stopifnot(
    "Arg df is not an object of class 'data.frame', or 'matrix'" =
      data.class(df) %in% "data.frame" |
      data.class(df) %in% "matrix"
  )

  stopifnot(
    "df does not have a year column? (colname sould start with 'y' or 'Y')" =
      any(substr(colnames(df), 1, 1) %in% c("Y", "y")) == TRUE
  )

  stopifnot(
    "df does not have a month column? (colname sould start with 'm' or 'M')" =
      any(substr(colnames(df), 1, 1) %in% c("M", "m")) == TRUE
  )

  match.test <- clim.var %in% colnames(df)
  stopifnot("Arg clim.var must match one unique column name in df" =
              length(match.test[match.test == TRUE]) == 1)


  stopifnot("Arg agg.fun must be either 'mean' or 'sum'" =
              agg.fun %in% "mean" |
              agg.fun %in% "sum")

  # make sure that "year" columns are labelled as such - ie standardize the names
  colnames(df)[which((substr(
    colnames(df), start = 1, stop = 1
  )
  %in% c("Y", "y")) == T)] <- "year"

  # same for month
  colnames(df)[which((substr(
    colnames(df), start = 1, stop = 1
  )
  %in% c("M", "m")) == T)] <- "month"

  stopifnot(
    "Year variable in climate data not numeric or integer" =
      is.numeric(df[, "year"]) |
      is.integer(df[, "year"]) |
      is.double(df[, "year"])
  )

  stopifnot(
    "Month variable in climate data not numeric or integer" =
      is.numeric(df[, "month"]) |
      is.integer(df[, "month"]) |
      is.double(df[, "month"])
  )

  # make sure year and month are integers from here on out
  df[, "year"] <- as.integer(df[, "year"])
  df[, "month"] <- as.integer(df[, "month"])

  # must ensure that the years and months are ordered correctly
  df <- df[order(df$year, df$month, decreasing = FALSE), ]

  # moving_win_multi assumes continuity and regularity of all years and months.
  # If even one month or year is missing somewhere, this will mess up everything that follows.

  all.yrs <- df[, "year"] |> unique()
  stopifnot("One or more years are missing in the climate data df" =
              length(all.yrs) == length(min(all.yrs):max(all.yrs)))

  mo.count <- aggregate(month ~ year, data = df, length)
  stopifnot("One or more years are missing months in the climate data df" =
              all(mo.count[, "month"] == 12))

  ## Run through the meat of the function
  # Run cumsum() on the clim.var
  cs <- cumsum(c(0, df[, clim.var])) # all windows start with cumsum(), so just do this once

  # the summed moving windows are just lagged differences of the cumsum
  all.windows <- lapply(win.lens, FUN = \(n) {
    result <- cs[(n + 1):length(cs)] - cs[1:(length(cs) - n)]
    this.n <- data.frame(
      year = df[, "year"],
      start.month = df[, "month"],
      agg.fun = agg.fun,
      win.len = n,
      result = c(result, rep(NA, n - 1))
    )
    colnames(this.n)[which(colnames(this.n) %in% "result")] <- clim.var
    this.n
  }) |> do.call(what = "rbind")

  # the averaged moving windows are just lagged differences of the cumsum divided by the
  # window length
  if (agg.fun == "mean") {
    all.windows <- lapply(win.lens, FUN = \(n) {
      result <- cs[(n + 1):length(cs)] - cs[1:(length(cs) - n)]
      this.n <- data.frame(
        year = df[, "year"],
        start.month = df[, "month"],
        agg.fun = agg.fun,
        win.len = n,
        result = c(result, rep(NA, n - 1)) / n
      )
      colnames(this.n)[which(colnames(this.n) %in% "result")] <- clim.var
      this.n
    }) |> do.call(what = "rbind")
  }
  # # attach the moving windows to the original data.frame containing the climate variable
  orig.data <- data.frame(
    year = df[, "year"],
    start.month = df[, "month"],
    agg.fun = agg.fun,
    win.len = 1,
    result = df[, clim.var]
  )
  colnames(orig.data)[which(colnames(orig.data) %in% "result")] <- clim.var

  output.data <- rbind(orig.data, all.windows)
  output.data[, "win.len"] <- factor(output.data[, "win.len"], levels = 1:max(win.lens))
  output.data
}
