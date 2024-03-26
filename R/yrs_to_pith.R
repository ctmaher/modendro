#' Estimate years to pith using estimated distance to pith and innermost measured rings of series
#'
#' @description
#' Function to estimate the number of years to pith from estimated distance to pith, like that
#' determined in CooRecorder. The way you determine distance to pith is not important, as long as
#' the input data.frame is formatted in the same way (see `d2pith` argument below). The `"rings"`
#' method (the default) uses the arithmetic mean of n.rings to apply over the estimated
#' distance, rounded to the nearest whole year.
#'
#' @param rwl A rwl-type data.frame (e.g., read in by \code{\link[dplR]{read.rwl}}). Essentially a
#' data.frame with columns names as series IDs and years as row names.
#' @param d2pith A data.frame with at least 2 columns: 1) "series", representing the series names,
#' 2) "d2pith", representing the estimated distances to pith from innermost measured rings.
#' @param method Character vector of length 1 indicating which method to use, `"rings"`, or
#' `"dist"`. Default is "rings". See details below.
#' @param n.rings Numeric vector of length 1 representing the number of rings (years) you want to
#' use as an aggregate. For `method = "rings"` only.
#' @param diam.df A data.frame containing series names and the outside diameter measurements
#' associated with the tree those series belong to. For `method = "modeled"` only.
#' @param plot.hist Logical vector that turns on or off a histogram of the estimated years to pith.
#' This is an important guide to choosing the method or setting `n.rings` for method `"rings"`. See
#' details below.
#'
#' @details
#' There are a number of methods to estimate distance to and years to pith in core samples that
#' don't contain pith, from the old transparency method (citation?) to the Duncan (1989) method. I
#' assume that the user has some familiarity with these methods before using `yrs_to_pith()`.
#' CooRecorder (https://www.cybis.se/forfun/dendro/helpcoorecorder7/distanceToPith/index.htm)
#' has it's own intuitive method that is like a flexible version of the transparency method.
#' CooRecorder has the additional benefit that distance to pith estimates are essentially seamless
#' with ring measurements. Estimating years to pith in CooRecorder, however, is somewhat cumbersome.
#' This function is designed to save time by doing the calculations necessary for each series all
#' at once.
#'
#' Method "rings" is what most researchers will be familiar with - this is simply the mean of the
#' specified innermost `n.rings`. The default is 5 rings, but this may be too few for some samples
#' or datasets. Use the 'plot.hist' output to check this - if you have some extreme values (large
#' values of y2pith), you may need to increase `n.rings` or switch to method "dist". Extreme y2pith
#' values can arise from large d2pith values and/or if the innermost n.rings happen to be very
#' narrow. The take home message is that the value of `n.rings` is not allowing a robust enough
#' estimate of average ring width to apply over d2pith. Method "dist" attempts to circumvent the
#' estimation of a "Goldilocks" `n.rings` value by using instead the number of rings that add up to
#' d2pith instead (`n.rings` is variable). I.e., we use the mean ring width over the innermost
#' \emph{distance} equal to d2pith. Defaults to the mean of 5 innermost rings if d2pith yields less
#' than 5 rings.
#'
#'
#' @return A data.frame with 3 columns: 1) "series", 2) "d2pith", 3) "y2pith".
#'
#' @references
#' Duncan, R. P. (1989). An Evaluation of Errors in Tree Ring Age Estimates Based on Increment
#' Cores in Kahikatea (Dacrycarpus dacrydioides). \emph{New Zealand Natural Sciences} \strong{16},
#' 31–37.
#'
#'
#' @import stats
#' @import mgcv
#' @import dplR
#' @import ggplot2
#' @export
#'
#' @examples
#' # Examples to come
#'
#'
#'

yrs_to_pith <- function(rwl = NULL,
                        d2pith = NULL,
                        method = "rings",
                        n.rings = 5,
                        diam.df = NULL,
                        plot.hist = TRUE) {
  ## Error catching & warnings

  #
  stopifnot(
    "rwl is not an object of class 'rwl', 'data.frame', or 'matrix'" =
      data.class(rwl) %in% "rwl" |
      data.class(rwl) %in% "data.frame" |
      data.class(rwl) %in% "matrix"
  )

  #
  stopifnot(
    "d2pith is not an object of class 'data.frame', or 'matrix'" =
      data.class(d2pith) %in% "data.frame" |
      data.class(d2pith) %in% "matrix"
  )

  #
  stopifnot(
    "d2pith must have columns 'series' & 'd2pith'" =
      any(colnames(d2pith) %in% "series") |
      any(colnames(d2pith) %in% "d2pith")
  )

  #
  stopifnot(
    "rwl has no rownames (must be years only) or no colnames (must be series IDs only)" =
      !is.null(rownames(rwl)) |
      !is.null(colnames(rwl))
  )

  #
  stopifnot(
    "method must be either 'rings' or 'dist'" =
      method %in% "rings" |
      method %in% "dist"
  )

  #
  stopifnot(
    "rwl and d2pith series names do not completely match" =
      all((colnames(rwl) %in% d2pith$series) == TRUE)
  )

  #
  stopifnot(
    "n.rings must be a numeric vector of length 1" =
      data.class(n.rings) == "numeric" |
      data.class(n.rings) == "double" &
      length(n.rings) == 1
  )


  xdf <- rwl_longer(rwl,
                    series.name = "series",
                    dat.name = "rw",
                    trim = TRUE)
  # Split by series
  xdf.list <- split(xdf, f = xdf$series)

  if (method %in% "rings") {
  # Get aggregates of n.rings for each series
  n.rings.agg <- lapply(xdf.list, FUN = \(z) {
  z <- z[order(z$year),] # Make sure the order is correct
  mean(z[1:n.rings, "rw"], na.rm = TRUE)
  }) |> do.call(what = "rbind") |> as.data.frame()
  colnames(n.rings.agg) <- paste("mean.rw", n.rings, sep = ".")
  n.rings.agg$series <- rownames(n.rings.agg)

  # merge the attributes and the mean.rw by series
  merged.att <- merge(d2pith, n.rings.agg, by = "series")
  merged.att$y2pith <- (merged.att$d2pith / merged.att[, paste("mean.rw", n.rings, sep = ".")]) |>
    round(digits = 0)

  merged.att

  } else { # for method == 'dist'
    # Find how many rings add up to the missing d2pith
    dist.agg <- lapply(xdf.list, FUN = \(z) {
      # Make sure the order is correct
      z <- z[order(z$year),]
      # Get this series's d2pith
      this.series.d2pith <- d2pith[d2pith$series %in% z[,"series"],]
      # calculate the cumulative ring width
      z$cum.rw <- cumsum(z[,"rw"])
      # Control for small d2pith
      these.rings <- z[z$cum.rw < this.series.d2pith$d2pith, "rw"]
      if (length(these.rings) < 5) {
        these.rings <- z[1:5, "rw"] # default to innermost 5 rings if dist yields fewer
      }
      # Get the mean of all the ring widths while the cumulative width is less than d2pith
      mean(these.rings, na.rm = TRUE)
    }) |> do.call(what = "rbind") |> as.data.frame()
    colnames(dist.agg) <- "mean.rw"
    dist.agg$series <- rownames(dist.agg)

    # merge the attributes and the mean.rw by series
    merged.att <- merge(d2pith, dist.agg, by = "series")
    merged.att$y2pith <- (merged.att$d2pith / merged.att$mean.rw) |> round(digits = 0)

    merged.att
  }

  # Make a histogram
  if (plot.hist == TRUE) {
    ggplot(merged.att, aes(y = y2pith)) +
      geom_histogram(binwidth = 5) +
      xlab(paste("Estimated years to pith;n/",
                 ifelse(method %in% "rings",
                        "method = 'rings', mean of innermost", n.rings, "rings",
                        "method = 'dist'")))
  }

}
