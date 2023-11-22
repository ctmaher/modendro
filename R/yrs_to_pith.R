#' Estimate years to pith using estimated distance to pith and innermost measured rings of series
#'
#' @description
#' Function to estimate the number of years to pith from estimated distance to pith, like that determined in CooRecorder. Other methods are fine too,
#' as long as the input data.frame is formatted in the same way (see `d2pith` argument below). The simple method (the default) uses the arithmetic mean of n.rings
#' to apply over the estimated distance, rounded to the nearest whole year.
#'
#' @param rwl A rwl-type data.frame (e.g., read in by \code{\link[dplR]{read.rwl}}). Essentially a data.frame with columns names as series IDs and years as rownames.
#' @param d2pith A data.frame with at least 2 columns: 1) "series", representing the series names, 2) "d2pith", representing the estimated distances to pith from innermost measured rings.
#' If `method = "modeled"` , a 3rd "pith.offset" column is required where this is an initial estimate of years to pith.
#' This is only used to identify those series that had pith in the core (taken as those with pith.offset = 1).
#' @param method Character vector of length 1 indicating which method to use, "simple", or "modeled". Default is "simple". See details below.
#' @param n.rings Numeric vector of length 1 representing the number of rings (years) you want to use as an aggregate. For `method = "simple"` only.
#' @param diam.df A data.frame containing series names and the outside diameter measurements associated with the tree those series belong to. For `method = "modeled"` only.
#'
#' @details
#' There are a number of methods to estimate distance to and years to pith in core samples that don't contain pith,
#' from the old transparency method (citation?) to the Duncan (1989) method. CooRecorder (https://www.cybis.se/forfun/dendro/helpcoorecorder7/distanceToPith/index.htm)
#' has it's own intuitive method that is like a flexible version of the transparency method. CooRecorder has the additional benefit that
#' distance to pith estimates are essentially seamless with ring measurements. Estimating years to pith in CooRecorder, however, is somewhat cumbersome.
#' This function is designed to save time by doing the calculations necessary for each series all at once.
#'
#' The "modeled" method is a somewhat experimental method that uses a generalized additive model (GAM) and a number of
#' series from your site that have pith.offset = 1 (i.e., from cores that contained pith). More is better, and hopefully
#' coverage of the diameter range is adequate to predict for the series that don't contain pith.
#'
#'
#' @return A data.frame with 3 columns: 1) "series", 2) "d2pith", 3) "pith.offset".
#'
#' @references
#' Duncan, R. P. (1989). An Evaluation of Errors in Tree Ring Age Estimates Based on Increment Cores in Kahikatea (Dacrycarpus dacrydioides). \emph{New Zealand Natural Sciences} \strong{16}, 31–37.
#'
#'
#' @import stats
#' @import mgcv
#' @import dplR
#' @export
#'
#' @examples
#' # Examples to come
#'
#'
#'

yrs_to_pith <- function(rwl = NULL,
                        d2pith = NULL,
                        method = "simple",
                        n.rings = 5,
                        diam.df = NULL) {
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
    "rwl has no rownames (must be years only) or no colnames (must be series IDs only)" =
      !is.null(rownames(rwl)) |
      !is.null(colnames(rwl))
  )

  stopifnot(
    "rwl and d2pith series names do not completely match" =
      all((colnames(rwl) %in% d2pith$series) == TRUE)
  )



  xdf <- rwl_longer(rwl,
                    omit.NAs = TRUE)
  # Split by series
  xdf.list <- split(xdf, f = xdf$series)

  if (method %in% "simple") {
  # Get aggregates of n.rings for each series
  n.rings.agg <- lapply(xdf.list, FUN = \(z) {
  z <- z[order(z$year),] # Make sure the order is correct
  mean(z[1:n.rings, "rw"], na.rm = TRUE)
  }) |> do.call(what = "rbind") |> as.data.frame()
  colnames(n.rings.agg) <- "mean.rw"
  n.rings.agg$series <- rownames(n.rings.agg)

  # merge the attributes and the mean.rw by series
  merged.att <- merge(d2pith, n.rings.agg, by = "series")
  merged.att$y2pith <- (merged.att$d2pith / merged.att$mean.rw) |> round(digits = 0)

  merged.att

  } else {
    message("method 'modeled' not available yet")
  }
}
