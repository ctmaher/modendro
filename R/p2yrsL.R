#' Highlight interannual variability in tree ring series by computing ring width as the proportion
#' of the last 2 years, limited by a variance cut-off
#'
#' @description
#' This function implements Lars-Ake Larsson's (author of the popular CooRecorder and CDendro
#' programs) proportion of the last two years (limited) P2YrsL standardizing. This is calculated
#' as i = rw(c)/rw(c) + rw(c-1), with a cut off of 2.6x the SD for the whole series
#' (for details see: https://cdendro.se/wiki/index.php/Normalization).
#' The cut off is imposed to reduce (not eliminate!) the chances that you will get spurious
#' correlations because of the alignment of relatively tiny rings in a reference and the series
#' you are evaluating.
#'
#' P2YrsL is similar in concept to computing AR residuals (a form of "prewhitening") - the default
#' in dplR. P2YrsL is the default in CDendro and CooRecorder. The essential element of both
#' approaches is that we remove the mid- to low-frequency variation and highlight the high-frequency
#' for cross-dating purposes only. The high-frequency variation facilitates checking both
#' statistical (correlations) and visual (squiggly line plots) correspondence of two series.
#'
#' @param rwl A rwl object (read in by dplR's \code{\link[dplR]{read.rwl}}). Essentially a
#' data.frame with columns names as series IDs and years as rownames.
#' @param limit A numeric vector specifying the number of standard deviations beyond the mean to
#' limit (cut off) the resulting series. Default is 2.6, the same value used in
#' CDendro and CooRecorder.
#'
#' @details
#' See https://cdendro.se/wiki/index.php/Normalization
#'
#' @return A rwl-style data.frame with the transformed ring width series.
#'
#' @references
#' Larsson & Larsson (2023) \emph{CDendro and CooRecorder programs of the CDendro package},
#'  Cybis Elektronik & Data AB. https://www.cybis.se/forfun/dendro/index.htm
#'
#' @import dplR
#'
#' @export
#'
#' @examples
#' library(dplR)
#' data(ca533)
#' # before
#' ca533 |> spag.plot()
#' # after
#' ca533 |> p2yrsL() |> spag.plot()
#'

p2yrsL <- function(rwl, limit = 2.6) {

  # Error catching
  stopifnot("rwl is not an object of class 'rwl', 'data.frame', or 'matrix'" =
              data.class(rwl) %in% "rwl" |
              data.class(rwl) %in% "data.frame" |
              data.class(rwl) %in% "matrix"
  )

  stopifnot("rwl has no rownames (must be years only) or no colnames (must be series IDs only)" =
              !is.null(rownames(rwl)) |
              !is.null(colnames(rwl))
  )

  if (apply(rwl, MARGIN = 2, FUN = \(x) all(is.na(x))) |> any() == TRUE) {
    these_are_NA <- colnames(rwl)[which(apply(rwl, MARGIN = 2, FUN = \(x) all(is.na(x))) == TRUE)]
    stop("The following series have no values (all NAs): " , paste(these_are_NA, collapse = ", "))
  }

  stopifnot("limit must be a positive numeric between 0:10" =
              data.class(limit) %in% "numeric" &
              limit > 0 & limit < 10
  )

  # Calculate rolling (2 year) sums
  # Add an NA at the beginning for proper alignment with the actual ring widths,
  # such that p2yrs = rw(i) / (rw(i) + rw(i-1))
  sum2yrs <- apply(rwl, MARGIN = 2, \(x) {
    c(NA, base::rowSums(stats::embed(x, dimension = 2)))
    }) |>
    as.data.frame()

  # sum2yrs <- rbind(rep(NA, ncol(sum2yrs)), sum2yrs)
  # Calculate the proportion of 2 years for each year
  p2yrs <- rwl / sum2yrs

  ## Set up outer margins for "limited" trimming
  # Calculate means
  means <- colMeans(p2yrs, na.rm = TRUE)
  # colMeans(na.rm = T) is equivalent to, but faster than:
  # apply(p2yrs, MARGIN = 2, FUN = \(x) mean(x, na.rm = TRUE))

  # Limits are defined in terms of a threshold * the SD from the means
  limits <- limit * apply(p2yrs, MARGIN = 2, FUN = \(x) sd(x, na.rm = TRUE))
  hi_lim <- means + limits # Upper limits
  lo_lim <- means - limits # Lower limits

  # Cut off the high extremes
  p2yrsL1 <- mapply(FUN = \(x, y) {
    ifelse(x > y, y, x)
  }, x = p2yrs, y = hi_lim) |> as.data.frame()

  # Cut off the low extremes
  p2yrsL <- mapply(FUN = \(x, y) {
    ifelse(x < y, y, x)
  }, x = p2yrsL1, y = lo_lim) |> as.data.frame()

  # Transfer the rownames (years)
  rownames(p2yrsL) <- rownames(rwl)
  # Return the final rwl-data.frame
  p2yrsL
}
