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

p2yrsL <- function(rwl = NULL, limit = 2.6, out.format = "long") {
  . <- j_variable <- NULL
  #### data format detection and conversion ####
  # Also does some quality checking. Leverages rwl_longer() and longer_rwl()
  rwl.long <- rwl_dt_conv(data = rwl,
                          series.name = "series",
                          dat.name = "rw",
                          new.val.internal.na = 0)

  ## More error catching

  # Check for negative values (not allowed)
  stopifnot(
    "rwl has negative values. p2yrsL is intended for raw ring width values, which must be > 0." =
      #all(apply(rwl, MARGIN = 2, FUN = \(x) all(as.numeric(na.omit(x)) >= 0)))
      all(rwl.long[, .(neg.val = any(rw >= 0)), by = series]$neg.val)
  )

  # Also cannot allow zero values - if there are several in a row we divide by 0
  # Replace ≤0 values in the rwl with the minimum possible non-zero value given the resolution of
  # the data
  # & make sure they match
  min.value <- ifelse(rwl.long[, .(n.char = nchar(sub(".", "", rw, fixed=TRUE))), by = series]$n.char |>
                        max() <= 3,
                      0.01, 0.001)

  # Replace 0s with min.value
  rwl.long[rw < min.value, rw := min.value,]

  # limit arg
  stopifnot("limit must be a positive numeric between 0:10" =
              data.class(limit) %in% "numeric" &
              limit > 0 & limit < 10
  )

  # out.format arg
  stopifnot(
    "out.format options are c('rwl', 'long')" =
      is.character(out.format) &&
      length(out.format) == 1 &&
      out.format %in% c("rwl", "long")
  )

  # Calculate rolling (2 year) sums
  # Add an NA at the beginning for proper alignment with the actual ring widths,
  # such that p2yrs = rw(i) / (rw(i) + rw(i-1))
  data.table::setorder(rwl.long, series, year)

  rwl.long[, sum2yrs := c(NA, base::rowSums(stats::embed(rw, dimension = 2))), by = series]

  # Calculate the proportion of 2 years for each year by dividing.
  rwl.long[, p2yrs := rw / sum2yrs]

  ## Set up outer margins for "limited" trimming
  rwl.long[, `:=`(hi_lim = mean(p2yrs, na.rm = TRUE) + limit * sd(p2yrs, na.rm = TRUE),
                  lo_lim = mean(p2yrs, na.rm = TRUE) - limit * sd(p2yrs, na.rm = TRUE)),
           by = series]

  # Replace the values that go outside of the limits with the limit values
  rwl.long[, p2yrsL := data.table::fifelse(p2yrs > hi_lim, hi_lim, p2yrs)]
  rwl.long[, p2yrsL := data.table::fifelse(p2yrs < lo_lim, lo_lim, p2yrs)]

  # Return the final product
  if (out.format %in% "rwl") {
    longer_rwl(rwl.long, series.name = "series", dat.name = "p2yrsL")
  } else { # return the long format data.table
    rwl.long
  }
} ## End of function
