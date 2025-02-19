#' Replace NA values when they occur within a rwl collection of tree-ring series
#'
#' @description
#' Simple utility function to replace NA values when they occur within a rwl. Chains of NA values
#' are common at the beginnings and ends of tree-ring series within rwl files, and by some
#' conventions missing or locally absent rings are marked with NA. This can cause problems for many
#' analyses.
#'
#' This is a wrapper function for the simpler \code{\link{replace_internal_NAs}} that works on any
#' simple numeric vector. One important difference is that \code{\link{rwl_replace_internal_NAs}}
#' will use the rownames (which represent years) to order the rwl before finding and replacing
#' internal NAs.
#'
#' Warning: for the default of new.val = 0, make sure that the NAs do indeed represent missing or
#' locally absent rings!
#'
#' Also, for some analyses, 0 values may cause issues. In these cases it might make sense to
#' provide some nominally small value, e.g., 0.001, 0.0001.
#'
#' @param rwl A rwl-type data.frame (e.g., read in by \code{\link[dplR]{read.rwl}}). Essentially a
#' data.frame with columns names as series IDs and years as rownames.
#' @param new.val a numeric vector of length 1. The default is zero
#'
#' @export
#'
#' @examples
#' data(PerkinsSwetnam96)
#' # replace some values with NAs
#' PS.NAs <- PerkinsSwetnam96
#' # Give at least one series an NA value
#' PS.NAs[sample(1:nrow(PS.NAs), 1),] <- NA
#'
#' # Try to run a common dplR function that will throw an error if there are missing values in the
#' # middle of series
#' tryCatch(
#' dplR::detrend(PS.NAs, method = "Spline"),
#' error = \(e) e$message
#' )
#'
#' PS.noNAs <- rwl_replace_internal_NAs(rwl = PS.NAs, new.val = 0)
#'
#' # No more error
#' dplR::detrend(PS.noNAs, method = "Spline")

rwl_replace_internal_NAs <- function(rwl = NULL,
                                     new.val = NULL) {

  ###### rwl
  stopifnot(
    "Arg rwl is not an object of class 'rwl' or 'data.frame'" =
      data.class(rwl) %in% "rwl" |
      data.class(rwl) %in% "data.frame"
  )

  stopifnot("new.val must be a numeric vector of length == 1" =
              is.numeric(new.val) &
              length(new.val) == 1
  )

  # Sort the RWL by year first
  rwl <- rwl[order(as.numeric(rownames(rwl)), decreasing = FALSE),]

  # Use apply to run the simple replace_internal_NAs on each series
  rwl.out <- apply(rwl, MARGIN = 2, FUN = \(series) {
    replace_internal_NAs(series, new.val = new.val)
  }, simplify = FALSE) |>
    do.call(what = "cbind") |>
    as.data.frame()

  rownames(rwl.out) <- rownames(rwl)
  class(rwl.out) <- c("rwl","data.frame")

  return(rwl.out)

}


