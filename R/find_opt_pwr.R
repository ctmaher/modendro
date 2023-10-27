#' Find the optimal power for transforming a tree ring series to approximate homoscedasticity
#'
#' @description
#' Use a spread versus level method to estimate the optimal power of transformation for a tree ring series
#'
#'
#' @details
#' This function estimates the optimal power of transformation via the local spread versus level relationship of each series,
#' where local spread is defined as the absolute value of the first
#' differences, S, (|rwt - rwt-1|) and the local level is the arithmetic mean of each pair of adjacent values, M, (rwt + rwt-1)/2.
#' The spread versus level relationship is then modeled in a simple linear regression as log10(S) ~ log10(M). The optimal power of transformation
#' is then estimated as p = |1-slope| of this regression model.
#' See ?cp_detrend and Cook and Peters (1997) for more details.
#'
#' Note that this is a very simple function that only implements the spread vs. level estimation
#' of an optimal power for making a tree ring series ± homoscedastic. This function does not select which of these estimated
#' powers to use, that is performed in the  \code{\link{pwr_t_rwl}} function.
#'
#' @param rwl A rwl object (read in by dplR's  \code{\link[dplR]{read.rwl}}). Essentially a data.frame with columns names as series IDs and years as rownames.
#'
#' @return A named numeric vector with the series IDs (colnames) and the estimated optimal power of transformation
#'
#' @references
#' Cook, E. R., and Peters, K. (1997) Calculating unbiased tree-ring indices for the study of climatic and environmental change.
#' \emph{The Holocene}, \strong{7}(3), 361-370.
#'
#' @seealso \code{\link{pwr_t_rwl}}, \code{\link{cp_detrend}}
#'
#' @export
#'
#' @examples
#' library(dplR)
#' data("ca533")
#' find_opt_pwr(rwl = ca533)


find_opt_pwr <- function(rwl) {
  # Error catching
  stopifnot(
    "rwl is not an object of class 'rwl', 'data.frame', or 'matrix'" =
      data.class(rwl) %in% "rwl" |
      data.class(rwl) %in% "data.frame" |
      data.class(rwl) %in% "matrix"
  )

  stopifnot(
    "rwl has no rownames (must be years only) or no colnames (must be series IDs only)" =
      !is.null(rownames(rwl)) |
      !is.null(colnames(rwl))
  )

  if (apply(rwl, MARGIN = 2, FUN = \(x) all(is.na(x))) |> any() == TRUE) {
    these_are_NA <-
      colnames(rwl)[which(apply(rwl, MARGIN = 2, FUN = \(x) all(is.na(x))) == TRUE)]
    stop("The following series have no values (all NAs): " ,
         paste(these_are_NA, collapse = ", "))
  }

  orig.IDs <- colnames(rwl) # original series names in original order

  # absolute value of 1st differences.
  # set up lagged vector for 1st differences
  rwl.lag <-
    rbind(rep(NA, ncol(rwl)), rwl[1:(nrow(rwl) - 1), , drop = FALSE])

  rwl.lag <- rwl.lag[, orig.IDs]

  # absolute value of 1st differences
  # Cook & Peters equate this as the standard deviation (it is not the same, but 1st differences is what they use).
  diffs <- abs(rwl - rwl.lag)
  lmean <- (rwl + rwl.lag) / 2

  diffs <- diffs[, orig.IDs]
  lmean <- lmean[, orig.IDs]

  # Any zero values in diffs or means should be avoided
  diffs.ind <- which(diffs < 0.001, arr.ind = TRUE)
  lmean.ind <- which(lmean < 0.001, arr.ind = TRUE)

  diffs[diffs.ind] <- NA
  diffs[lmean.ind] <- NA
  lmean[diffs.ind] <- NA
  lmean[lmean.ind] <- NA

  # Make the dfs into lists - use apply(MARGIN = 2, simplify = F) & an identity function
  diffs <- apply(diffs,
                 MARGIN = 2,
                 FUN = \(x) x,
                 simplify = FALSE)
  lmean <- apply(lmean,
                 MARGIN = 2,
                 FUN = \(x) x,
                 simplify = FALSE)

  # Get the slopes
  slopes <- mapply(FUN = \(x, y) {
    lm(log10(y) ~ log10(x), na.action = "na.exclude")$coef[2][[1]]
  }, x = lmean[orig.IDs], y = diffs[orig.IDs])

  # Make sure the order matches the original rwl
  slopes <- slopes[orig.IDs]

  # Determine optimal power of transformation by subtracting slope from 1
  abs(1 - slopes)
} # End of function
