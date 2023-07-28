#' Curve intervention detection - statistical identification and removal of outliers in tree ring series
#'
#' @description
#' This function implements an analysis based on the curve intervention detection techniques described by Druckenbrod et al,
#' Rydval et al, and Rydval et al. The
#'
#'
#' @param rwl A rwl object (read in by dplR's `read.rwl()`). Essentially a data.frame with columns names as series IDs and years as rownames.
#' @param detrend.method A character vector specifying the detrending method to use. Passes to `dplR::detrend()`. Default is "AgeDepSpline".
#' @param nyrs A numeric vector that determines the flexibility of the "AgeDepSpline" or the "Spline" detrending methods. Passes to `dplR::detrend()`. The default is 50 years.
#' @param min.win The minimum outlier length in years (i.e., a moving window) to search for. The default is 9.
#' @param max.win The maximum outlier length in years (i.e., a moving window) to search for. The default is 30.
#' @param thresh The outlier detection threshold, corresponding to the number of devations from the mean. The default is 3.29, following Druckenbrod et al. 2013.
#' @param out.span Parameter to determine the wiggliness of the loess splines fit to outlier periods. Higher numbers = more wiggles. Passes to `base::loess()`. The default is 1.
#' @param max.iter The maximum number of iterations to run the outlier detection and removal processes. The default is 15.
#'
#' @details
#'
#'
#' @return A 3-element list containing the "disturbance-free" series, the disturbance index, and a record of the outlier
#' detection and removal iterations.
#'
#' @references
#' Larsson & Larsson (2023) \emph{CDendro and CooRecorder programs of the CDendro package},
#'  Cybis Elektronik & Data AB. https://www.cybis.se/forfun/dendro/index.htm
#'
#' @import dplR
#' @import zoo
#' @import stats
#' @export
#'
#' @examples
#' library(dplR)
#' data("co021")
#' # before
#' ci_detect(co021, max.iter = 5)


ci_detect <- function(rwl,
                      detrend.method = "AgeDepSpline",
                      nyrs = NULL,
                      min.win = 9,
                      max.win = 30,
                      thresh = 3.29,
                      out.span = 1,
                      max.iter = 10) {

  ## Run cp_detrend to power transform and detrend the rwl
  cp_out <- cp_detrend(rwl, detrend.method = detrend.method, nyrs = nyrs)
  # 1st element of cp_out is a rwl-data.frame of the residual transformed and detrended series-
  # take this and turn it into a list, with NAs removed from each series
  # Simplify = FALSE keeps the rownames (which are the years)
  cp_list <- apply(cp_out[[1]], MARGIN = 2, simplify = FALSE, FUN = \(x) na.omit(x))

  # Run the out_det_rem() process
  # The for loop makes sense here - after the 1st iteration (based on the original data),
  # each subsequent iteration uses the values generated in the previous iteration.
  # Operations in the out_det_rem() function are vectorized, so this is reasonably efficient.
  out_iter <- vector("list", length = max.iter + 1)
  # The 1st element in the outlier iteration list is the initial data. This has to be a 2-element list
  # to match the output of out_det_rem(). It is just an empty filler.
  start.list <- list(cp_list, vector("list", length = length(cp_list)))
  names(start.list) <- c("Original RWI","Empty filler")
  out_iter[[1]] <- start.list
  names(out_iter) <- 0:max.iter
  for (i in 2:(max.iter+1)) {
    out_iter[[i]] <- out_det_rem(out_iter[[i-1]][[1]], min.win = min.win, max.win = max.win, span = out.span)
  }

  ## Take the last element of out_iter as the final output series...
  # add back the detrend curve
  retrended <- mapply(FUN = \(x, y) {
    x + na.omit(y)
  }, x = out_iter[[length(out_iter)]][[1]], y = cp_out[["Detrend curves"]])

  # & undo any transformation - this results in a "disturbance-free" series in original units
  untransformed <- mapply(FUN = \(x, y) {
    # Undo the transformations
    if (x["action"] %in% "log10 transformed") {
      10^y
    } else {
      if (x["action"] %in% "Power transformed") {
        y^(1/as.numeric(x[["optimal.pwr"]]))
      } else {
        y
      }
    }
  }, x = cp_out[["Metadata about transformations"]], y = retrended)

  # Calculate the disturbance index -
  # this is the difference between the disturbance-free series & the original series
  dis_index <- mapply(FUN = \(x, y) {
    na.omit(x) - y
  }, x = rwl, y = untransformed, SIMPLIFY = FALSE)

  # Restore the rwl-data.frame format for the disturbance-free and disturbance index series
  # (with NAs for the years that aren't covered)
  # 1st add a year column to each element of the lists
  series_names <- names(untransformed)
  untransformed_rwl0 <- mapply(FUN = \(x, y) {
    x <- as.data.frame(x)
    colnames(x) <- y
    x[,"year"] <- rownames(x) |> as.numeric()
    x
  }, x = untransformed, y = series_names, SIMPLIFY = FALSE)

  dis_index_rwl0 <- mapply(FUN = \(x, y) {
    x <- as.data.frame(x)
    colnames(x)[1] <- y
    x[,"year"] <- rownames(x) |> as.numeric()
    x
  }, x = dis_index, y = series_names, SIMPLIFY = FALSE)

  # Can apply merge() iteratively with Reduce():
  untransformed_rwl <- Reduce(f = \(x, y) merge(x, y, by = "year", all = TRUE), untransformed_rwl0)

  dis_index_rwl <- Reduce(f = \(x, y) merge(x, y, by = "year", all = TRUE), dis_index_rwl0)

  # Make years as rownames
  rownames(untransformed_rwl) <- untransformed_rwl$year
  untransformed_rwl <- untransformed_rwl[,!c(colnames(untransformed_rwl) %in% "year")]

  rownames(dis_index_rwl) <- dis_index_rwl$year
  dis_index_rwl <- dis_index_rwl[,!c(colnames(dis_index_rwl) %in% "year")]

  ## Last steps are to output the main results (disturbance-free, disturbance index, original series)
  # and all of the iterations of the outlier removal process, with the trend curves, etc.
  # all of these can then be plotted in the plot_ci_detect() function
  ci_output_list <- list(untransformed_rwl, dis_index_rwl, out_iter, rwl)
  names(ci_output_list) <- c("Disturbance-free series", "Disturbance index", "Outlier removal iterations", "Original series")
  ci_output_list

} # End of the ci_detect() function
