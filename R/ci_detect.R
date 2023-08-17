#' Curve intervention detection - statistical identification and removal of disturbances in tree ring series
#'
#' @description
#' This function implements an analysis based on the curve intervention detection techniques described by Druckenbrod et al. 2013,
#' Rydval et al. 2016, and Rydval et al. 2018. The modendro implementation differs from the version described by these authors. See Details for more information.
#' The basic motivations for performing this analysis are to identify, quantify, & remove abrupt changes in radial growth. You might use
#' this information to help make inferences about forest stand dynamics, insect outbreaks, or isolate a stronger climate signal, depending on your dendroecological or dendroclimatological tendencies.
#'
#'
#' @param rwl A rwl object (read in by \code{\link[dplR]{read.rwl}}). Essentially a data.frame with columns names as series IDs and years as rownames.
#' @param detrend.method A character vector specifying the detrending method to use. Passes to \code{\link[dplR]{detrend}}. Default is "Mean".
#' @param nyrs A numeric vector that determines the flexibility of the "AgeDepSpline" or the "Spline" detrending methods. Passes to \code{\link[dplR]{detrend}}. The default is 50 years or 1/3 the series length.
#' @param min.win The minimum disturbance length in years (i.e., a moving window) to search for. The default is 9.
#' @param max.win The maximum disturbance length in years (i.e., a moving window) to search for. The default is the smallest of 30 years or 1/3 of the series length.
#' @param thresh The disturbance detection threshold, corresponding to the number of deviations from the robust mean. The default is 3.29, following Druckenbrod et al. 2013.
#' @param out.span Parameter to determine the wiggliness of the loess splines fit to disturbances periods. Higher numbers = less wiggles. Passes to \code{\link[stats]{loess}}. The default is 1.25
#' @param max.iter The maximum number of iterations to run the disturbance detection and removal processes. The default is 10.
#'
#' @details
#' Intervention detection is a statistical time series approach to identify and remove abrupt changes in radial growth. The "curve" part
#' of curve intervention detection (CID) describes the type of line fitted to each period that is identified as a disturbance.
#' Following Rydval et al. (2016, 2018), the implementation here uses a version of the equation presented by Warren & MacWilliam (1981), with
#' slightly different coefficients than the Hugershoff curve used by Rydval et al. Namely, the following values are fixed: b = 1, and d = 0.
#' The Warren & MacWilliam curve has a 't' coefficient which we use here. These changes improve the robustness of the fits. Fitting is done
#' via \code{\link[stats]{nls}}. A major difference for the `modendro` implementation of CID is that in the cases where the \code{\link[stats]{nls}} fits fail, a \code{\link[stats]{loess}} spline is fit
#' instead. See \code{\link{out_det_rem}} for more details.
#'
#'
#' \code{\link{function}}
#'
#' @return A 5-element list containing the "disturbance-free" series, the disturbance index,
#' a data.frame containing basic data on all the detected disturbances,
#' a record of the outlier detection and removal iterations, and the output from the cp_detrend() process.
#'
#' @references
#' Druckenbrod, D. L., N. Pederson, J. Rentch, and E. R. Cook. 2013. A comparison of times series approaches for dendroecological reconstructions of past canopy disturbance events. \emph{Forest Ecology and Management}, \strong{302}, 23–33.
#'
#' Rydval, M., D. Druckenbrod, K. J. Anchukaitis, and R. Wilson. 2016. Detection and removal of disturbance trends in tree-ring series for dendroclimatology. \emph{Canadian Journal of Forest Research}, \strong{401}, 387–401.
#'
#' Rydval, M., D. L. Druckenbrod, M. Svoboda, V. Trotsiuk, P. Janda, M. Mikoláš, V. Čada, R. Bače, M. Teodosiu, and R. Wilson. 2018. Influence of sampling and disturbance history on climatic sensitivity of temperature limited conifers. \emph{The Holocene}, \strong{28}(10), 1574-1587.
#'
#' Warren, W. G., and S. L. MacWilliam. 1981. Test of a new method for removing the growth trend from dendrochronological data. \emph{Tree Ring Bulletin}, \strong{41}, 55–66.
#'
#'
#' @import dplR
#' @importFrom zoo rollmean
#' @import stats
#' @import DescTools
#' @export
#'
#' @examples
#' library(dplR)
#' data("co021")
#' # before
#' ci_detect(co021, max.iter = 5)


ci_detect <- function(rwl,
                      detrend.method = "None",
                      nyrs = NULL,
                      min.win = 9,
                      max.win = 30,
                      var.type = "s_bi",
                      thresh = 3.29,
                      out.span = 1.25,
                      max.iter = 10,
                      add.recent.rwi = TRUE
) {

  ## Error catching & warnings
  #
  stopifnot("rwl is not an object of class 'rwl', 'data.frame', or 'matrix'" =
              data.class(rwl) %in% "rwl" |
              data.class(rwl) %in% "data.frame" |
              data.class(rwl) %in% "matrix"
  )

  #
  stopifnot("rwl has no rownames (must be years only) or no colnames (must be series IDs only)" =
              !is.null(rownames(rwl)) |
              !is.null(colnames(rwl))
  )

  #
  if (apply(rwl, MARGIN = 2, FUN = \(x) all(is.na(x))) |> any() == TRUE) {
    these_are_NA <- colnames(rwl)[which(apply(rwl, MARGIN = 2, FUN = \(x) all(is.na(x))) == TRUE)]
    stop("The following series have no values (all NAs): " , paste(these_are_NA, collapse = ", "))
  }

  #
  stopifnot("min.win is too small for effective curve fitting. Choose a value of 5 years or larger." =
              min.win >= 5
  )


  ## Run cp_detrend to power transform and detrend the rwl
  cp_out <- cp_detrend(rwl, detrend.method = detrend.method, nyrs = nyrs, pos.slope = TRUE)
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
    out_iter[[i]] <- out_det_rem(out_iter[[i-1]][[1]],
                                 min.win = min.win,
                                 max.win = max.win,
                                 thresh = thresh,
                                 var.type = var.type,
                                 out.span = out.span,
                                 add.recent.rwi = add.recent.rwi
    )
  }
  # Remove the 0 iteration
  out_iter <- out_iter[2:(max.iter+1)]

  ## Take the last element of out_iter as the final output series...
  # add back the detrend curve (or the mean)

  retrended <- mapply(FUN = \(x, y) {
      x + na.omit(y)
    }, x = out_iter[[length(out_iter)]][["Corrected RWI"]], y = cp_out[["Detrend curves"]])

  # Undo any transformation - this results in a "disturbance-free" series in original units
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
  }, x = cp_out[["Transformation metadata"]], y = retrended)

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

  # Make years as rownames & then remove the year column
  rownames(untransformed_rwl) <- untransformed_rwl$year
  untransformed_rwl <- untransformed_rwl[,!c(colnames(untransformed_rwl) %in% "year")]

  rownames(dis_index_rwl) <- dis_index_rwl$year
  dis_index_rwl <- dis_index_rwl[,!c(colnames(dis_index_rwl) %in% "year")]

  # Collapse the detection & removal iterations into a data.frame that contains all the detected disturbances
  dist_det <- lapply(out_iter, FUN = \(x) {
    x1 <- x[["Outlier curves"]]

    x2 <- lapply(x1, FUN = \(y) {
      if (!is.character(y)){
        y[1, c("year", "dir", "dur", "eq")]
      }
    }) |> do.call(what = "rbind")

    x2$series <- rownames(x2)
    x2[, c("series", "year", "dir", "dur", "eq")]

  }) |> do.call(what = "rbind")

  ## Last steps are to output the main results (disturbance-free, disturbance index, original series)
  # and all of the iterations of the outlier removal process, with the trend curves, etc.
  # all of these can then be plotted in the plot_ci_detect() function
  ci_output_list <- list(untransformed_rwl, dis_index_rwl, dist_det, out_iter, cp_out)
  names(ci_output_list) <- c("Disturbance-free series", "Disturbance index",
                             "Detected disturbances", "Disturbance removal iterations",
                             "Cook & Peters detrend")
  ci_output_list

} # End of the ci_detect() function
