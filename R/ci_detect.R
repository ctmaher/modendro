#' Curve intervention detection - statistical identification and removal of disturbances in tree
#' ring series
#'
#' @description
#' This function implements an analysis based on the curve intervention detection techniques
#' described by Druckenbrod et al. 2013, Rydval et al. 2016, and Rydval et al. 2018. The modendro
#' implementation differs from the version described by these authors. See Details for more
#' information.
#' The basic motivations for performing this analysis are to identify, quantify, & remove abrupt
#' changes in radial growth. You might use this information to help make inferences about forest
#' stand dynamics, insect outbreaks, or isolate a stronger climate signal, depending on your
#' dendroecological or dendroclimatological tendencies.
#' The main output are "disturbance-free" series that are in original ring width units (e.g., mm).
#'
#'
#' @param rwl A rwl-type data.frame (e.g., read in by \code{\link[dplR]{read.rwl}}). Essentially a
#' data.frame with columns names as series IDs and years as rownames.
#' @param detrend.method A character vector specifying the detrending method to use. Passes to
#' \code{\link[dplR]{detrend}} via \code{\link{cp_detrend}}. Default is "Mean".
#' @param nyrs A numeric vector that determines the flexibility of the `"AgeDepSpline"` or
#' the `"Spline"` detrending methods. Passes to \code{\link[dplR]{detrend}}. The default is 50
#' years or 1/3 the series length.
#' @param min.win The minimum disturbance length in years (i.e., a moving window) to search for.
#' The default is 9.
#' @param max.win The maximum disturbance length in years (i.e., a moving window) to search for.
#' The default is the smallest of 30 years or 1/3 of the series length.
#' @param thresh The disturbance detection threshold, corresponding to the number of deviations
#' from the robust mean. The default is 3.29, following Druckenbrod et al. 2013.
#' @param dist.span Parameter to determine the wiggliness of the loess splines fit to disturbances
#'  periods. Higher numbers = less wiggles. Passes to \code{\link[stats]{loess}}.
#'  The default is 1.25
#' @param max.iter The maximum number of iterations to run the disturbance detection and removal
#' processes. The default is 10.
#'
#'
#' @details
#' Intervention detection is a statistical time series approach to identify and remove abrupt
#' changes in radial growth. Disturbances are detected using moving averages of autoregressive
#' residuals of transformed detrended ring width series. The "curve" part of curve intervention
#' detection (CID) describes the type of line fitted to each period that is identified as a
#' disturbance. The implementation here uses a version of the equation presented by Warren &
#' MacWilliam (1981), with slightly different coefficients than the Hugershoff curve used by
#' Rydval et al. Namely, the following values are fixed: b = 1, and d = 0. The Warren & MacWilliam
#' curve has a 't' coefficient which we use here. These changes improve the robustness of the fits.
#' Fitting is done via \code{\link[stats]{nls}}. A major difference for the `modendro`
#' implementation of CID is that in the cases where the \code{\link[stats]{nls}} fits fail,
#' a \code{\link[stats]{loess}} spline is fit instead. See \code{\link{dist_det_rem}} for more
#' details.
#'
#' The choice of initial detrending can have an effect on the subsequent disturbance detection and
#'  removal, so choose wisely. All detrending is done with inside the function with
#'  \code{\link{cp_detrend}} so that the resulting series are transformed to homogenize variance
#'  and are residuals with a mean of 0.
#'
#' The final steps of the process involve adding back the initial detrending curves and reversing
#' the transformations so that the resulting "disturbance-free" series are in the original ring
#' width units (typically mm). These series can then be treated as any other ring width series,
#' with the type of detrending/standardization depending on the research goals.
#'
#' @return A named list with the following elements:
#' 1) "Disturbance-free series" - an rwl-data.frame of the final 'disturbance-free' series in
#' original units (ring width; typically mm).
#' 2) "Disturbance index" - an rwl-data.frame of the differences (in original units) between the
#' original ring widths and the 'disturbance-free' series.
#' 3) "Detected disturbances" - a list of data.frames containing the AR residuals, a moving average
#' series, series robust means, and detection thresholds for the max disturbance detected.
#' Mostly for plotting disturbance iterations in \code{\link{plot_ci_detect}}.
#' 4) "Disturbance removal iterations" - a list of lists containing the details on each disturbance
#' detected and removed for each iteration and each series.
#' 5) "Cook & Peters detrend" - the output from the \code{\link{cp_detrend}} function.
#'
#'
#' @references
#' Druckenbrod, D. L., N. Pederson, J. Rentch, and E. R. Cook. (2013) A comparison of times series
#' approaches for dendroecological reconstructions of past canopy disturbance events.
#' \emph{Forest Ecology and Management}, \strong{302}, 23–33.
#'
#' Rydval, M., D. Druckenbrod, K. J. Anchukaitis, and R. Wilson. (2016) Detection and removal of
#' disturbance trends in tree-ring series for dendroclimatology.
#' \emph{Canadian Journal of Forest Research}, \strong{401}, 387–401.
#'
#' Rydval, M., D. L. Druckenbrod, M. Svoboda, V. Trotsiuk, P. Janda, M. Mikoláš, V. Čada, R. Bače,
#' M. Teodosiu, and R. Wilson. (2018) Influence of sampling and disturbance history on climatic
#' sensitivity of temperature limited conifers.
#' \emph{The Holocene}, \strong{28}(10), 1574-1587.
#'
#' Warren, W. G., and S. L. MacWilliam. 1981. Test of a new method for removing the growth trend
#' from dendrochronological data.
#' \emph{Tree Ring Bulletin} \strong{41}, 55–66.
#'
#' @seealso \code{\link{dist_det_rem}}, \code{\link{plot_ci_detect}}, \code{\link{cp_detrend}},
#'  \code{\link{plot_cp_detrend}}
#'
#' @import dplR
#' @import stats
#' @export
#'
#' @examples
#' library(dplR)
#' data("co021")
#'
#'
#' co021.ci.out <- ci_detect(co021,
#' detrend.method = "AgeDepSpline",
#' max.iter = 5)
#'
#' names(co021.ci.out)
#'
#' # dplR::spag.plot() can give you a quick look at some of the key outputs
#' spag.plot(co021.ci.out[["Disturbance-free series"]])
#' spag.plot(co021.ci.out[["Disturbance index"]])
#'
#' # You can get more detailed plots from plot_ci_detect
#' co021.ci.out.plots <- plot_ci_detect(co021.ci.out)
#' co021.ci.out.plots[["Disturbance detection & removal plots"]][["641121"]]
#' # outputs 5 plots (more if max.iter is increased)
#' co021.ci.out.plots[["Final disturbance-free series plots"]][["641121"]]
#'
#' # The power transformation and detrending steps can also be retrieved using plot_cp_detrend:
#' co021.detrend.plots <- plot_cp_detrend(co021.ci.out[["Cook & Peters detrend"]])
#' co021.detrend.plots[["641121"]]
#'


ci_detect <- function(rwl = NULL,
                      detrend.method = "Mean",
                      nyrs = NULL,
                      min.win = 9,
                      max.win = 30,
                      thresh = 3.29,
                      dist.span = 1.25,
                      max.iter = 10
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
  stopifnot("min.win is too small for effective curve fitting.
            Choose a value of 5 years or larger."
            =
              min.win >= 5
  )

  stopifnot("thresh must be > 0" =
              thresh > 0
  )

  # save the original order of the names to realign at each step
  orig.IDs <- colnames(rwl)
  ## Run cp_detrend to power transform and detrend the rwl
  cp_out <- cp_detrend(rwl,
                       detrend.method = detrend.method,
                       nyrs = nyrs,
                       pos.slope = FALSE,
                       standardize = FALSE)
  # 1st element of cp_out is a rwl-data.frame of the residual transformed and detrended series-
  # take this and turn it into a list, with NAs removed from each series
  # Simplify = FALSE keeps the rownames (which are the years)
  cp_list <- apply(cp_out[["Resid. detrended series"]],
                   MARGIN = 2,
                   simplify = FALSE,
                   FUN = \(x) {
                     x1 <- as.numeric(na.omit(x))
                     names(x1) <- names(x[!is.na(x)])
                     x1
                   })
  cp_list <- cp_list[orig.IDs] # make sure names are ordered
  # Run the dist_det_rem() process
  # The for loop makes sense here - after the 1st iteration (based on the original data),
  # each subsequent iteration uses the values generated in the previous iteration.
  # Operations in the dist_det_rem() function are vectorized, so this is reasonably efficient.
  dist_iter <- vector("list", length = max.iter + 1)
  # The 1st element in the disturbance iteration list is the initial data.
  # This has to be a 2-element list
  # to match the output of dist_det_rem(). It is just an empty filler.
  start.list <- list(cp_list, vector("list", length = length(cp_list)))
  names(start.list) <- c("Original RWI","Empty filler")
  dist_iter[[1]] <- start.list
  names(dist_iter) <- 0:max.iter
  for (i in 2:(max.iter+1)) {
    dist_iter[[i]] <- dist_det_rem(dist_iter[[i-1]][[1]],
                                   min.win = min.win,
                                   max.win = max.win,
                                   thresh = thresh,
                                   dist.span = dist.span
    )
  }
  # Remove the 0 iteration
  dist_iter <- dist_iter[2:(max.iter+1)]

  ## Take the last element of dist_iter as the final output series...
  # add back the detrend curve (or the mean)

  # Make the detrend curves a list
  detrend.curves <- apply(cp_out[["Detrend curves"]][orig.IDs],
                          MARGIN = 2,
                          FUN = \(x) {as.numeric(na.omit(x))},
                          simplify = FALSE)

  retrended <- mapply(FUN = \(x, y) {
    x + as.numeric(na.omit(y))
  }, x = dist_iter[[length(dist_iter)]][["Corrected RWI"]][orig.IDs],
  y = detrend.curves[orig.IDs])

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
  }, x = cp_out[["Transformation metadata"]][orig.IDs], # Make sure they share the same order!
  y = retrended[orig.IDs])

  # Calculate the disturbance index -
  # this is the difference between the disturbance-free series & the original series
  # Note that there is an assumption that series from na.omit(x) and y are in the same order of
  # years

  # The original rwl as a list, with no NAs
  orig.rwl <- apply(rwl,
                    MARGIN = 2,
                    FUN = \(x) {as.numeric(na.omit(x))},
                    simplify = FALSE)

  dis_index <- mapply(FUN = \(x, y) {
    as.numeric(na.omit(x)) - y
  }, x = orig.rwl[orig.IDs], y = untransformed[orig.IDs], SIMPLIFY = FALSE)

  # Restore the rwl-data.frame format for the disturbance-free and disturbance index series
  # (with NAs for the years that aren't covered)
  # 1st add a year column to each element of the lists
  untransformed_rwl0 <- mapply(FUN = \(x, y) {
    x <- as.data.frame(x)
    colnames(x) <- y
    x[,"year"] <- rownames(x) |> as.numeric()
    x
  }, x = untransformed[orig.IDs], y = orig.IDs, SIMPLIFY = FALSE)

  dis_index_rwl0 <- mapply(FUN = \(x, y) {
    x <- as.data.frame(x)
    colnames(x)[1] <- y
    x[,"year"] <- rownames(x) |> as.numeric()
    x
  }, x = dis_index[orig.IDs], y = orig.IDs, SIMPLIFY = FALSE)

  # Can apply merge() iteratively with Reduce():
  untransformed_rwl <- Reduce(f = \(x, y) merge(x, y, by = "year", all = TRUE),
                              untransformed_rwl0[orig.IDs])

  dis_index_rwl <- Reduce(f = \(x, y) merge(x, y, by = "year", all = TRUE),
                          dis_index_rwl0[orig.IDs])

  # Make years as rownames & then remove the year column
  rownames(untransformed_rwl) <- untransformed_rwl$year
  untransformed_rwl <- untransformed_rwl[,orig.IDs]

  rownames(dis_index_rwl) <- dis_index_rwl$year
  dis_index_rwl <- dis_index_rwl[,orig.IDs]

  # Collapse the detection & removal iterations into a data.frame that contains all the detected
  # disturbances
  dist_det <- lapply(dist_iter, FUN = \(x) {
    x1 <- x[["Disturbance curves"]][orig.IDs]

    x2 <- lapply(x1, FUN = \(y) {
      if (!is.character(y)){
        y[1, c("year", "dir", "dur", "eq")]
      }
    }) |> do.call(what = "rbind")

    x2$series <- rownames(x2)
    x2[, c("series", "year", "dir", "dur", "eq")]

  }) |> do.call(what = "rbind")

  # Assign data classes
  class(untransformed_rwl) <- c("rwl", "data.frame")
  class(dis_index_rwl) <- c("rwl", "data.frame")

  ## Last steps are to output the main results (disturbance-free, disturbance index,
  # original series) and all of the iterations of the disturbance removal process, with the
  # trend curves, etc.
  # all of these can then be plotted in the plot_ci_detect() function

  ci_output_list <- list(untransformed_rwl, dis_index_rwl, dist_det, dist_iter, cp_out)
  names(ci_output_list) <- c("Disturbance-free series", "Disturbance index",
                             "Detected disturbances", "Disturbance removal iterations",
                             "Cook & Peters detrend")
  ci_output_list

} # End of the ci_detect() function
