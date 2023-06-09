#' Detrend raw ring widths using Cook & Peters' (1997) method.
#'
#' @description
#' Compute variance-stabilized & residual detrended ring width indices that minimize the series end effect problems
#' that can result from using ratios to derive detrended indices.
#'
#'
#' @details
#' For decades, the default method of removing long-term size/age trends from tree ring width series has been to fit a curve
#' (e.g., a negative exponential equation) to the data, then divide the raw data by the predicted (trend) values.
#' The resulting indices have ± stable variance for their duration and a mean of 1, and are ready to be aggregated
#' together into chronologies with other series standardized via the same method. However, this method can introduce artificial
#' trends in locations along the series where the trend line has a poor fit with the data, especially when the predicted values drop below 0.5.
#' This is common at the end of series. The risk is that the researcher might obtain biased estimates of growth trends - i.e., enhanced growth
#' increases or decreases.
#'
#' Cook and Peters (1997) proposed the method employed here as a fix this problem.
#' The basic steps are to estimate the optimal power of transformation via the local spread versus level relationship of each series,
#' where local spread is defined as the absolute value of the first
#' differences, S, (|rwt - rwt-1|) and the local level is the arithmetic mean of each pair of adjacent values, M, (rwt + rwt-1)/2.
#' The spread versus level relationship is then modeled in a simple linear regression as log10(S) ~ log10(M). The optimal power of transformation
#' is then estimated as p = |1-slope| of this regression model. If p ≤ 0.1 (i.e., slope is near 1), a log10() transformation is applied.
#' If p > 1 (i.e., slope is near 0), no transformation is applied. Otherwise, p is applied to all the raw ring widths as rw^p.
#'
#' With no special accommodation, this power transformation process can produce odd results for 0 values (missing rings).
#' Log transformation will obviously fail in these circumstances. To reduce these problems, 0 ring width values are replaced
#' with the minimum value possible given the resolution of the data, which for tree ring data is typically 0.01 or 0.001 mm.
#' If you have many 0 ring width values in your data, make sure to take a look at the output plots (you should do this anyway) to
#' check for any weirdness.
#'
#'
#' @param rwl A rwl object (read in by dplR's `read.rwl()`)
#' @param detrend.method Character string of the detrending method to use. Passes to dplR's `detrend()` function.
#' @param nyrs Numeric vector, used in dplR's `detrend()` function for the "Spline" and "AgeDepSpline" methods.
#'
#' @return A list with the following elements:
#' 1) the residual detrended power transformed series.
#' 2) the curves fitted to the transformed series.
#' 3) the power transformed series (before detrending).
#' 4) a data.frame containing each series about the type of transformation applied.
#' If no detrending is performed, then only elements 3-5 are in the output list.
#'
#' @references
#' Cook, E. R., and Peters, K. (1997) Calculating unbiased tree-ring indices for the study of climatic and environmental change.
#' \emph{The Holocene}, \strong{7}(3), 361-370.
#'
#' @import dplR
#' @import tidyr
#' @export
#'
#' @examples
#' library(dplR)
#' # Bristlecone pine tree ring collection from Campito mountain, White Mountains, CA
#' # Many 0 value rings in this collection
#' data("ca533")
#' ca533_cp <- cp_detrend(rwl = ca533,
#'                        detrend.method = "AgeDepSpline",
#'                        nyrs = 50)
#' ca533_cp[[1]][[1]] # The output of the first series
#'
#' # You can conveniently make plots using the plot_cp_detrend() function
#' ca533_cp_plots <- plot_cp_detrend(ca533_cp)
#' # You can also save these plots (they are ggplot objects)
#' # to disk with something like this:
#' # Create a new directory for the plots (there could be a lot!)
#' cp_plot_dir <- paste0(getwd(),"/CP_plots")
#' dir.create(cp_plot_dir)
#' # Apply `ggsave()` to each plot:
#' mapply(FUN = \(x, y) {
#'         ggsave(filename = paste0(cp_plot_dir, "/", y, ".pdf"),
#'         plot = x,
#'         device = "pdf",
#'         width = 20,
#'         height = 14,
#'         units = "cm")
#' }, x = ca533_cp_plots, y = names(ca533_cp_plots))

cp_detrend <-
  function(rwl,
           detrend.method = "none",
           nyrs = NULL) {
    # Power transform series prior to detrending using methods of Cook and Peters (1997)

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

    stopifnot(
      "detrend.method is not a known option. See ?dplR::detrend.series for options" =
        detrend.method %in% c(
          "none",
          "None",
          "Spline",
          "ModNegExp",
          "Mean",
          "Ar",
          "Friedman",
          "ModHugershoff",
          "AgeDepSpline"
        )
    )

    if (detrend.method %in% c("none", "None")) {
      message(
        "Proceding with the default of no detrending - series will be transformed only.
            Make sure you want this, otherwise chose a detrending method wisely"
      )
    }

    if (any(rwl < 0, na.rm = TRUE)) {
      rwl[which(rwl < 0),] <- NA
      warning(
        "One or more negative values detected in the rwl file. These were replaced with NAs.
              Check your data to make sure this is correct (e.g., sometimes -9999 is used as a place-holder value for 0 rings."
      )
    }

    # Find the optimal power of transformation and transform the series
    trans.list <- pwr_t_rwl(rwl) # Returns a list
    trans <- trans.list[[1]] # The transformed series
    mess.df <- trans.list[[2]] # Info about transformations

    # detrend the transformed series -or- don't
    if (!c(is.null(detrend.method) |
           detrend.method %in% "None" |
           detrend.method %in% "none")) {
      detr.result <-
        dplR::detrend(
          trans,
          method = detrend.method,
          make.plot = FALSE,
          difference = TRUE,
          nyrs = nyrs,
          return.info = TRUE
        )

      curv <- detr.result$curves
      detr <- detr.result$series


      out.list <- list(detr,
                       curv,
                       trans,
                       mess.df,
                       rwl)
      names(out.list) <-
        c(
          "Resid. detrended series",
          "Detrend curves",
          "Transformed ring widths",
          "Metadata about transformations",
          "Raw ring widths"
        )

    } else {
      # i.e., no detrending

      out.list <- trans.list

    }
    out.list

} ## End of function


