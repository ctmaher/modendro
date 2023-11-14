#' Detrend raw ring widths using Cook & Peters' (1997) method.
#'
#' @description
#' Compute variance-stabilized & residual detrended ring width indices that minimize the series end effect problems
#' that can result from using ratios to derive detrended indices.
#'
#' @param rwl A rwl-type data.frame (e.g., read in by \code{\link[dplR]{read.rwl}}). Essentially a data.frame with columns names as series IDs and years as rownames.
#' @param detrend.method Character string of the detrending method to use. Passes to \code{\link[dplR]{detrend}}.
#' @param nyrs Numeric vector, used in dplR's \code{\link[dplR]{detrend}} function for the `"Spline"` and `"AgeDepSpline"` methods.
#' @param pos.slope Should positive slopes be allowed in the detrending curves? Generally this should be FALSE (the default),
#' but when used in \code{\link{ci_detect}} it is TRUE to detect deviations from any long term trend. Passes to \code{\link[dplR]{detrend}}.
#' @param standardize Logical vector indicating whether to standardize the output detrended & transformed residuals (default is TRUE). See details for more information.
#'
#' @details
#' For decades, the most common method of removing long-term size/age trends from tree ring width series has been to fit a curve
#' (e.g., a negative exponential equation) to the data, then divide the raw data by the predicted (trend) values.
#' The resulting indices have ± stable variance for their duration and a mean of 1, and are ready to be aggregated
#' together into chronologies with other series standardized via the same method. However, this method can introduce artificial
#' trends in locations along the series where the trend line has a poor fit with the data, especially when the predicted values drop below 0.5.
#' This is common at the end of series. The risk is that the researcher might obtain biased estimates of growth trends - i.e., enhanced growth
#' increases or decreases.
#'
#' Cook & Peters (1997) proposed the method employed here as a fix this problem.
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
#' The `standardize` argument indicates whether a standardization of the residual detrended series should be performed. This is
#' an important consideration because of the variable transformations applied here can result in series with dramatically different
#' variances - which will influence the any tree means or chronologies generated from these data. The standardization method we apply here is the same type found in
#' the `dplR` function \code{\link[dplR]{detrend}} with the `method = "Ar"` option. This is done by adding the means of the original ring width series
#' to the residual detrended series, then dividing the resulting series by their means. The results are series with a mean of 1 and similar variance.
#'
#' @return A named list with the following elements:
#' 1) "Resid. detrended series" - the residual detrended power transformed series.
#' 2) "Detrend curves" - the curves fitted to the transformed series.
#' 3) "Transformed ring widths" - the power transformed series (before detrending).
#' 4) "Transformation metadata" - a list of data.frames for each series about the type of transformation applied.
#' 5) "Detrending metadata" - a list of data.frames for each series about the detrending applied (from \code{\link[dplR]{detrend}}).
#' 6) "Raw ring widths" - the original rwl-data.frame
#'
#' @references
#' Cook, E. R., and Peters, K. (1997) Calculating unbiased tree-ring indices for the study of climatic and environmental change.
#' \emph{The Holocene}, \strong{7}(3), 361-370.
#'
#' @seealso \code{\link{pwr_t_rwl}}, \code{\link{find_opt_pwr}}
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
#' # Apply ggsave() to each plot:
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
           detrend.method = "Mean",
           nyrs = NULL,
           pos.slope = FALSE,
           standardize = TRUE) {
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

    if (detrend.method %in% c("none", "None", "Mean")) {
      detrend.method <- "Mean"
      message(
        "Proceding with the default of no detrending - series will be transformed and their means subtracted.
            Make sure you want this, otherwise chose a detrending method wisely"
      )
    }

    # catch any negative values here and set to 0
    if (any(rwl < 0, na.rm = TRUE)) {
      rwl[which(rwl < 0), ] <- 0
      warning(
        "One or more negative values detected in the rwl file. These were replaced with 0s.
              Check your data to make sure this is correct (e.g., sometimes -9999 is used as a place-holder value for 0 rings."
      )
    }

    orig.IDs <- names(rwl) # capture the original names in their original order

    # Find the optimal power of transformation and transform the series
    trans.list <- pwr_t_rwl(rwl) # Returns a list
    trans <- trans.list[["Transformed ring widths"]][,orig.IDs] # The transformed series
    mess.df <- trans.list[["Transformation metadata"]] # Info about transformations



    # detrend the transformed series
    detr.result <-
      dplR::detrend(
        trans + 2, # add 2 here so that fits are more likely to be positive values
        method = detrend.method,
        make.plot = FALSE,
        difference = TRUE,
        nyrs = nyrs,
        return.info = TRUE,
        pos.slope = pos.slope
      )

    # Extract the detrending info.
    detr.info <- detr.result[["model.info"]][orig.IDs]
    detr.info <- Map(f = \(d, n) {
      df <- do.call("rbind", d) |> as.data.frame()
      df[, "method"] <- names(d)
      df[, "series"] <- n
      df
    },
    d = detr.info[orig.IDs],
    n = names(detr.info[orig.IDs]))


    curv <- detr.result$curves - 2 # subtract the 2 we added above
    detr <- detr.result$series

    if (standardize == TRUE) {
      orig.means <- colMeans(rwl, na.rm = TRUE)
      resid.plus.orig <- sweep(detr, 2, orig.means, "+")
      resid.plus.orig.means <- colMeans(resid.plus.orig, na.rm = TRUE)
      # Control for the potential for negative series means. Probably impossible with tree ring data.
      resid.plus.orig.means[resid.plus.orig.means < 0] <- ((resid.plus.orig.means[resid.plus.orig.means < 0]) * -1)
      detr <- sweep(resid.plus.orig, 2, resid.plus.orig.means, "/")
    }

    out.list <- list(detr,
                     curv,
                     trans,
                     mess.df,
                     detr.info,
                     rwl)
    names(out.list) <-
      c(
        "Resid. detrended series",
        "Detrend curves",
        "Transformed ring widths",
        "Transformation metadata",
        "Detrending metadata",
        "Raw ring widths"
      )

    out.list

  } ## End of function
