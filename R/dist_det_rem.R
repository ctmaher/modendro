#' Supporting function for ci_detect - statistical detection and removal of disturbances in tree ring series
#'
#' @description
#' This function implements the disturbance detection and removal steps for the curve intervention detection techniques described by Druckenbrod et al. 2013,
#' Rydval et al. 2016, and Rydval et al. 2018. The function works on a collection of transformed ring width series. Iterations
#' are performed with \code{\link{ci_detect}}.
#'
#' @param rwi A rwl object of detrended residual series produced by \code{\link{cp_detrend}}.
#' @param min.win The minimum disturbance length in years (i.e., a moving window) to search for. The default is 9.
#' @param max.win The maximum disturbance length in years (i.e., a moving window) to search for. The default is 30.
#' @param thresh The disturbance detection threshold, corresponding to the number of deviations from the robust mean. The default is 3.29, following Druckenbrod et al. 2013.
#' @param dist.span Parameter to determine the wiggliness of the loess splines fit to disturbance periods (when Hugershoff fits fail). Higher numbers = more wiggles. Passes to \code{\link[stats]{loess}}. The default is 1.25.
#'
#' @details
#' The basic process that `dist_det_rem` performs is to take standardized tree ring series (specifically the output from the Cook & Peters (1997) process, as implemented by \code{\link{cp_detrend}} in \code{\link{ci_detect}}),
#' fit AR models (using the Burg method) to them and extract the residuals. The AR models are also "backcasted" to estimate residuals for the 1st n years which are NAs (n being equal to the AR order).
#' Backcasting involves reversing each series then fitting AR models of the same order as the original. Then, we compute moving averages of various window lengths. Disturbances are defined as moving average values that fall outside
#' the threshold of number of deviations from the robust mean (3.29), computed for each series and each moving average window length. The largest deviation - the largest absolute value of the difference between
#' the moving average value & the (negative or positive) threshold - for each series is selected first. Thus the largest magnitude disturbance could be a suppression or a release.
#' This defines the onset year and duration of the disturbance. A modified version of a Hugershoff curve (Warren & MacWilliam 1981) is then fit (via \code{\link[stats]{nls}}) to the disturbance period and beyond (all the way to the end of the series)
#' transformed, detrended residual series (not the AR residuals), with slightly different coefficients than the Hugershoff curve used by Rydval et al. Namely, the following values are fixed: \emph{b} = 1, and \emph{d} = 0. The Warren & MacWilliam curve has a \emph{t} coefficient which we use here.
#' These modifications allow for more robust fitting (reduced parameters require fewer degrees of freedom). The \emph{d} term acts like an intercept, and since all series have an overall mean of 0 (they are residuals),
#' this constraint is reasonable. Setting the \emph{b} term to 1 allows the beginning y values of the Hugershoff curve to go above or below 0 (more directly determined by \emph{t}), which helps in minimizing artifacts from poor fit in the early years of a disturbance period.
#'
#' While the Hugershoff curve fitting is reasonably robust, \code{\link[stats]{nls}} occasionally fails to converge. In these cases a \code{\link[stats]{loess}} spline is fit to the disturbance period only instead (not the whole remainder of the series as in the Hugershoff).
#' The wiggliness of these splines can be adjusted using the `dist.span` argument. For both methods, the resulting fitted curve is subtracted from the series. The recent value of the series before the disturbance period - a robust mean of the period before
#' the disturbance, or, if there is not an adequate period before (at least as long as the disturbance period window), the period after, is added back to the curve-series difference. This is done to avoid large artifacts from the disturbance removal process.
#' The output series with disturbances removed is now "corrected".
#'
#'
#' @return a 3-element list that contains the corrected RWI, the fitted disturbance curves, and data on disturbance detection iterations
#'
#' @references
#' Druckenbrod, D. L., N. Pederson, J. Rentch, and E. R. Cook. (2013) A comparison of times series approaches for dendroecological reconstructions of past canopy disturbance events.
#' \emph{Forest Ecology and Management}, \strong{302}, 23–33.
#'
#' Rydval, M., D. Druckenbrod, K. J. Anchukaitis, and R. Wilson. (2016) Detection and removal of disturbance trends in tree-ring series for dendroclimatology.
#' \emph{Canadian Journal of Forest Research}, \strong{401}, 387–401.
#'
#' Rydval, M., D. L. Druckenbrod, M. Svoboda, V. Trotsiuk, P. Janda, M. Mikoláš, V. Čada, R. Bače, M. Teodosiu, and R. Wilson. (2018) Influence of sampling and disturbance history on climatic sensitivity of temperature limited conifers.
#' \emph{The Holocene}, \strong{28}(10), 1574-1587.
#'
#' Cook, E. R., and Peters, K. (1997) Calculating unbiased tree-ring indices for the study of climatic and environmental change.
#' \emph{The Holocene}, \strong{7}(3), 361-370.
#'
#' Warren, W. G., and S. L. MacWilliam. 1981. Test of a new method for removing the growth trend from dendrochronological data.
#' \emph{Tree Ring Bulletin} \strong{41}, 55–66.
#'
#' @import dplR
#' @importFrom zoo rollmean
#' @import stats
#' @import DescTools
#' @export
#'
#'
#'

dist_det_rem <- function(rwi,
                         min.win = 9,
                         max.win = 30,
                         thresh = 3.29,
                         dist.span = 1.25,
                         var.type = "s_bi",
                         add.recent.rwi = TRUE
) {

  ## Error catching
  stopifnot("var.type is not valid. Must be 's_bi' or 'mad'." =
              var.type %in% "s_bi" |
              var.type %in% "mad")

  stopifnot("Value for thresh is not valid. Must be '>=1'." =
              thresh >= 1)

  ## Find the best AR model for each series
  ar_fits <-
    lapply(rwi, FUN = \(x) ar(
      x,
      method = "burg",
      aic = TRUE,
      na.action = na.omit
    ))

  # Get the resids
  ar_resids <- lapply(ar_fits, FUN = \(x) x$resid)

  # Reverse the cp series
  cp_rev <- lapply(rwi, rev)

  ## "Backcast" the NAs (due to the ar order lag) at the beginning of each series
  # Fit ar models of the same order as those above to the reversed data

  comp_resids <- mapply(
    FUN = \(x, y) {
      if (x$order > 0) {
        # If best fit order was 0, we don't need to do backcasting
        br <-
          ar(y,
             order.max = x$order,
             method = "burg",
             aic = FALSE)$resid |>
          rev() # flip it back to the correct chronological order
        # subset just the ones that need to be filled in & attach to beginning of rest of resids
        c(br[1:x$order], na.omit(x$resid)) |>
          as.numeric()
      } else {
        # If best fit order was 0, just return the AR resids (no backcasting).
        x$resid
      }
    },
    x = ar_fits,
    y = cp_rev,
    SIMPLIFY = FALSE
  )

  # add the names back (these are the years) as rownames of data.frames
  comp_resids <- mapply(
    FUN = \(x, y) {
      rn <- names(x)
      x <- as.data.frame(x)
      colnames(x) <- "value"
      rownames(x) <- rn
      x[, "value"] <- y
      x
    },
    x = rwi,
    y = comp_resids,
    SIMPLIFY = FALSE
  )

  ## Compute moving window averages of various lengths - the length of the window corresponds to disturbance length
  # the min.win sets the smallest possible size, and tradition has it that 20 is the longest period. We can
  # reassess that later - maybe the user can just specify this, within reasonable limits (max = 1/3 of series?).

  mov_avgs <- lapply(comp_resids, \(x) {
    # Cap max.win at 1/3 the series length. If not, detection becomes oversensitive for short series
    max.win <- min(max.win, round(nrow(x)/3))

    win_lens <- min.win:max.win
    ma_list <- lapply(win_lens, \(w) {
      # Extract the values from the data frame
      values <- x[,"value"]

      # Compute the moving average using rollmean from the zoo package
      ma <- rollmean(values,
                     k = w,
                     align = "left",
                     fill = NA)

      # Return the result
      data.frame(value = ma)
    })
    # Pass on the window lengths as names
    names(ma_list) <- win_lens
    ma_list
  })

  # The outer level of mov_avgs are each moving window length, with the series inside each of those

  ## Define upper and lower disturbance thresholds for each moving window length
  # The threshold is by default 3.29 robust scales from the robust mean of each moving average series
  tbrms <- lapply(mov_avgs, FUN = \(x) {
    lapply(x, FUN = \(x) {
      TukeyBiweight(x$value, const = 9, na.rm = TRUE)
      #tbrm(x$value, C = 9) # The dplR version
    })
  })

  # The robust scale is from Hoaglin et al 1983, p417.
  if (var.type %in% "s_bi") {
    var <- lapply(mov_avgs, FUN = \(x) {
      lapply(x, FUN = \(x) {
        s_bi(na.omit(x$value))$s_bi
      })
    })
  } else {
    if (var.type %in% "mad") {
      var <- lapply(mov_avgs, FUN = \(x) {
        lapply(x, FUN = \(x) {
          mad(x$value, na.rm = TRUE)
        })
      })
    }
  }

  lo_vals <- mapply(
    FUN = \(x, y) {
      mapply(
        FUN = \(a, b) {
          a - b * thresh
        },
        a = x,
        b = y,
        SIMPLIFY = FALSE
      )
    },
    x = tbrms,
    y = var,
    SIMPLIFY = FALSE
  )

  hi_vals <- mapply(
    FUN = \(x, y) {
      mapply(
        FUN = \(a, b) {
          a + b * thresh
        },
        a = x,
        b = y,
        SIMPLIFY = FALSE
      )
    },
    x = tbrms,
    y = var,
    SIMPLIFY = FALSE
  )


  ## Identify disturbances for each series & each moving window size, select the largest among window sizes
  # disturbance size is defined as the difference between the moving average and the threshold
  # There is potential to bias detection toward the smallest window sizes if I use only the "raw" value of the
  # moving average disturbances. The better metric needs to be relative to the threshold for each particular window size.
  # So that would be the abs of the difference between ma value & the (lo or hi) threshold.

  # 1st step is to determine which ma values are outside of lo.vals & hi.vals, if any, for all ma window sizes
  # These steps return index values for all disturbances
  pos_out <- mapply(
    FUN = \(x, y) {
      mapply(
        FUN = \(a, b) {
          # Find the index position of any disturbances
          index_pos <- which(a$value > b)

          # Compute how big the disturbance is
          dev <- abs(a$value[index_pos] - b)

          # Put the results in a data.frame
          dist_df <- data.frame(index_pos = index_pos, dev = dev)

          # return just the largest disturbance within each window size
          dist_df[which.max(dist_df$dev), ]

        },
        a = x,
        b = y,
        SIMPLIFY = FALSE
      )
    },
    x = mov_avgs,
    y = hi_vals,
    SIMPLIFY = FALSE
  )

  neg_out <- mapply(
    FUN = \(x, y) {
      mapply(
        FUN = \(a, b) {
          # Find the index position of any disturbances
          index_pos <- which(a$value < b)

          # Compute how big the disturbance is
          dev <- abs(a$value[index_pos] - b)

          # Put the results in a data.frame
          dist_df <- data.frame(index_pos = index_pos, dev = dev)

          # return just the largest disturbance within each window size
          dist_df[which.max(dist_df$dev), ]

        },
        a = x,
        b = y,
        SIMPLIFY = FALSE
      )
    },
    x = mov_avgs,
    y = lo_vals,
    SIMPLIFY = FALSE
  )


  # Now find the largest disturbances for each series among all window sizes and record the window size and the index value.
  # The index value will correspond to year.

  # First the whole list for all window lengths
  max_pos_outs <- lapply(pos_out, FUN = \(x) {
    all_win_df <- do.call("rbind", x)
    all_win_df$win_len <- as.numeric(rownames(all_win_df))
    if (nrow(all_win_df) > 0) {
      all_win_df$dist_dir <- "pos"
    }
    all_win_df
  })

  max_neg_outs <- lapply(neg_out, FUN = \(x) {
    all_win_df <- do.call("rbind", x)
    all_win_df$win_len <- as.numeric(rownames(all_win_df))
    if (nrow(all_win_df) > 0) {
      all_win_df$dist_dir <- "neg"
    }
    all_win_df
  })

  # Second the max among all window lengths (i.e., just 1 or no disturbance for each series)
  max_pos_out <- lapply(max_pos_outs, FUN = \(x) {
    x[which.max(x$dev), ]
  })

  max_neg_out <- lapply(max_neg_outs, FUN = \(x) {
    x[which.max(x$dev), ]
  })

  # Third find which is the largest of the pos and neg disturbances
  max_out <- mapply(
    FUN = \(x, y) {
      if (nrow(x) > 0 &
          nrow(y) > 0) {
        # If there are both pos & neg disturbances...
        if (x$dev > y$dev) {
          # Chose largest of the two
          x
        } else {
          y
        }
      } else {
        # if there aren't both
        if (nrow(x) > 0) {
          # if there is only a pos disturbance
          x
        } else {
          if (nrow(y) > 0) {
            # if there is only a neg disturbance
            y
          } else {
            "No disturbances detected"
          }
        }
      }

    },
    x = max_pos_out,
    y = max_neg_out,
    SIMPLIFY = FALSE
  )
  # max_out contains some basic disturbance info - the index of the starting year of the disturbance (ie, year),
  # the dev value of the disturbance (magnitude), the window length, and the direction of the disturbance.
  # index_pos & win_len determine the subset over which to fit a curve (dist_curves, below).
  # Can also use win_len to choose the relevant mov_avgs, tbrms, and lo_ & hi_vals out thresholds
  # The AR residual series can be extracted...

  dist_mov_avgs <- Map(
    f = \(a, b, c, d, e, x) {
      if (is.character(a)) {
        # For the series with no disturbances detected
        a
      } else {
        # create a data.frame that makes plotting easy (using ggplot2)
        win_len <- a[, "win_len"]
        # Residual series and moving averages are data.frames already
        ar_resid_vals <- b
        ar_resid_vals[, "type"] <- "AR residuals"
        ar_ma_vals <- c[[as.character(win_len)]]
        ar_ma_vals[, "type"] <- paste0(win_len, "-year mean")

        # The following are just single values - they must be repeated
        tbrm_vals <-
          data.frame(value = rep(d[[as.character(win_len)]], nrow(ar_ma_vals)),
                     type = "TBRM")
        lo_thresh_vals <-
          data.frame(value = rep(e[[as.character(win_len)]], nrow(ar_ma_vals)),
                     type = "Detection thresh.")
        hi_thresh_vals <-
          data.frame(value = rep(x[[as.character(win_len)]], nrow(ar_ma_vals)),
                     type = " ")

        # bind them together & return the result
        rbind(ar_resid_vals,
              ar_ma_vals,
              tbrm_vals,
              lo_thresh_vals,
              hi_thresh_vals)

      }

    },
    a = max_out,
    b = comp_resids,
    c = mov_avgs,
    d = tbrms,
    e = lo_vals,
    x = hi_vals
  )


  ## Fit curves to the largest disturbance from each series (the resid detrended series in cp_list),
  # subtract the disturbance curve, store the records of everything.
  dist_curves <- mapply(FUN = \(x, y) {
    # The index in max_out corresponds rows/elements in cp_list
    # Control for the ones with no disturbances
    if (is.character(y)) {
      dist_period <- y
    } else {

      # Make a data.frame from the series
      series_df <- as.data.frame(x)
      colnames(series_df) <- "rwi"
      series_df$year <-
        rownames(series_df) |>
        as.numeric()

      # isolate just the disturbance period
      dist_period <-
        series_df[y[, "index_pos"]:(y[, "index_pos"] + y$win_len - 1), ]

      # Attach the rest of series if there is one
      if ((max(dist_period$year) + 1) < max(series_df$year)) {
        rest_of_years <- (max(dist_period$year) + 1) : max(series_df$year)
        rest_of_series <- data.frame(rwi = series_df[series_df$year %in% rest_of_years, "rwi"],
                                     year = rest_of_years
        )
        # Rbind it
        dist_period <- rbind(dist_period, rest_of_series)
      }

      # Record the direction of the disturbance
      dist_period$dir <- y$dist_dir

      # Record the actual duration of the disturbance
      dist_period$dur <- y$win_len

      ## Curve fitting
      # Hugershoff - fits to the detected period and the reminder of the series too
      # The formula is modified. t is an added parameter that controls how far above/below the
      # initial fit can go beyond the asymptote. b = 1, always, to allow t to work. d = 0, always.
      # d mainly controls the asymptote value. We get a better chance at a successful fit if
      # we set these parameters here.
      hug_form0 <- formula(rwi ~ a * ((x - t)^1) * exp(-c*(x - t)) + 0)

      dist_period$x <- 1:(nrow(dist_period))

      # Set up some start values & constraints for a
      a_start <- ifelse(y$dist_dir %in% "pos", 0.1, -0.1)
      if (y$dist_dir %in% "pos") {
        a_const <- c(0.005, 5)
      } else {
        a_const <- c(-5, -0.005)
      }

      lower_const <- list(a = a_const[1],
                          c = -5,
                          t = -10)
      upper_const <- list(a = a_const[2],
                          c = 5,
                          t = 10)

      hug_fit <- try(
        nls(hug_form0,
            data = dist_period,
            start = list(a = a_start,
                         c = 0.1,
                         t = 1.5),
            algorithm = "port",
            lower = lower_const,
            upper = upper_const,
            control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = FALSE)),
        silent = TRUE)

      if (data.class(hug_fit) %in% "try-error") { # If the Hugershoff fit failed, fit a spline instead.
        # if this option, we should only fit & subtract the disturbance period itself
        dist_period <- dist_period[1:y$win_len,]
        spline_fit <-
          loess(rwi ~ year,
                data = dist_period,
                span = dist.span,
                family = "symmetric")

        dist_period$curve <- spline_fit$fitted
        dist_period$eq <- "loess spline"

      } else {
        dist_period$curve <- predict(hug_fit, newdata = dist_period)
        hug_coef <- coef(hug_fit) |> round(4)
        plus_minus <- ifelse(hug_coef[[3]] > 0, "-", "+")
        dist_period$eq <- paste0("y == ", hug_coef[[1]],
                                 " * (x ", plus_minus, " ", abs(hug_coef[[3]]), ")",
                                 " * e^(", -1*hug_coef[[2]],
                                 " * (x ", plus_minus, " ", abs(hug_coef[[3]]), "))")
      }

      # "Correct" the rwi values for the disturbance period by subtracting the fitted curve (aka, the residuals)
      # add back the recent value of the series before the disturbance period - a robust mean of the period before
      # the disturbance, or, if there is not an adequate period before, do the period after.

      if (add.recent.rwi == TRUE) {
        per_len <- ifelse(y$win_len < 15, 15, y$win_len)
        ba_dist <- ((min(dist_period$year, na.rm = TRUE) - 1) - per_len):(min(dist_period$year, na.rm = TRUE) - 1)
        if (min(ba_dist) < min(series_df$year) | # if there is no before or
            length(ba_dist[ba_dist %in% series_df$year]) < per_len){ # the before is < than the disturbance or 15 years
          # select the period after
          ba_dist <- (max(ba_dist)+1):(max(ba_dist) + per_len)
          #ba_dist <- (max(ba_dist) + 1):(max(ba_dist) + max(series_df$year))
        } # Otherwise, proceed with the before period already determined

        tbrm_recent_rwi <- series_df$rwi[series_df$year %in% ba_dist] |>
          TukeyBiweight()

      } else { # just use the raw difference for the mean of 0 if add.recent.rwi is FALSE
        tbrm_recent_rwi <- 0
      }

      dist_period$rwi.cor <-
        (dist_period$rwi - dist_period$curve) + tbrm_recent_rwi
    }
    # Return the results
    dist_period

  }, x = rwi, y = max_out, SIMPLIFY = FALSE)


  # Now the corrected rwi values can be inserted into the series to remove the disturbances
  rwi2 <- mapply(FUN = \(x, y) {
    # Control for the ones with no disturbances
    if (is.character(y)) {
      x # return the original series
    } else {
      x[which(names(x) %in% y[,"year"])] <- y[,"rwi.cor"] # substitute the corrected section
      x # Return the series
    }
  }, x = rwi, y = dist_curves, SIMPLIFY = FALSE)


  det_rem_dist_list <- list(rwi2, dist_curves, dist_mov_avgs)
  names(det_rem_dist_list) <-
    c("Corrected RWI", "Disturbance curves", "Disturbance detection")
  det_rem_dist_list

} # End of the dist_det_rem() function
