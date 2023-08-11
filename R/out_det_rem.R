#' Supporting function for `ci_detect()` - statistical identification and removal of outliers in tree ring series
#'
#' @description
#' This function implements the outlier detection and premoval steps for the curve intervention detection techniques described by Druckenbrod et al,
#' Rydval et al, and Rydval et al.
#'
#' @param rwi A rwl object of detrended standardized series produced by `cp_detrend()`.
#' @param min.win The minimum outlier length in years (i.e., a moving window) to search for. The default is 9.
#' @param max.win The maximum outlier length in years (i.e., a moving window) to search for. The default is 30.
#' @param out.span Parameter to determine the wiggliness of the loess splines fit to outlier periods. Higher numbers = more wiggles. Passes to `stats::loess()`. The default is 1.
#'
#' @import dplR
#' @importFrom zoo rollapply
#' @import stats
#' @import DescTools
#' @export

out_det_rem <- function(rwi,
                        min.win = 9,
                        max.win = 30,
                        var.type = "s_bi",
                        thresh = 3.29,
                        out.span = 1.25,
                        add.recent.rwi = TRUE
) {

  ## Error catching
  stopifnot("var.type is not valid. Must be 's_bi' or 'mad'." =
              var.type %in% "s_bi" |
              var.type %in% "mad")

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

  # mapply(FUN = \(x,y){
  #
  #   mess <- any(is.na(y))
  #   mess
  # }, x = ar_fits,
  # y = cp_rev)
  #
  # ar_fits[["2079"]]
  # rwi[["2079"]]
  # cp_rev
  # gb["2079"]

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

  ## Compute moving window averages of various lengths - the length of the window corresponds to outlier length
  # the min.win sets the smallest possible size, and tradition has it that 20 is the longest period. We can
  # reassess that later - maybe the user can just specify this, within reasonable limits (max = 1/3 of series?).

  mov_avgs <- lapply(comp_resids, \(x) {
    # Cap max.win at 1/3 the series length. If not, detection becomes oversensitive for short series
    max.win <- min(max.win, round(nrow(x)/3))

    win_lens <- min.win:max.win
    ma_list <- lapply(win_lens, \(w) {
      # Extract the values from the data frame
      values <- x$value

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

  ## Define upper and lower outlier thresholds for each moving window length
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


  ## Identify outliers for each series & each moving window size, select the largest among window sizes
  # Outlier size is defined as the difference between the moving average and the threshold
  # There is potential to bias detection toward the smallest window sizes if I use only the "raw" value of the
  # moving average outliers. The better metric needs to be relative to the threshold for each particular window size.
  # So that would be the abs of the difference between ma value & the (lo or hi) threshold.

  # 1st step is to determine which ma values are outside of lo.vals & hi.vals, if any, for all ma window sizes
  # These steps return index values for all outliers
  pos_out <- mapply(
    FUN = \(x, y) {
      mapply(
        FUN = \(a, b) {
          # Find the index position of any outliers
          index_pos <- which(a$value > b)

          # Compute how big the outlier is
          dev <- abs(a$value[index_pos] - b)

          # Put the results in a data.frame
          out_df <- data.frame(index_pos = index_pos, dev = dev)

          # return just the largest outlier within each window size
          out_df[which.max(out_df$dev), ]

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
          # Find the index position of any outliers
          index_pos <- which(a$value < b)

          # Compute how big the outlier is
          dev <- abs(a$value[index_pos] - b)

          # Put the results in a data.frame
          out_df <- data.frame(index_pos = index_pos, dev = dev)

          # return just the largest outlier within each window size
          out_df[which.max(out_df$dev), ]

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


  # Now find the largest outliers for each series among all window sizes and record the window size and the index value.
  # The index value will correspond to year.

  # First the whole list for all window lengths
  max_pos_outs <- lapply(pos_out, FUN = \(x) {
    all_win_df <- do.call("rbind", x)
    all_win_df$win_len <- as.numeric(rownames(all_win_df))
    if (nrow(all_win_df) > 0) {
      all_win_df$out_dir <- "pos"
    }
    all_win_df
  })

  max_neg_outs <- lapply(neg_out, FUN = \(x) {
    all_win_df <- do.call("rbind", x)
    all_win_df$win_len <- as.numeric(rownames(all_win_df))
    if (nrow(all_win_df) > 0) {
      all_win_df$out_dir <- "neg"
    }
    all_win_df
  })

  # Second the max among all window lengths (i.e., just 1 or no outlier for each series)
  max_pos_out <- lapply(max_pos_outs, FUN = \(x) {
    x[which.max(x$dev), ]
  })

  max_neg_out <- lapply(max_neg_outs, FUN = \(x) {
    x[which.max(x$dev), ]
  })

  # Third find which is the largest of the pos and neg outliers
  max_out <- mapply(
    FUN = \(x, y) {
      if (nrow(x) > 0 &
          nrow(y) > 0) {
        # If there are both pos & neg outliers...
        if (x$dev > y$dev) {
          # Chose largest of the two
          x
        } else {
          y
        }
      } else {
        # if there aren't both
        if (nrow(x) > 0) {
          # if there is only a pos outlier
          x
        } else {
          if (nrow(y) > 0) {
            # if there is only a neg outlier
            y
          } else {
            "No outliers detected"
          }
        }
      }

    },
    x = max_pos_out,
    y = max_neg_out,
    SIMPLIFY = FALSE
  )
  # max_out contains some basic outlier info - the index of the starting year of the outlier (ie, year),
  # the dev value of the outlier (magnitude), the window length, and the direction of the outlier.
  # index_pos & win_len determine the subset over which to fit a curve (out_curves, below).
  # Can also use win_len to choose the relevant mov_avgs, tbrms, and lo_ & hi_vals out thresholds
  # The AR residual series can be extracted...

  out_mov_avgs <- Map(
    f = \(a, b, c, d, e, x) {
      if (is.character(a)) {
        # For the series with no outliers detected
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


  ## Fit curves to the largest outlier from each series (the resid detrended series in cp_list),
  # subtract the outlier curve, store the records of everything.
  out_curves <- mapply(FUN = \(x, y) {
    # The index in max_out corresponds rows/elements in cp_list
    # Control for the ones with no outliers
    if (is.character(y)) {
      out_period <- y
    } else {

      # Make a data.frame from the series
      series_df <- as.data.frame(x)
      colnames(series_df) <- "rwi"
      series_df$year <-
        rownames(series_df) |>
        as.numeric()

      # isolate just the outlier period
      out_period <-
        series_df[y[, "index_pos"]:(y[, "index_pos"] + y$win_len - 1), ]

      # Add a rest of series if there is one
      if ((max(out_period$year) + 1) < max(series_df$year)) {
        rest_of_years <- (max(out_period$year) + 1) : max(series_df$year)
        rest_of_series <- data.frame(rwi = series_df[series_df$year %in% rest_of_years, "rwi"],
                                     year = rest_of_years
        )
        # Rbind it
        out_period <- rbind(out_period, rest_of_series)
      }

      # Record the direction of the disturbance
      out_period$dir <- y$out_dir

      # Record the actual duration of the disturbance
      out_period$dur <- y$win_len

      ## Curve fitting

      #if (fit.type == "Hugershoff") {
      # Hugershoff - fits to the detected period and the reminder of the series too
      # The formula is modified. z is an added parameter that controls how far above/below the
      # initial fit can go beyond the asymptote. b = 1, always, to allow z to work. d = 0, always.
      # d mainly controls the asymptote value. We get a better chance at a successful fit if
      # we set these parameters here.
      hug_form0 <- formula(rwi ~ a * ((x - z)^1) * exp(-c*(x - z)) + 0)

      out_period$x <- 1:(nrow(out_period))

      #hug_fit <- NULL
      hug_fit <- try(
        nls(hug_form0,
            data = out_period,
            start = list(a = ifelse(y$out_dir %in% "pos", 0.1, -0.1),
                         c = 0.1,
                         z = 1.5),
            control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = FALSE)),
        silent = TRUE)

      if (data.class(hug_fit) %in% "try-error") { # If the Hugershoff fit failed, fit a spline instead.
        # if this option, we should only plot & subtract the disturbance period itself
        out_period <- out_period[1:y$win_len,]
        spline_fit <-
          loess(rwi ~ year,
                data = out_period,
                span = out.span,
                family = "symmetric")

        out_period$curve <- spline_fit$fitted

      } else {
        out_period$curve <- predict(hug_fit, newdata = out_period)
      }

      # "Correct" the rwi values for the outlier period by subtracting the fitted curve (aka, the residuals)
      # add back the recent value of the series before the outlier period - a robust mean of the period before
      # the disturbance, or, if there is not an adequate period before, do the period after.

      if (add.recent.rwi == TRUE) {
        ba_dist <- (min(out_period$year, na.rm = TRUE) - y$win_len):min(out_period$year, na.rm = TRUE)
        if (min(ba_dist) < min(series_df$year) | # if there is no before or the before is < than the disturbance
            length(ba_dist[ba_dist %in% series_df$year]) < y$win_len){
          # select the period after
          ba_dist <- (max(ba_dist)+1):(max(ba_dist) + y$win_len)
          out_period$rwi.cor <- # or just use the raw difference for the mean of 0
            out_period$rwi - out_period$curve
        } else {# Otherwise, proceed with the before period already determined
          tbrm_recent_rwi <- series_df$rwi[series_df$year %in% ba_dist] |>
            TukeyBiweight()

          out_period$rwi.cor <-
            (out_period$rwi - out_period$curve) + tbrm_recent_rwi
        }
      } else {
        out_period$rwi.cor <-
          (out_period$rwi - out_period$curve)
      }

    }
    # Return the results
    out_period

  }, x = rwi, y = max_out)


  # Now the corrected rwi values can be inserted into the series to remove the outliers
  rwi2 <- mapply(FUN = \(x, y) {
    # Control for the ones with no outliers
    if (is.character(y)) {
      x # return the original series
    } else {
      x[which(names(x) %in% y$year)] <- y$rwi.cor # substitute the corrected section
      x # Return the series
    }
  }, x = rwi, y = out_curves)

  #rwi2

  det_rem_out_list <- list(rwi2, out_curves, out_mov_avgs)
  names(det_rem_out_list) <-
    c("Corrected RWI", "Outlier curves", "Outlier detection")
  det_rem_out_list

} # End of the out_det_rem() function
