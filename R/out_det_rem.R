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
#' @import zoo
#' @import stats
#' @import DescTools
#' @export

out_det_rem <- function(rwi,
                        min.win = 9,
                        max.win = 30,
                        span = 1) {
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

  ## Compute moving window averages of various lengths - the length of the window corresponds to outlier length
  # the min.win sets the smallest possible size, and tradition has it that 20 is the longest period. We can
  # reassess that later - maybe the user can just specify this, within reasonable limits (max = 1/3 of series?).

  mov_avgs <- lapply(comp_resids, \(x) {
    # Cap max.win at 1/3 the series length. If not, detection becomes oversensitive for short series
    max.win <- min(max.win, round(nrow(x) / 3))

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

  # The robust scale is from Hoaglin et al 1983, p417. Here I use the MAD as a simple substitute for now.
  mads <- lapply(mov_avgs, FUN = \(x) {
    lapply(x, FUN = \(x) {
      mad(x$value, na.rm = TRUE)
    })
  })

  lo_vals <- mapply(
    FUN = \(x, y) {
      mapply(
        FUN = \(a, b) {
          a - b * 3.29
        },
        a = x,
        b = y,
        SIMPLIFY = FALSE
      )
    },
    x = tbrms,
    y = mads,
    SIMPLIFY = FALSE
  )

  hi_vals <- mapply(
    FUN = \(x, y) {
      mapply(
        FUN = \(a, b) {
          a + b * 3.29
        },
        a = x,
        b = y,
        SIMPLIFY = FALSE
      )
    },
    x = tbrms,
    y = mads,
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
          out_df[which.max(out_df$dev),]

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
          out_df[which.max(out_df$dev),]

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
    x[which.max(x$dev),]
  })

  max_neg_out <- lapply(max_neg_outs, FUN = \(x) {
    x[which.max(x$dev),]
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
      # extract the outlier period and fit a curve
      series_df <- data.frame(x)
      colnames(series_df) <- "rwi"
      series_df$year <-
        rownames(series_df) # keep this as a character for now
      out_period <-
        series_df[y[, "index_pos"]:(y[, "index_pos"] + y$win_len - 1),]

      # Fit loess splines
      # Note: you can adjust the weight of each point using the weights argument
      # This could help with mimicking the the Hugershoff-type curves (which are more flexible at the start).
      # It might also make more sense to weight the last value more so that the residual is minimized.
      # This would limit sharp jumps in the resulting series.
      wts <- rep(1, nrow(out_period))
      # wts[nrow(out_period)] <- 4
      curve <-
        loess(rwi ~ year,
              data = out_period,
              span = span,
              weights = wts)
      out_period$curve <- curve$fitted
      # "Correct" the rwi values for the outlier period by subtracting the fitted curve (aka, the residuals)
      out_period$rwi.cor <- curve$residuals
      out_period$dir <- y$out_dir
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
      x[y$year] <- y$rwi.cor # substitute the corrected section
      x # Return the series
    }
  }, x = rwi, y = out_curves)

  #rwi2

  det_rem_out_list <- list(rwi2, out_curves, out_mov_avgs)
  names(det_rem_out_list) <-
    c("Corrected RWI", "Outlier curves", "Outlier detection")
  det_rem_out_list

} # End of the out_det_rem() function
