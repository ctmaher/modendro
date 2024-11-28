#' Calculate multiple growth-climate relationships
#'
#' @description
#' This function is a workhorse function behind \code{\link{n_mon_corr}}. It splits and merges the
#' climate data and the rwl data, then calculates all the correlations between them.
#'
#' Internal function for modendro.
#'
#' @param rwl.group A rwl-type data.frame (e.g., read in by \code{\link[dplR]{read.rwl}}).
#' Essentially a data.frame with columns names as series IDs and years as rownames.
#' @param clim.group a `data.frame` with at least 3 columns: year, month (numeric), and a
#' climate variable.
#' @param clim.var character vector - the colname of the climate variable of interest in the `clim`
#' data.frame.
#' @param common.years numeric vector - a sequence of years to subset the clim and rwl data before
#' running correlation analyses. Must be at least 25 years long and years must exist in clim and
#' rwl.
#' @param agg.fun character vector specifying the function to use for aggregating monthly
#' climate combinations. Options are "mean" or "sum", e.g., for temperature or precipitation data,
#' respectively. Default is "mean".
#' @param max.win integer vector specifying how long, in months, the longest (or widest) moving
#' window is. Values limited to between 2 and 12. Default is 6 months.
#' @param win.align the alignment of the moving windows. Options are "left" or "right". If "left",
#' month will indicate the starting month of each moving window and NAs will appear at the end of
#' the series. If "right", month will indicate the ending month of each moving window and NAs appear
#' at the beginning. The default is "right".
#' @param max.lag integer vector specifying how many years of lag to calculate calculations for.
#' Default (and minimum) is 1 year.
#' @param hemisphere a character vector specifying which hemisphere your tree ring data - &
#' climate data - comes from ("N" or "S"). Conventions for assigning growth years - and thus
# aligning tree ring and climate data - are different for N and S hemisphere (see Details below).
# In the current version, this simply adds a "lag +1" year to the correlation calculations.
#' @param prewhiten logical vector specifying whether or not to convert tree ring & climate time
#' series to ARIMA residuals (aka "prewhitening"). A "best fit" ARIMA model is automatically
#' selected using \code{\link[forecast]{auto.arima}}.
#' This removes most autocorrelation in a time series, leaving only the high-frequency variation.
#' This is common practice before using standard methods for cross-correlations. Default is TRUE.
#' @param corr.method character vector specifying which correlation method to use. Options are
#' `c("spearman", "kendall", "pearson")`. Default is `"spearman"` using the
#' \code{\link[corTESTsrd]{corTESTsrd}} function (also used for `corr.method = "kendall"`). This
#' method reduces the type I error rate associated with autocorrelated series. CAUTION: Currently
#' `corr.method = "pearson"` doesn't make any adjustments for autocorrelation. See Details below.
#' @param group.var a character vector specifying a grouping variable (e.g., site, plot).
#' Must match a column in clim.


multi_clim_gro_corr <- function(rwl.group = NULL,
                                clim.group = NULL,
                                clim.var = NULL,
                                common.years = NULL,
                                agg.fun = "mean",
                                max.win = 6,
                                win.align = "right",
                                max.lag = 1,
                                prewhiten = TRUE,
                                hemisphere = NULL,
                                corr.method = "spearman",
                                group.var = NULL) {
  ## Compute the moving windows
  clim1 <- moving_win_multi(
    clim.group,
    clim.var = clim.var,
    win.lens = 2:max.win,
    win.align = win.align,
    agg.fun = agg.fun
  )

  if (!is.null(group.var)) {
    clim1[, group.var] <- unique(clim.group[, group.var])
  }

  # Split by month
  clim1.mo.split <- split(clim1, f = clim1$month)


  # Now is a good time to get the arima residuals - before we create copies of the data when we
  # set up the lagged datasets.
  if (prewhiten == TRUE) {
    # Climate data first
    clim1.mo.split <- lapply(clim1.mo.split, FUN = \(mo) {
      # Split into the different win.lens
      win.split <- split(mo, f = mo$win.len)
      lapply(win.split, FUN = \(this.win) {
        # Make sure we have years in order!
        this.win <- this.win[order(this.win[, "year"], decreasing = FALSE), ]
        arima.resids <- forecast::auto.arima(this.win[, clim.var], seasonal = FALSE) |>
          residuals() |>
          na.omit() |>
          as.numeric()
        # Get the NAs
        notNAs <- which(!is.na(this.win[, clim.var]))
        # put the resids into the "body"
        this.win[min(notNAs):max(notNAs), clim.var] <- arima.resids
        this.win
      }) |> do.call(what = "rbind")

    })
    # Now the tree ring data
    rwl.group <- rwl.group[order(rwl.group[, "year"], decreasing = FALSE), ]
    rwl.arima.resids <- apply(
      rwl.group[, !(colnames(rwl.group) %in% "year"), drop = FALSE],
      MARGIN = 2,
      FUN = \(this.series) {
        these.vals <- which(!is.na(this.series))
        trim.series <- this.series[min(these.vals):max(these.vals)]
        these.resids <- forecast::auto.arima(trim.series,
                                             seasonal = FALSE) |>
          residuals() |>
          as.numeric()
        this.series[min(these.vals):max(these.vals)] <- these.resids
        this.series
      },
      simplify = FALSE
    ) |> do.call(what = "cbind")
    # Put the data back together
    rwl.group <- cbind(rwl.arima.resids, rwl.group[, "year", drop = FALSE]) |>
      as.data.frame()
    colnames(rwl.group)[ncol(rwl.group)] <- "year"
    rownames(rwl.group) <- rwl.group[, "year"]
  }

  ## Set up the lags
  # basic idea - add n lag to actual year, then merge tree-ring & clim data based on year and
  # lagn.year, respectively
  # It is faster to rbind these together & do only a few merges over larger data.frames
  # than to merge over and over again for each month.

  if (hemisphere == "S") {
    # The lags are a bit different for the S Hemisphere - we need a lag+1 to capture all
    # the relevant climate periods

    clim.lag.bind <- lapply(clim1.mo.split, FUN = \(mo) {
      lag.years <- sapply(-1:max.lag, FUN = \(n.lag) {
        this.lag <- data.frame(mo[, "year"] + n.lag)
        colnames(this.lag) <- paste0("lag", n.lag, ".year")
        this.lag
      }) |> do.call(what = "cbind")
      cbind(mo, lag.years)
    }) |> do.call(what = "rbind")

    # Merge all the series at once for each lag

    lag.list <- sapply(-1:max.lag, FUN = \(n.lag) {
      merge(
        clim.lag.bind,
        rwl.group,
        by.x = paste0("lag", n.lag, ".year"),
        by.y = "year"
      )
    }, simplify = FALSE)
    names(lag.list) <- c(paste0("lag+", -1 * (-1:0)), paste0("lag-", 1:max.lag))

  } else {
    # For northern hemisphere

    clim.lag.bind <- lapply(clim1.mo.split, FUN = \(mo) {
      lag.years <- sapply(0:max.lag, FUN = \(n.lag) {
        this.lag <- data.frame(mo[, "year"] + n.lag)
        colnames(this.lag) <- paste0("lag", n.lag, ".year")
        this.lag
      }) |> do.call(what = "cbind")
      cbind(mo, lag.years)
    }) |> do.call(what = "rbind")

    # Merge all the series at once for each lag

    lag.list <- sapply(0:max.lag, FUN = \(n.lag) {
      merge(
        clim.lag.bind,
        rwl.group,
        by.x = paste0("lag", n.lag, ".year"),
        by.y = "year"
      )
    }, simplify = FALSE)
    names(lag.list) <- c(paste0("lag+", 0), paste0("lag-", 1:max.lag))

  }

  # set up all the possible combinations of series and moving window lengths
  all.corr.combos <- expand.grid(series.name = colnames(rwl.group)[!(colnames(rwl.group)
                                                                     %in% "year")],
                                 win.len = unique(clim.lag.bind[, "win.len"]))


  ## Run the correlations, which involves several layers of splitting & applying
  cor.res.df <- mapply(
    FUN = \(lag.x, which.lag) {
      # Need to split by month again
      lag.x.mo.split <- split(lag.x, f = lag.x[, "month"])

      # Take each month split and do all the possible correlations
      this.lag.corr <- lapply(lag.x.mo.split, FUN = \(mo) {
        # Take each possible correlation out of mo and run the correlation
        apply(all.corr.combos, MARGIN = 1, FUN = \(combo) {
          series.name <- combo[["series.name"]]
          win.len <- combo[["win.len"]]

          # The correlations
          # This subset eliminates NA values in either series.
          # NAs in the middle of the series should be a rare occurrence
          subs.for.corr <- mo[mo[, "win.len"] == win.len, c("year",
                                                            "month",
                                                            "win.len",
                                                            as.character(series.name),
                                                            clim.var)] |>
            na.omit()

          # Subset based on the common years specified by the user
          subs.for.corr <- subs.for.corr[subs.for.corr$year %in% common.years,]

          # Control for ties by simply removing them when they occur. This should be
          # relatively rare.
          subs.for.corr <-
            subs.for.corr[which(!duplicated(subs.for.corr[, clim.var])), ]
          subs.for.corr <-
            subs.for.corr[which(!duplicated(subs.for.corr[, as.character(series.name)])), ]


          if (corr.method == "pearson") {
            cor.res <- cor.test(
              subs.for.corr[, clim.var],
              subs.for.corr[, as.character(series.name)],
              method = "pearson",
              alternative = "two.sided"
            )
            out.df <- data.frame(
              series = as.character(series.name),
              month = unique(subs.for.corr[, "month"]),
              win.len = unique(subs.for.corr[, "win.len"]),
              coef = cor.res[["estimate"]][[1]],
              p = cor.res[["p.value"]][[1]],
              dir = ifelse(cor.res[["estimate"]][[1]] < 0, "Neg.", "Pos."),
              min.year = min(na.omit(subs.for.corr[, "year"])),
              max.year = max(na.omit(subs.for.corr[, "year"])),
              overlap = length(min(na.omit(subs.for.corr[, "year"])):
                                 max(na.omit(subs.for.corr[, "year"]))),
              clim.var = clim.var,
              corr.method = corr.method,
              hemisphere = hemisphere

            )
            out.df[, group.var] <- unique(mo[, group.var])
            out.df

          } else {
            # for "spearman" or "kendall" rank methods
            cor.res <- corTESTsrd::corTESTsrd(
              subs.for.corr[, clim.var],
              subs.for.corr[, as.character(series.name)],
              method = corr.method,
              # Add a switch in here for prewhitening?
              # Perhaps it is best to always be safe here.
              iid = FALSE,
              alternative = "two.sided"
            )
            out.df <- data.frame(
              series = as.character(series.name),
              month = unique(mo[, "month"]),
              win.len = unique(subs.for.corr[, "win.len"]),
              coef = cor.res[["rho"]],
              p = cor.res[["pval"]],
              dir = ifelse(cor.res[["rho"]] < 0, "Neg.", "Pos."),
              min.year = min(na.omit(subs.for.corr[, "year"])),
              max.year = max(na.omit(subs.for.corr[, "year"])),
              overlap = length(min(na.omit(subs.for.corr[, "year"])):
                                 max(na.omit(subs.for.corr[, "year"]))),
              clim.var = clim.var,
              corr.method = corr.method,
              hemisphere = hemisphere
            )
            out.df[, group.var] <- unique(mo[, group.var])
            out.df

          }

        }) |> do.call(what = "rbind")
      }) |> do.call(what = "rbind")

      this.lag.corr$lag <- ifelse(substr(unique(which.lag), 5, 6) == "0", "0",
                                  substr(unique(which.lag), 4, 6))
      this.lag.corr

    },
    lag.x = lag.list,
    which.lag = names(lag.list),
    SIMPLIFY = FALSE
  ) |> do.call(what = "rbind")

  # Assign a class to the rwl data
  rwl.dat <- rwl.group[, !(colnames(rwl.group) %in% "year"), drop = FALSE]
  class(rwl.dat) <- c("rwl", "data.frame")

  if (prewhiten == TRUE) {
    res.list <- list(cor.res.df,
                     do.call(what = "rbind", clim1.mo.split),
                     clim1,
                     rwl.dat)
    names(res.list) <- c("cor.res.dat", "clim.dat.pw", "clim.dat", "rwl.dat")
  } else {
    res.list <- list(cor.res.df, clim1, rwl.dat)
    names(res.list) <- c("cor.res.dat", "clim.dat", "rwl.dat")
  }

  return(res.list)
} # End of function
