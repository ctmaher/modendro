#' Calculate multiple growth-climate relationships
#'
#' @description
#' This function is a workhorse function behind \code{\link{n_mon_corr}}. It splits and merges the
#' climate data and the rwl data, then calculates all the correlations between them.
#'
#' Internal function for modendro.


multi_clim_gro_corr <- function(rwl.group = NULL,
                                clim.group = NULL,
                                clim.var = NULL,
                                group.IDs.df = NULL,
                                group.var = NULL,
                                gro.period.end = NULL,
                                agg.fun = "mean",
                                max.lag = 1,
                                prewhiten = TRUE,
                                corr.method = "spearman",
                                hemisphere = NULL) {
  ## Compute the moving windows
  clim1 <- moving_win_multi(
    clim.group,
    clim.var = clim.var,
    win.lens = 2:12,
    agg.fun = agg.fun
  )

  # Split by month
  clim1.mo.split <- split(clim1, f = clim1$start.month)


  # Now is a good time to get the arima residuals - before we create copies of the data when we
  # set up the lagged datasets.
  if (prewhiten == TRUE) {
    # Climate data first
    clim1.mo.split <- lapply(clim1.mo.split, FUN = \(start.mo) {
      # Split into the different win.lens
      start.mo.split <- split(start.mo, f = start.mo$win.len)
      lapply(start.mo.split, FUN = \(this.win) {
        # Make sure we have years in order!
        this.win <- this.win[order(this.win[, "year"], decreasing = FALSE), ]
        this.win[, clim.var] <- forecast::auto.arima(this.win[, clim.var], seasonal = FALSE) |>
          residuals() |>
          as.numeric()
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

    clim.lag.bind <- lapply(clim1.mo.split, FUN = \(start.mo) {
      lag.years <- sapply(-1:max.lag, FUN = \(n.lag) {
        this.lag <- data.frame(start.mo[, "year"] + n.lag)
        colnames(this.lag) <- paste0("lag", n.lag, ".year")
        this.lag
      }) |> do.call(what = "cbind")
      cbind(start.mo, lag.years)
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

    clim.lag.bind <- lapply(clim1.mo.split, FUN = \(start.mo) {
      lag.years <- sapply(0:max.lag, FUN = \(n.lag) {
        this.lag <- data.frame(start.mo[, "year"] + n.lag)
        colnames(this.lag) <- paste0("lag", n.lag, ".year")
        this.lag
      }) |> do.call(what = "cbind")
      cbind(start.mo, lag.years)
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
      lag.x.mo.split <- split(lag.x, f = lag.x[, "start.month"])

      # Take each month split and do all the possible correlations
      this.lag.corr <- lapply(lag.x.mo.split, FUN = \(mo.split) {
        # Take each possible correlation out of mo.split and run the correlation
        apply(all.corr.combos, MARGIN = 1, FUN = \(combo) {
          series.name <- combo[["series.name"]]
          win.len <- combo[["win.len"]]

          # The correlations
          # This subset eliminates NA values in either series.
          # NAs in the middle of the series should be a rare occurrence
          subs.for.corr <- mo.split[mo.split[, "win.len"] == win.len, c("year",
                                                                        "start.month",
                                                                        "win.len",
                                                                        as.character(series.name),
                                                                        clim.var)] |>
            na.omit()

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
            data.frame(
              series = as.character(series.name),
              start.month = unique(subs.for.corr[, "start.month"]),
              win.len = unique(subs.for.corr[, "win.len"]),
              coef = cor.res[["estimate"]][[1]],
              p = cor.res[["p.value"]][[1]],
              dir = ifelse(cor.res[["estimate"]][[1]] < 0, "Neg.", "Pos."),
              overlap = min(length(na.omit(
                mo.split[, as.character(series.name)]
              )), length(na.omit(
                mo.split[mo.split[, "win.len"] == win.len, clim.var]
              ))),
              corr.method = corr.method
            )
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
            data.frame(
              series = as.character(series.name),
              start.month = unique(mo.split[, "start.month"]),
              win.len = unique(subs.for.corr[, "win.len"]),
              coef = cor.res[["rho"]],
              p = cor.res[["pval"]],
              dir = ifelse(cor.res[["rho"]] < 0, "Neg.", "Pos."),
              overlap = min(length(na.omit(
                mo.split[, as.character(series.name)]
              )), length(na.omit(
                mo.split[mo.split[, "win.len"] == win.len, clim.var]
              ))),
              corr.method = corr.method
            )
          }

        }) |> do.call(what = "rbind")
      }) |> do.call(what = "rbind")

      this.lag.corr$lag <- ifelse(substr(unique(which.lag), 5, 5) == "0", "0",
                                  substr(unique(which.lag), 4, 5))
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
