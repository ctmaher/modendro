#' Growth-climate correlations via moving windows
#'
#' @description
#' Function to take a tree ring chronology and a monthly climate variable to compute
#' monthly aggregate correlations for every combination of consecutive months for calendar years and water years (Oct-Sept)
#' The function computes the correlations for a specified number lagged years (i.e., years before ring formation).
#' This function is designed for exploratory data analysis on growth-climate relationships in trees.
#' 3 years is the default and is a sensible maximum for most analyses.
#'
#' @details
#' This
#'
#'
#' @param chrono a chron object (such as that produced by dplR's `chron()`).
#' @param clim a `data.frame` with at least 3 columns: year, month (numeric), and a climate variable.
#' @param var character vector - the colname of the climate variable of interest in the `clim` data.frame.
#' @param chrono.col character vector - the colname of the chronology series (default is "xxxstd", which is the defualt produced by dplR's `chron()`).
#' @param agg.fun character vector specifying the function to use for aggregating monthly climate combinations.
#' Options are "mean" or "sum", e.g., for temperature or precipitation data, respectively. Default is "mean".
#' @param max.lag numeric vector specifying how many years of lag to calculate calculations for. Default is 3 years.
#' @param corr.method character vector specifying which correlation method to use. Passes to `cor.test()`.
#' @param chrono.name character vector - the name of your chronology. This is used in the title of your plot.
#' @param plots logical vector indicating whether or not to produce plots.
#'
#' @return A 1-2 element list containing data.frames of the correlation results and the default plots of the same data.
#'
#' @import plyr
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' Will add some later
#'


n_mon_corr <- function(chrono = NULL, clim = NULL,
                       var = NULL, clim.rel.per.begin = NULL,
                       chrono.col = "std", agg.fun = "mean",
                       max.lag = 2, corr.method = "pearson",
                       chrono.name = NULL, plots = TRUE){

  if (chrono.col %in% colnames(chrono)){

  } else {
    stop("Please supply the chrono.col argument with the correct name
         of the growth variable in the chronology data.frame")
  }

  stopifnot("No valid climatically relevant period begin month provided (must be a numeric month)" =
              is.numeric(clim.rel.per.begin))

  # make sure that "year" columns are labelled as such
  colnames(clim)[which((substr(colnames(clim), start = 1, stop = 2)
                        %in% c("Ye","ye")) == T)] <- "year"
  colnames(chrono)[which((substr(colnames(chrono), start = 1, stop = 2)
                          %in% c("Ye","ye")) == T)] <- "year"
  # same for month
  colnames(clim)[which((substr(colnames(clim), start = 1, stop = 2)
                        %in% c("Mo","mo")) == T)] <- "month"

  # make sure year and month are integers
  clim$year <- as.integer(clim$year)
  # I think some older versions of dplR included a year column, now it does not
  # Get the year from row names
  chrono$year <- rownames(chrono) |> as.integer()

  # Make month an integer for now
  clim$month <- as.integer(clim$month)

  # Create a chronological sequence of months starting with the numeric month
  # given as the clim.rel.per.begin argument
  mon.seq <- rep(1:12, 2) # A list of months representing 2 whole calendar years.
  mon.seq <- mon.seq[clim.rel.per.begin:length(mon.seq)][1:12] # 12 months in a row
  # This seq of months will then be relevant to a particular current growth year, though the
  # months span 2 calendar years. According to standard practice, for N hemisphere tree rings, the
  # growth year is the later calendar year. In the S hemisphere, it is the earlier calendar year.

  # This should skip the water year concept in favor of the concept of a "growth year"
  #clim$month <- factor(clim$month, levels = mon.seq)

  clim$growyear <- if (clim.rel.per.begin %in% 1:6) {
    offset <- clim.rel.per.begin - 1
    c(rep(min(clim$year) - 1, offset), clim$year[1:(length(clim$year) - offset)])
  } else {
    offset <- 12 - clim.rel.per.begin + 1
    c(clim$year[(offset + 1):length(clim$year)], rep(max(clim$year) + 1, offset))
  }

  # add water year for Oct-Dec
  #clim$water.year <- ifelse(clim$month %in% c(10,11,12), clim$year + 1, clim$year)

  # Find the complete years in the climate data
  clim.complete <- aggregate(month ~ growyear, data = clim, length)
  clim.complete <- clim.complete[clim.complete$month == 12,]
  clim <- clim[clim$growyear %in% clim.complete$growyear,]

  # Find min and max complete years for the correlations - for the plots later
  min.y <- max(min(chrono$year), min(clim$growyear))
  max.y <- min(max(chrono$year), max(clim$growyear))

  # Create a vector of all possible combinations of months
  # Regular calendar year first
  mos.mat <- expand.grid(mon.seq, mon.seq)
  mos <- apply(mos.mat, MARGIN = 1, FUN = \(x){
    seq(from = which(mon.seq %in% x[1]),
        to = which(mon.seq %in% x[2]))
  })
  # remove the non-ascending sequences
  mos <- lapply(mos, FUN = \(x) {
    if (x[1] > x[length(x)]) {
      x <- NA
    } else {
      x
    }
  })

  mos <- mos[!is.na(mos)]

  # These are indices, so now apply them to the mon.seq
  mos <- lapply(mos, FUN = \(x) {
    mon.seq[x]
  })

  # Annual lags?
  # Hold the data.frame of results in a list, with the length of the list being equal to
  # the max lag?
  lag.seq <- 0:max.lag
  #lag.list <- vector("list", length(lag.seq))
  #for (y in seq_along(lag.seq)){
  lag.list <- ldply(lag.seq, .fun = function(l){
    # Challenge: vectorize the code below. The order of month groupings doesn't really matter,
    # So I could make this a lot faster by not doing a loop.

    # calendar year
    cor.results <- ldply(mos, .fun = \(x){
      # Aggregate the variable of interest for the given month sequence
      clim.mo <- aggregate(formula(paste(var, "growyear", sep = "~")),
                           data = clim[clim$month %in% x, ],
                           FUN = \(z){ifelse(agg.fun %in% "mean", mean(z), sum(z))})
      # attach the chronology to the climate data - enter the current lag before this step
      chrono.new <- chrono
      chrono.new$year <- chrono.new$year - l
      clim.mo <- merge(clim.mo, chrono.new, by.x = "growyear", by.y = "year")
      # Run the correlation test between climate and the chronology
      ct <- cor.test(clim.mo[,var], clim.mo[, chrono.col], method = corr.method)
      # put the results together in a data.frame
      if (corr.method %in% "pearson") {

        # which.months <- if (x[1] > x[length(x)]) {
        #   paste(x[1], x[length(x)], sep = ":")
        # }
        result <- data.frame(months = ifelse(length(x) > 1,
                                             paste(x[1], x[length(x)], sep = ":"),
                                             paste(x)),
                             coef = round(ct$estimate[[1]], 3),
                             p = ct$p.value[[1]],
                             ci.lo = ct$conf.int[1],
                             ci.hi = ct$conf.int[2])
      } else {
        result <- data.frame(months = ifelse(length(x) > 1,
                                             paste(x[1], x[length(x)], sep = ":"),
                                             paste(x)),
                             coef = round(ct$estimate[[1]], 3),
                             p = ct$p.value[[1]])
      }
      # return the result
      result
    })

    cor.results$sig <- ifelse(cor.results$p > 0.05, "",
                              ifelse(cor.results$p <= 0.05 & cor.results$p >= 0.01, "*",
                                     ifelse(cor.results$p <= 0.001, "***", "**")))

    cor.results$lag <- ifelse(l == 0, paste(l), paste0("-", l))
    cor.results
  })

  lag.list$lag <- factor(lag.list$lag, levels = lag.seq*-1)
  # sort by correlation coef
  lag.list <- lag.list[order(lag.list$coef, decreasing = T),]


  if (plots == TRUE) {

    plot.df <- lag.list
    for (i in 1:nrow(plot.df)){
      plot.df$start.mo[i] <- as.numeric(strsplit(plot.df$months, ":")[[i]][1])
    }
    for (i in 1:nrow(plot.df)){
      plot.df$end.mo[i] <- as.numeric(ifelse(length(strsplit(plot.df$months, ":")[[i]]) == 2,
                                             strsplit(plot.df$months, ":")[[i]][2],
                                             strsplit(plot.df$months, ":")[[i]][1]))
    }

    plot.df$start.mo <- unlist(plot.df$start.mo) |> as.numeric()
    plot.df$end.mo <- unlist(plot.df$end.mo) |> as.numeric()

    # plot.df$start.mo <- ifelse(plot.df$seas.win %in% "Water year" & plot.df$start.mo %in% 12, 0, plot.df$start.mo)
    # plot.df$start.mo <- ifelse(plot.df$seas.win %in% "Water year" & plot.df$start.mo %in% 11, -1, plot.df$start.mo)
    # plot.df$start.mo <- ifelse(plot.df$seas.win %in% "Water year" & plot.df$start.mo %in% 10, -2, plot.df$start.mo)
    # plot.df$end.mo <- ifelse(plot.df$seas.win %in% "Water year" & plot.df$end.mo %in% 12, 0, plot.df$end.mo)
    # plot.df$end.mo <- ifelse(plot.df$seas.win %in% "Water year" & plot.df$end.mo %in% 11, -1, plot.df$end.mo)
    # plot.df$end.mo <- ifelse(plot.df$seas.win %in% "Water year" & plot.df$end.mo %in% 10, -2, plot.df$end.mo)

    plot.df$simp.sig <- ifelse(plot.df$sig %in% c("*","**","***"), "Sig.","Not sig.")

    title <- ifelse(is.null(chrono.name), paste0("Chronology correlations with ",var),
                    paste0(chrono.name, " chronology correlations with ",var))

    out.plot <- ggplot(plot.df, aes(factor(start.mo, levels = mon.seq), coef, color = simp.sig)) +
      geom_point(shape = 5, size = 1.5) +
      geom_point(data = plot.df, aes(factor(end.mo, levels = mon.seq), coef, color = simp.sig), shape = 124, size = 2) +
      scale_x_discrete(breaks = mon.seq, labels = mon.seq) +
      scale_y_continuous(breaks = seq(-1,1, by = 0.1)) +
      geom_segment(aes(xend = end.mo, yend = coef)) + xlab("Month") +
      ylab(paste0("Correlation coefficient\n(",corr.method,")")) +
      scale_color_manual(name = "", values = c("grey80","black")) +
      theme_bw() + facet_wrap(~lag, ncol = 1) +
      ggtitle(label = title,
              subtitle = paste0("Monthly climate ",
                                agg.fun, "s with annual lag 0:",
                                max.lag," (", min.y,"-",max.y,")"))

    out.list <- list(lag.list, out.plot)
    return(out.list)
  } else {
    return(lag.list)
  }

} # End of function
