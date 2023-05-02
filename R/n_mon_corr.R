#' Growth-climate correlations via moving windows
#'
#' @description
#' Exploratory Data Analysis function to take a tree ring chronology and a monthly climate variable to compute
#' static monthly aggregate correlations for every combination of consecutive months for calendar years and "water years"
#' that are putatively relevant climatically to each year of tree ring formation (e.g., Oct-Sept for N hemisphere, Apr-Mar for S hemisphere).
#' The function computes the correlations for a specified number lagged years (i.e., years before ring formation).
#' 2 years is the default and is a sensible maximum for most analyses. The correlations are "static" because there
#' is no moving window analysis - the correlations are for the entire period of maximum overlap between your chronology
#' and climate time series data.
#'
#' Because this is EDA, the "climatically relevant period" is defined as a an inclusive 12 month period beginning
#' on the month you specify with `clim.rel.per.begin`. There is no sense in trying to shorten this period because
#' all combinations of consecutive months are run by default - this includes all the possible shorter periods.
#' Longer periods are essentially covered by the lags, so periods longer than 12 months will be irrelevant for
#' the majority of cases (& thus are not possible to specify). Chose for `clim.rel.per.begin` the first month in
#' after which you are reasonably certain that no tree ring growth is occurring at your site. Note that this analysis
#' assumes that there is a seasonality of growth.
#'
#' Fair warning: this is a basic function that will accept any tree ring chronology and climate data in the proper format.
#' It is the user's responsibility to make sure that the chronology is properly constructed -
#' properly dealing with the "dark arts" of detrending & standardization and also accounting for temporal autocorrelation in the data.
#' If you don't know what this means or you just took a chronology directly from the ITRDB without knowing how it was made,
#'  you have some homework to do.
#'
#'
#' @param chrono a `chron` object (such as that produced by dplR's `chron()`).
#' @param clim a `data.frame` with at least 3 columns: year, month (numeric), and a climate variable.
#' @param var character vector - the colname of the climate variable of interest in the `clim` data.frame.
#' @param clim.rel.per.begin an integer month representing the beginning of the climatically relevant period to the growth year (always a 12 month period).
#' This will include the "water year" of the calendar year before growth. E.g., 10 for N hemisphere (the default), 4 for S hemisphere. See details below for more info.
#' @param chrono.col character vector - the colname of the chronology series (default is "std", which is the defualt produced by dplR's `chron()`).
#' @param agg.fun character vector specifying the function to use for aggregating monthly climate combinations.
#' Options are "mean" or "sum", e.g., for temperature or precipitation data, respectively. Default is "mean".
#' @param max.lag numeric vector specifying how many years of lag to calculate calculations for. Default is 2 years.
#' @param corr.method character vector specifying which correlation method to use. Passes to `cor.test()`.
#' @param chrono.name character vector - the name of your chronology (optional). This is used in the title of your plot.
#' If you produce many plots, this helps keep them identifiable.
#' @param plots logical vector indicating whether or not to produce plots. Default is TRUE.
#'
#' @details
#' Exploring a wide range of plausible growth-climate relationships can be a useful first step once you have
#' a collection of cross-dated tree ring series and have properly detrended them, standardized them, and dealt
#' with temporal autocorrelation. There aren't universal answers to these tasks - each dataset is unique.
#'
#' A note on choosing tree ring analyses based in the Southern hemisphere:
#' `n_mon_corr()` is designed to work in both the Northern and Southern hemispheres. Hemisphere matters
#' for tree ring growth-climate relationships because tree ring formation in the Southern hemisphere typically
#' spans two calendar years (e.g., starting in Nov 2000 and ending in Mar of 2001).
#' It was Schulman's (1956) protocol to assign the earlier calendar year to the tree rings in the Southern hemisphere,
#' i.e., the calendar year in which growth began. `n_mon_corr()` assumes your data follows this standard as well.
#' This has implications for how the climate data is aligned with the treering data.
#' The current implementation handles this implicitly by assuming that if `clim.rel.per.begin` is between 1:6,
#' this is a S. hemisphere analysis and the current "growth year" is the same as the calendar year of `clim.rel.per.begin`.
#' If `clim.rel.per.begin` is between 7:12, it is assumed that the is a N. hemisphere analysis and the current "growth year"
#' is the calendar year following current `clim.rel.per.begin`. E.g., if `clim.rel.per.begin = 4`, the climatically relevant period
#' will be defined as months `c(4,5,6,7,8,9,10,11,12,1,2,3)` with the calendar year of the FIRST 9 months as the
#' "growth year". If `clim.rel.per.begin = 10`, the climatically relevant period
#' will be defined as months `c(10,11,12,1,2,3,4,5,6,7,8,9)` with the calendar year of the LAST 9 months as the
#' "growth year".
#'
#' Interpreting the plots:
#'
#'
#'
#' @return A 1-2 element list containing data.frames of the correlation results and the default plots of the same data.
#'
#' @references
#' Schulman, E. (1956) \emph{Dendroclimatic changes in semiarid America}, University of Arizona Press.
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
                       var = NULL, clim.rel.per.begin = 10,
                       chrono.col = "std", agg.fun = "mean",
                       max.lag = 2, corr.method = "pearson",
                       chrono.name = NULL, plots = TRUE){


  # Error catching
  stopifnot("Arg chrono or clim are not an object of class 'chron', 'data.frame', or 'matrix'" =
              data.class(chrono) %in% "chron" |
              data.class(chrono) %in% "data.frame" |
              data.class(chrono) %in% "matrix" |
              data.class(clim) %in% "data.frame" |
              data.class(clim) %in% "matrix"
  )

  match.test <- var %in% colnames(clim)
  stopifnot("Arg var must match one unique column name in clim" =
              length(match.test[match.test == TRUE]) == 1
  )

  match.test <- chrono.col %in% colnames(chrono)
  stopifnot("Arg chrono.col must match the name
         of the growth variable in the chrono data.frame" =
              length(match.test[match.test == TRUE]) == 1
  )

  stopifnot("Arg agg.fun must be either 'mean' or 'sum'" =
              agg.fun %in% "mean" |
              agg.fun %in% "sum"
  )


  stopifnot("No valid climatically relevant period begin month provided (must be an integer month)" =
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
  clim$month <- as.integer(clim$month)

  # Get the year from row names
  chrono$year <- rownames(chrono) |> as.integer()

  # Create a chronological sequence of months starting with the numeric month
  # given as the clim.rel.per.begin argument
  mon.seq <- rep(1:12, 2) # A list of months representing 2 whole calendar years.
  mon.seq <- mon.seq[clim.rel.per.begin:length(mon.seq)][1:12] # 12 months in a row
  # This seq of months will then be relevant to a particular current growth year, though the
  # months span 2 calendar years. According to standard practice, for N hemisphere tree rings, the
  # growth year is the later calendar year. In the S hemisphere, it is the earlier calendar year.

  # Define the "growth year" based on the clim.rel.per.begin

  clim$growyear <- if (clim.rel.per.begin %in% 1:6) {
    offset <- clim.rel.per.begin - 1
    c(rep(min(clim$year) - 1, offset), clim$year[1:(length(clim$year) - offset)])
  } else {
    offset <- 12 - clim.rel.per.begin + 1
    c(clim$year[(offset + 1):length(clim$year)], rep(max(clim$year) + 1, offset))
  }


  # Find the complete years in the climate data
  clim.complete <- aggregate(month ~ growyear, data = clim, length)
  clim.complete <- clim.complete[clim.complete$month == 12,]
  clim <- clim[clim$growyear %in% clim.complete$growyear,]

  # Find min and max complete years for the correlations, i.e., the complete overlap
  min.y <- max(min(chrono$year), min(clim$growyear))
  max.y <- min(max(chrono$year), max(clim$growyear))

  # Create a vector of all possible combinations of months
  # Regular calendar year first
  mos.mat <- expand.grid(mon.seq, mon.seq)
  mos <- apply(mos.mat, MARGIN = 1, FUN = \(x){
    seq(from = which(mon.seq %in% x[1]),
        to = which(mon.seq %in% x[2]))
  })
  # remove the non-ascending sequences - these are not sensical in this analysis
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

  # Annual lags
  # Hold the data.frame of results in a list, with the length of the list being equal to
  # 1 + the max lag
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

  # Make the lag a factor
  lag.list$lag <- factor(lag.list$lag, levels = lag.seq*-1)
  # sort by correlation coef
  lag.list <- lag.list[order(lag.list$coef, decreasing = T),]

  # Set up a data.frame for the plots - add some variables that are useful for plotting, but not
  # for the main output.
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

    # Unlist these and make sure they are numeric.
    plot.df$start.mo <- unlist(plot.df$start.mo) |> as.numeric()
    plot.df$end.mo <- unlist(plot.df$end.mo) |> as.numeric()

    # Simplified version of significance
    plot.df$simp.sig <- ifelse(plot.df$sig %in% c("*","**","***"), "Sig.","Not sig.")

    # Build a nice title for the plots
    title <- ifelse(is.null(chrono.name), paste0("Chronology correlations with ",var),
                    paste0(chrono.name, " chronology correlations with ",var))

    # Make the plot
    out.plot <- ggplot(plot.df) +
      geom_point(aes(factor(start.mo, levels = mon.seq), # start points
                     coef, color = simp.sig),
                 shape = 18, size = 2) +
      geom_point(aes(factor(end.mo, levels = mon.seq), # end points
                     coef, color = simp.sig),
                 shape = 124, size = 2) +
      geom_segment(aes(x = factor(start.mo, levels = mon.seq), # lines connecting
                       xend = factor(end.mo, levels = mon.seq),
                       y = coef, yend = coef, color = simp.sig)) +
      scale_x_discrete(breaks = mon.seq, labels = mon.seq) +
      scale_y_continuous(breaks = seq(-1,1, by = 0.1)) +
      xlab("Month") +
      ylab(paste0("Correlation coefficient\n(", corr.method, ")")) +
      scale_color_manual(name = "", values = c("grey80","black")) +
      theme_bw() +
      facet_wrap(~lag, ncol = 1, strip.position	= "right") +
      ggtitle(label = title,
              subtitle = paste0("Monthly climate ",
                                agg.fun, "s with annual lag 0:",
                                max.lag, " (years ", min.y, "-", max.y, ")"))

    out.list <- list(lag.list, out.plot)
    return(out.list)
  } else {
    return(lag.list)
  }

} # End of function
