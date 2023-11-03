#' Flexible monthly aggregate growth-climate cross correlations for exploratory data analysis
#'
#' @description
#' Exploratory data analysis (EDA) function to compute correlations between a tree ring chronology and a monthly climate variable
#' aggregated for every combination (lengths 1:12) of consecutive months inside of a 12-month long "relevant climate period" that could conceivably be relevant
#' to growth in any given year.
#'
#' The ideas behind the relevant climate period are that climate no longer has an effect on radial growth after growth has stopped &
#' that current year's growth could have been influenced by climate in any month AFTER the previous year's growth stopped.
#'
#' The user specifies the beginning month of the relevant climate period - this is the 1st month after radial growth stops
#' (i.e., the 1st month of the fall season). For example, if radial growth typically terminates sometime in September, the user would enter
#' `rel.per.begin = 10` to specify a relevant climatic period that starts in October of the previous year and ends in September of
#' the next year (i.e., the calendar year of growth). In this way the relevant climatic period covers the "water year" months leading
#' up to the growing period and only extends through the months when growth occurs. It is nonsensical to include months after
#' radial growth stops - doing so may result in spurious correlations.
#'
#' Compared to methods that force rigidly-defined seasons of a fixed length, this approach should allow greater discovery of meaningful growth-climate relationships.
#'
#' Fair warning: this is a basic function that will accept any tree ring chronology and climate data in the proper format.
#' It is the user's responsibility to make sure that the chronology is properly constructed -
#' properly dealing with the "dark arts" of detrending & standardization (see \code{\link[dplR]{detrend.series}}) and also accounting for temporal autocorrelation in the data.
#' If you don't know what this means or you just took a chronology directly from the ITRDB without knowing how it was made,
#'  you have some homework to do. There aren't universal answers to these tasks - each dataset is somewhat unique.
#'
#' @param chrono a `chron` object (such as that produced by dplR's \code{\link[MASS]{chron}}). Make sure this has a `year` variable.
#' @param chrono.col character vector - the colname of the chronology series (default is "std", which is the defualt produced by dplR's \code{\link[MASS]{chron}}).
#' @param clim a `data.frame` with at least 3 columns: year, month (numeric), and a climate variable.
#' @param clim.var character vector - the colname of the climate variable of interest in the `clim` data.frame.
#' @param rel.per.begin an integer month representing the beginning of the climatically relevant period to the growth year (always a 12 month period).
#' This will include the "water year" of the calendar year before growth. E.g., 10 for N hemisphere, 4 for S hemisphere. See details below for more info.
#' @param hemisphere a character vector specifying which hemisphere your chronology - & climate data - comes from ("N" or "S").
#' Conventions for assigning growth years - and thus aligning tree ring and climate data - are different for N and S hemisphere.
#' @param agg.fun character vector specifying the function to use for aggregating monthly climate combinations.
#' Options are "mean" or "sum", e.g., for temperature or precipitation data, respectively. Default is "mean".
#' @param max.lag numeric vector specifying how many years of lag to calculate calculations for. Default is 1 year.
#' @param prewhiten logical vector specifying whether or not to convert climate time series to AR residuals (aka "prewhitening").
#' This removes autocorrelation in a time series, leaving only the high-frequency variation. If you are doing this, you should also construct your chronology from AR residuals/prewhitened series. Default is FALSE.
#' @param auto.corr logical vector specifying whether there is temporal autocorrelation in either your tree ring chronology or climate time series (there typically is autocorrelation, unless both are "prewhitened").
#' If TRUE (the default), & corr.method is "spearman" or "kendall", then the \code{\link[corTESTsrd]{corTESTsrd}}function is used to compute modified significance testing to account for autocorrelation (From Lun et al. 2022).
#' Caution! Currently auto.corr = TRUE & corr.method = "Pearson" doesn't make any adjustments. This will be included soon.
#' @param corr.method character vector specifying which correlation method to use. Options are `c("pearson", "kendall", "spearman")`.
#'  Passes to \code{\link[stats]{cor.test}} or to \code{\link[corTESTsrd]{corTESTsrd}}.
#' @param chrono.name character vector - the name of your chronology (optional). This is used in the title of your plot.
#' If you produce many plots, this helps keep them identifiable.
#' @param plots logical vector indicating whether or not to produce plots. Default is TRUE.
#' @param silent logical vector indicating whether messages about relevant period and hemisphere conventions will be printed. Default is FALSE.
#'
#' @details
#' Exploring a wide range of plausible growth-climate relationships can be a useful first step once you have
#' a collection of cross-dated tree ring series and have properly detrended them, standardized them, and dealt
#' with temporal autocorrelation.
#'
#' The default correlation test method is Pearson product moment correlation coefficient, although this might
#' not be appropriate for your analysis.
#'
#' A note on tree ring analyses based in the Southern hemisphere:
#' \code{\link{n_mon_corr}} is designed to work in both the Northern and Southern hemispheres. Hemisphere matters
#' for tree ring growth-climate relationships because tree ring formation in the Southern hemisphere typically
#' spans two calendar years (e.g., starting in Nov 2000 and ending in Mar of 2001).
#' It was Schulman's (1956) protocol to assign the earlier calendar year to the tree rings in the Southern hemisphere,
#' i.e., the calendar year in which growth began. \code{\link{n_mon_corr}} assumes your data follows this standard as well.
#' This has implications for how the climate data is aligned with the treering data.
#' The current implementation handles this implicitly by assuming that if `rel.per.begin` is between 1:6,
#' this is a S. hemisphere analysis and the current "growth year" is the same as the calendar year of `rel.per.begin`.
#' If `rel.per.begin` is between 7:12, it is assumed that the is a N. hemisphere analysis and the current "growth year"
#' is the calendar year following `rel.per.begin`. E.g., if `rel.per.begin = 4`, the climatically relevant period
#' will be defined as months `c(4,5,6,7,8,9,10,11,12,1,2,3)` with the calendar year of the FIRST 9 months as the
#' "growth year". If `rel.per.begin = 10`, the climatically relevant period
#' will be defined as months `c(10,11,12,1,2,3,4,5,6,7,8,9)` with the calendar year of the LAST 9 months as the
#' "growth year".
#'
#' Interpreting the plots:
#' The plots show a 12-month sequence of consecutive months on the x-axis & the correlation coefficient on the y-axis.
#' The diamonds indicate the starting month of an n-month aggregate period, small vertical bars the end. Horizontal lines connect
#' the start and end months for periods > 1 month. Significant correlations (as determined by \code{\link[stats]{cor.test}}) are shown
#' in black, no significant ones in grey. Plot panel labels (right-hand side of plots) indicate lag years: 0 = current year,
#' -1 = previous year, -2 = 2 years back.
#'
#' @return A 2-3 element list containing data.frames of the correlation results, the data used in the correlations, and the default plots of the same data.
#'
#' @references
#' Schulman, E. (1956) \emph{Dendroclimatic changes in semiarid America}, University of Arizona Press.
#'
#' Lun, D., S. Fischer, A. Viglione, and G. Blöschl. (2022). Significance testing of rank cross-correlations between autocorrelated time series with short-range dependence, \emph{Journal of Applied Statistics}:1–17.
#'
#' @import ggplot2
#' @import corTESTsrd
#'
#' @export
#'
#' @examples
#' Will add some later
#'


n_mon_corr <- function(chrono = NULL,
                       chrono.col = "std",
                       clim = NULL,
                       clim.var = NULL,
                       rel.per.begin = NULL,
                       hemisphere = NULL,
                       agg.fun = "mean",
                       max.lag = 1,
                       prewhiten = FALSE,
                       auto.corr = FALSE,
                       corr.method = "pearson",
                       chrono.name = NULL,
                       plots = TRUE,
                       silent = FALSE){


  ## Initial error catching and interactive prompts

  stopifnot("Arg chrono or clim are not an object of class 'chron', 'data.frame', or 'matrix'" =
              data.class(chrono) %in% "chron" |
              data.class(chrono) %in% "data.frame" |
              data.class(chrono) %in% "matrix" |
              data.class(clim) %in% "data.frame" |
              data.class(clim) %in% "matrix"
  )

  # stopifnot("No year column in chronology dataframe, please add one (year may be contained in rownames)" =
  #             substr(colnames(chrono), 1, 1) %in% c("Y","y")
  #           )

  match.test <- clim.var %in% colnames(clim)
  stopifnot("Arg clim.var must match one unique column name in clim" =
              length(match.test[match.test == TRUE]) == 1
  )

  match.test <- colnames(chrono) %in% chrono.col
  stopifnot("Arg chrono.col must match the name
         of the growth variable in the chrono data.frame" =
              length(match.test[match.test == TRUE]) == 1
  )

  if (is.null(rel.per.begin)) {
    cat("You haven't specified the beginning month of the relevant climate period -\n",
        "this is the 1st month after radial growth typically stops in a year\n",
        "(i.e., the 1st month of the fall season at your site).\n",
        "This is approximate, but you should have a reasonable idea of what month this is for your study system.\n")
    rel.per.begin <- readline(prompt = "First month of relevant climate period = ") |> as.integer()
  }

  if (is.null(hemisphere)) {
    cat("You haven't specified the hemisphere from which your tree ring series comes from.\n",
        "This is important because there are different conventions for linking growth years\n",
        "to climate years for N vs. S hemispheres")
    hemisphere <- readline(prompt = "Enter hemisphere ('N' or 'S') = ")
  }

  stopifnot("Invalid climatically relevant period begin month provided (must be a single integer month)" =
              is.numeric(rel.per.begin) &
              length(rel.per.begin) == 1)

  stopifnot("Invalid hemisphere argument provided (must be a character vector & either 'S' or 'N')" =
              is.character(hemisphere) &
              substr(hemisphere, 1, 1) %in% c("s","S","N","n")) # actually more permissive than the error message suggests

  stopifnot("Arg agg.fun must be either 'mean' or 'sum'" =
              agg.fun %in% "mean" |
              agg.fun %in% "sum"
  )

  stopifnot("Arg max.lag must be a numeric vector of length = 1" =
              length(max.lag) == 1 |
              is.numeric(max.lag)
  )

  # accept max.lag inputs that have a negative in front
  if (max.lag < 0) {
    max.lag <- as.numeric(max.lag) |> abs()
  }

  stopifnot("Arg corr.method must be an exact match of one of these: c('pearson','kendall','spearman')" =
              corr.method %in% c("pearson", "kendall", "spearman")
  )


  # make sure that "year" columns are labelled as such
  colnames(clim)[which((substr(colnames(clim), start = 1, stop = 1)
                        %in% c("Y","y")) == T)] <- "year"

  # same for month
  colnames(clim)[which((substr(colnames(clim), start = 1, stop = 1)
                        %in% c("M","m")) == T)] <- "month"

  stopifnot("Month variable in climate data not numeric or integer" =
              is.numeric(clim[,"month"]) |
              is.integer(clim[,"month"]) |
              is.double(clim[,"month"])
  )

  # make sure year and month are integers from here on out
  clim[,"year"] <- as.integer(clim[,"year"])
  clim[,"month"] <- as.integer(clim[,"month"])

  # Get the year from row names
  if (any(substr(colnames(chrono), 1, 1) %in% c("Y","y")) == FALSE) {
    chrono[,"year"] <- rownames(chrono) |> as.numeric()
  } else {
    colnames(chrono)[which((substr(colnames(chrono), start = 1, stop = 1)
                            %in% c("Y","y")) == T)] <- "year"
    chrono[,"year"] <- as.numeric(chrono[,"year"])
  }

  # n_mon_corr assumes that all years have all 12 months! If even one month is missing somewhere, this will
  # mess up everything that follows.

  mon.count <- aggregate(month ~ year, data = clim, length)

  # stopifnot("Not all years in climate data have all 12 months represented. n_mon_corr() requires that all years have all 12 months." =
  #             all(mon.count$month == 12)
  # )
  if (all(mon.count$month != 12)) {

    stop(paste("Year", mon.count$year[mon.count$month < 12],
                   "does not have all 12 months represented."))
  }

  # n_mon_corr also assumes absolute regularity (this is true for some of the correlation tests too) in both chrono & clim
  chron.year.seq <- chrono[,"year"]
  chron.year.seq.diff <- chron.year.seq[order(chron.year.seq)] |> diff()
  clim.year.seq <- unique(clim[,"year"])
  clim.year.seq.diff <- clim.year.seq[order(clim.year.seq)] |> diff()
  if (any(chron.year.seq.diff != 1)) {
    paste("Year", chron.year.seq[which(chron.year.seq.diff > 1)], "is missing from chronology.")
  }
  stopifnot("Chronology does not have complete continuous years." =
              all(chron.year.seq.diff == 1)
  )

  if (any(clim.year.seq.diff != 1)) {
    paste("Year", clim.year.seq[which(clim.year.seq.diff > 1)], "is missing from climate data.")
  }
  stopifnot("Climate data does not have complete continuous years." =
              all(clim.year.seq.diff == 1)
  )


  # Give a warning & maybe stop the function if there is autocorrelation in the tree ring series
  # There should be a prompt (verbal or otherwise)
  ac.test <- ar(x = na.omit(chrono[order(chrono[,"year"]), chrono.col]))
  if (ac.test$order > 0 & auto.corr == FALSE & prewhiten == FALSE) {
    cat("Autocorrelation detected in your chronology, recommend choose auto.corr = TRUE &
        corr.method = c('spearman', 'kendall') to avoid spurious correlation results.
        You can also construct a chronology of AR residuals (aka prewhitening).\n")
    # auto.corr <- readline(prompt = "Enter auto.corr (TRUE or FALSE) = ")
    # corr.method <- readline(prompt = "Enter corr.method ('spearman' or 'kendall') = ")
  }

  # Clean up the hemisphere argument if needed
  hemisphere <- ifelse(substr(hemisphere, 1, 1) %in% c("n","N") , "N", "S")

  # Create a chronological sequence of months starting with the numeric month
  # given as the clim.rel.per.begin argument
  mon.seq <- rep(1:12, 2) # A list of months representing 2 whole calendar years.
  mon.seq <- mon.seq[rel.per.begin:length(mon.seq)][1:12] # 12 months in a row
  # This seq of months will then be relevant to a particular current growth year, though the
  # months span 2 calendar years. According to standard practice, for N hemisphere tree rings, the
  # growth year is the later calendar year. In the S hemisphere, it is the earlier calendar year.

  # Print the relevant climate period
  if (silent == FALSE) {
    cat("You have specified the following months for your relevant climate period in the\n",
        hemisphere, "hemisphere:", mon.seq, "\n")
  }

  # The operations below assume that the climate data is arranged by month, then year.
  # Let's ensure this is the case.
  clim <- clim[order(clim$year, clim$month),]

  # Define the "growth year" based on the hemisphere argument
  if (hemisphere %in% "S") {
    if (silent == FALSE) {
      message("\nAssuming Southern hemisphere conventions for linking growth years
    and climate years (see ?n_mon_corr for details)\n")
    }
    offset <- rel.per.begin - 1
    clim$growyear <- c(rep(min(clim[,"year"]) - 1, offset), clim[,"year"][1:(length(clim[,"year"]) - offset)])

  } else {
    if (silent == FALSE) {
      message("\nAssuming Northern hemisphere conventions for linking growth years
    and climate years (see ?n_mon_corr for details)\n")
    }
    offset <- 12 - rel.per.begin + 1
    clim$growyear <- c(clim[,"year"][(offset + 1):length(clim[,"year"])], rep(max(clim[,"year"]) + 1, offset))

  }

  # Give warnings - and stop the function - if someone uses unusual values for rel.per.begin for a given hemisphere
  # This is to provide guardrails for users who don't understand the rel.per & to minimize the chance that they will
  # be looking at potentially spurious correlations for months AFTER radial growth has ceased in a given year.
  # This can be done with a combination of

  # Find the complete years in the climate data
  clim.complete <- aggregate(month ~ growyear, data = clim, length)
  clim.complete <- clim.complete[clim.complete$month == 12,]
  clim <- clim[clim[,"growyear"] %in% clim.complete[,"growyear"],]

  # Find min and max complete years for the correlations, i.e., the complete overlap
  min.y <- max(min(chrono[,"year"]), min(clim[,"growyear"]))
  max.y <- min(max(chrono[,"year"]), max(clim[,"growyear"]))

  # Create a vector of all possible combinations of months
  # Regular calendar year first
  mos.mat <- expand.grid(mon.seq, mon.seq)

  mos <- apply(mos.mat, MARGIN = 1, FUN = \(x){
    seq(from = which(mon.seq %in% x[1]),
        to = which(mon.seq %in% x[2]))
  })
  # remove the non-ascending sequences - these are nonsensical in this analysis
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

  # Get the 1st month of each
  mon1 <- lapply(mos, FUN = \(x) {
    x[1]
  }) |> unlist()

  mos <- mos[order(lengths(mos), decreasing = TRUE)]

  mos.fac <- lapply(mos, FUN = \(x) {
    ifelse(length(x) > 1,
           paste(x[1], x[length(x)], sep = ":"),
           paste(x))
  }) |> unlist()

  # Annual lags
  # Hold the data.frame of results in a list, with the length of the list being equal to
  # 1 + the max lag
  lag.seq <- 0:max.lag

  lag.list <- lapply(lag.seq, FUN = function(l){

    # Run through all the month aggregates
    cor.res.df <- lapply(mos, FUN = \(x){
      # Aggregate the variable of interest for the given month sequence
      clim.mo <- aggregate(formula(paste(clim.var, "growyear", sep = "~")),
                           data = clim[clim$month %in% x, ],
                           FUN = \(z){ifelse(agg.fun %in% "mean", mean(z), sum(z))})

      # If we want to convert climate & chrono to AR residuals (aka, "prewhitening"), do it here
      # This follows what dplR's detrend.series(method = "Ar) does for tree ring series
      if (prewhiten == TRUE) {
        ar.mod <- ar(clim.mo[!is.na(clim.mo[,clim.var]),clim.var])
        ar.mean.resid <- ar.mod$resid + ar.mod$x.mean
        clim.mo[!is.na(clim.mo[,clim.var]),clim.var] <- ar.mean.resid/mean(ar.mean.resid, na.rm = TRUE)
        ar.mod2 <- ar(chrono[!is.na(chrono[,chrono.col]),chrono.col])
        ar.mean2.resid <- ar.mod2$resid + ar.mod2$x.mean
        chrono[!is.na(chrono[,chrono.col]),chrono.col] <- ar.mean2.resid/mean(ar.mean2.resid, na.rm = TRUE)
      }

      # The months vector in the desired order.
      month.vec <- ifelse(length(x) > 1,
                          paste(x[1], x[length(x)], sep = ":"),
                          paste(x))

      # attach the chronology to the climate data
      clim.mo.new <- clim.mo
      clim.mo.new$growyear <- clim.mo.new$growyear + l
      clim.mo.new <- merge(clim.mo.new, chrono, by.x = "growyear", by.y = "year")
      clim.mo.new$lag <- ifelse(l == 0, paste(l), paste0("-", l))
      clim.mo.new$months <- month.vec

      # Remove any ties from the data
      clim.mo.new <- clim.mo.new[which(!duplicated(clim.mo.new[,clim.var])),]
      clim.mo.new <- clim.mo.new[which(!duplicated(clim.mo.new[,chrono.col])),]

      # Run the correlation test between climate and the chronology
      if (auto.corr == TRUE) {
        if (corr.method %in% "pearson") {
          ct <- cor.test(clim.mo.new[, clim.var], clim.mo.new[, chrono.col],
                         method = corr.method, alternative = "two.sided")
          # put the results together in a data.frame
          result <- data.frame(months = month.vec,
                               coef = ct$estimate[[1]],
                               p = ct$p.value[[1]],
                               ci.lo = ct$conf.int[1],
                               ci.hi = ct$conf.int[2])
        } else { # if spearman or kendall
          ct <- corTESTsrd(clim.mo.new[,clim.var], clim.mo.new[, chrono.col],
                           method = corr.method,
                           iid = FALSE, alternative = "two.sided")
          # put the results together in a data.frame
          result <- data.frame(months = month.vec,
                               coef = ct[["rho"]],
                               p = ct[["pval"]])
        }
      } else {
        ct <- cor.test(clim.mo.new[,clim.var], clim.mo.new[, chrono.col], method = corr.method)

        # put the results together in a data.frame
        if (corr.method %in% "pearson") {

          result <- data.frame(months = month.vec,
                               coef = ct$estimate[[1]],
                               p = ct$p.value[[1]],
                               ci.lo = ct$conf.int[1],
                               ci.hi = ct$conf.int[2])
        } else {
          result <- data.frame(months = month.vec,
                               coef = ct$estimate[[1]],
                               p = ct$p.value[[1]])
        }
      }
      # return the correlation results and the data.frame of the merged climate and chronology data
      list(result, clim.mo.new)
    })

    cor.results <- lapply(cor.res.df, FUN = \(x) {
      x[[1]]
    }) |> do.call(what = "rbind")

    cor.df <- lapply(cor.res.df, FUN = \(x) {
      x[[2]]
    }) |> do.call(what = "rbind")

    cor.results$sig <- ifelse(cor.results$p > 0.05, "",
                              ifelse(cor.results$p <= 0.05 & cor.results$p >= 0.01, "*",
                                     ifelse(cor.results$p <= 0.001, "***", "**")))

    cor.results$lag <- ifelse(l == 0, paste(l), paste0("-", l))
    cor.df$lag <- ifelse(l == 0, paste(l), paste0("-", l))

    # Return all the results
    list(cor.results, cor.df)
  })

  lag.res <- lapply(lag.list, FUN = \(x) {
    x[[1]]
  }) |> do.call(what = "rbind")

  lag.df <- lapply(lag.list, FUN = \(x) {
    x[[2]]
  }) |> do.call(what = "rbind")

  # Make the lag a factor in both data.frames
  lag.res$lag <- factor(lag.res$lag, levels = lag.seq*-1)
  lag.df$lag <- factor(lag.df$lag, levels = lag.seq*-1)

  # sort correlation results by correlation coef
  lag.res <- lag.res[order(lag.res$coef, decreasing = T),]

  # order the months as a factor in both data.frames
  lag.res$months <- factor(lag.res$months, levels = mos.fac)
  lag.df$months <- factor(lag.df$months, levels = mos.fac)

  # Set up a data.frame for the plots - add some variables that are useful for plotting, but not
  # for the main output.
  if (plots == TRUE) {

    plot.df <- lag.res
    for (i in 1:nrow(plot.df)){
      plot.df$start.mo[i] <- as.numeric(strsplit(as.character(plot.df$months), ":")[[i]][1])
    }
    for (i in 1:nrow(plot.df)){
      plot.df$end.mo[i] <- as.numeric(ifelse(length(strsplit(as.character(plot.df$months), ":")[[i]]) == 2,
                                             strsplit(as.character(plot.df$months), ":")[[i]][2],
                                             strsplit(as.character(plot.df$months), ":")[[i]][1]))
    }

    # Unlist these and make sure they are numeric.
    plot.df$start.mo <- unlist(plot.df$start.mo) |> as.numeric()
    plot.df$end.mo <- unlist(plot.df$end.mo) |> as.numeric()

    # Simplified version of significance
    plot.df$simp.sig <- ifelse(plot.df$sig %in% c("*","**","***"), "Sig.","Not sig.")

    # Build a nice title for the plots
    if (prewhiten == TRUE) {
      title <- ifelse(is.null(chrono.name), paste0("Chronology correlations with ", clim.var, " (prewhitened)"),
                      paste0(chrono.name, " chronology correlations with ", clim.var, " (prewhitened)"))
    } else {
      title <- ifelse(is.null(chrono.name), paste0("Chronology correlations with ", clim.var),
                      paste0(chrono.name, " chronology correlations with ", clim.var))
    }

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
                                max.lag, " (overlapping years ", min.y, "-", max.y, ")"))

    out.list <- list(lag.res, lag.df, out.plot)
    names(out.list) <- c("Correlation results", "Correlation data", "Results plots")
    return(out.list)
  } else {
    out.list <- list(lag.res, lag.df)
    names(out.list) <- c("Correlation results", "Correlation data")
    return(out.list)
  }

} # End of function
