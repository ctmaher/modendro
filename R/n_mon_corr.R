#' Flexible monthly aggregate growth-climate cross correlations for exploratory data analysis
#'
#' @description
#' Exploratory data analysis (EDA) function to compute correlations between tree ring data and a
#' monthly climate variable aggregated for every combination (lengths 1:12) of consecutive months
#' inside of a 12-month long "relevant climate period" that could conceivably be relevant to growth
#' in any given year.
#'
#' The ideas behind the relevant climate period are that climate no longer has an effect on radial
#' growth after growth has stopped & that current year's growth could have been influenced by
#' climate in any month AFTER the previous year's growth stopped.
#'
#' The user specifies the beginning month of the relevant climate period - this is the 1st month
#' after radial growth stops (i.e., the 1st month of the fall season). For example, if radial
#' growth typically terminates sometime in September, the user would enter `rel.per.begin = 10`
#' to specify a relevant climatic period that starts in October of the previous year and ends in
#' September of the next year (i.e., the calendar year of growth). In this way the relevant
#' climatic period covers the "water year" months leading up to the growing period and only extends
#' through the months when growth occurs. It is generally nonsensical to include months after
#' radial growth stops - doing so may result in spurious correlations.
#'
#' Compared to methods that force rigidly-defined seasons of a fixed length, this approach should
#' facilitate discovery of potentially more meaningful growth-climate relationships.
#'
#' Fair warning: this is a basic function that will accept any tree ring data and climate data in
#' the proper format. It is the user's responsibility to make sure that your data is appropriate
#' to the analyses.
#'
#' @param rw a tree ring series in "long format" with at least two columns representing the year
#' and the tree ring data. Could be individual tree series or a chronology. Not technically
#' restricted to ring widths.
#' @param rw.col character vector - the colname of the tree ring data series.
#' @param clim a `data.frame` with at least 3 columns: year, month (numeric), and a
#' climate variable.
#' @param clim.var character vector - the colname of the climate variable of interest in the `clim`
#'  data.frame.
#' @param rel.per.begin an integer month representing the beginning of the climatically relevant
#' period to the growth year (always a 12 month period).
#' This will include the "water year" of the calendar year before growth. E.g., 10 for N
#' hemisphere, 4 for S hemisphere. See details below for more info.
#' @param hemisphere a character vector specifying which hemisphere your tree ring data - &
#' climate data - comes from ("N" or "S").
#' Conventions for assigning growth years - and thus aligning tree ring and climate data - are
#'  different for N and S hemisphere.
#' @param agg.fun character vector specifying the function to use for aggregating monthly
#' climate combinations. Options are "mean" or "sum", e.g., for temperature or precipitation data,
#' respectively. Default is "mean".
#' @param max.lag numeric vector specifying how many years of lag to calculate calculations for.
#' Default is 1 year.
#' @param prewhiten logical vector specifying whether or not to convert tree ring & climate time
#' series to ARIMA residuals (aka "prewhitening"). A "best fit" ARIMA model is automatically
#' selected using \code{\link[forecast]{auto.arima}}.
#' This removes autocorrelation in a time series, leaving only the high-frequency variation.
#' This is common practice before using standard methods for cross-correlations. Default is FALSE.
#' @param auto.corr logical vector specifying whether there is temporal autocorrelation in
#' either your tree ring or climate time series (there typically is autocorrelation,
#' unless both are "prewhitened").
#' If TRUE (the default), & corr.method is "spearman" or "kendall", then the
#' \code{\link[corTESTsrd]{corTESTsrd}} function is used to compute modified significance
#' testing to account for autocorrelation (From Lun et al. 2022).
#' Caution! Currently auto.corr = TRUE & corr.method = "Pearson" doesn't make any adjustments.
#' @param corr.method character vector specifying which correlation method to use. Default is
#' `"spearman"`. Options are `c("pearson", "kendall", "spearman")`.
#'  Passes to \code{\link[stats]{cor.test}} or to \code{\link[corTESTsrd]{corTESTsrd}}.
#' @param rw.name character vector - the name of your tree ring series (optional).
#' This is used in the title of your plot. If you produce many plots, this helps keep them
#' identifiable.
#' @param plots logical vector indicating whether or not to produce plots. Default is TRUE.
#' @param silent logical vector indicating whether messages about relevant period and hemisphere
#' conventions will be printed. Default is FALSE.
#'
#' @details
#' Exploring a wide range of plausible growth-climate relationships can be a useful first step once
#' you have a collection of cross-dated tree ring series and have properly detrended them,
#' standardized them, etc.
#'
#' The default correlation test method is Spearman rank correlation. This will be ±equivalent
#' to Pearson for linear relationships, but will also capture any non-linear relationships.
#'
#' A note on tree ring analyses based in the Southern hemisphere:
#' \code{\link{n_mon_corr}} is designed to work in both the Northern and Southern hemispheres.
#' Hemisphere matters for tree ring growth-climate relationships because tree ring formation in the
#' Southern hemisphere typically spans two calendar years (e.g., starting in Nov 2000 and ending
#' in Mar of 2001). It was Schulman's (1956) protocol to assign the earlier calendar year to the
#' tree rings in the Southern hemisphere, i.e., the calendar year in which growth began.
#' \code{\link{n_mon_corr}} assumes your data follows this standard as well.
#' This has implications for how the climate data is aligned with the treering data.
#' The current implementation handles this implicitly by assuming that if `rel.per.begin` is
#' between 1:6, this is a S. hemisphere analysis and the current "growth year" is the same as the
#' calendar year of `rel.per.begin`. If `rel.per.begin` is between 7:12, it is assumed that the
#' is a N. hemisphere analysis and the current "growth year" is the calendar year following
#' `rel.per.begin`. E.g., if `rel.per.begin = 4`, the climatically relevant period will be defined
#' as months `c(4,5,6,7,8,9,10,11,12,1,2,3)` with the calendar year of the FIRST 9 months as the
#' "growth year". If `rel.per.begin = 10`, the climatically relevant period will be defined as
#' months `c(10,11,12,1,2,3,4,5,6,7,8,9)` with the calendar year of the LAST 9 months as the
#' "growth year".
#'
#' Interpreting the plots:
#' The plots show a 12-month sequence of consecutive months on the x-axis & the correlation
#' coefficient on the y-axis. The diamonds indicate the starting month of an n-month aggregate
#' period, small vertical bars the end. Horizontal lines connect the start and end months for
#' periods > 1 month. Significant correlations (as determined by \code{\link[stats]{cor.test}})
#' are shown in black, no significant ones in grey. Plot panel labels (right-hand side of plots)
#' indicate lag years: 0 = current year, -1 = previous year, -2 = 2 years back.
#'
#' @return A 2-4 element list containing data.frames of the correlation results, the data used
#' in the correlations (both prewhitened and raw if prewhiten = TRUE), and the default plots
#'  of the same data.
#'
#' @references
#' Schulman, E. (1956) \emph{Dendroclimatic changes in semiarid America},
#' University of Arizona Press.
#'
#' Lun, D., S. Fischer, A. Viglione, and G. Blöschl. (2022). Significance testing of rank
#' cross-correlations between autocorrelated time series with short-range dependence,
#'  \emph{Journal of Applied Statistics}:1–17.
#'
#' @import ggplot2
#' @import corTESTsrd
#' @importFrom forecast auto.arima
#'
#' @export
#'
#' @examples
#' # Make some synthetic stand-in tree ring data
#' rw.ex <- data.frame(year = 1:50, rw.mm = runif(50, 0.1, 2))
#'
#'
#' # Make some synthetic stand-in climate data
#' mo <- 1:12
#' x <- seq(-4, 4, length.out = 12)
#' gauss.curv <- \(x) {(10/sqrt(2*pi*1.6))*exp(-((x^2)/(2*1.6^2)))}
#' clim.mo <- data.frame(month = mo, clim.var = gauss.curv(x))
#' clim.list <- vector("list", length = 50)
#' for (y in seq_along(clim.list)) {
#'   clim.y <- clim.mo
#'   clim.y[,"clim.var"] <- clim.y[,"clim.var"] + runif(1, min = 0, max = 4)
#'   clim.y$year <- y
#'   clim.list[[y]] <- clim.y
#' }
#' clim <- do.call("rbind", clim.list)
#' n_mon_corr.out <- n_mon_corr(rw = rw.ex,
#' rw.col = "rw.mm",
#' clim = clim,
#' clim.var = "clim.var",
#' rel.per.begin = 3,
#' hemisphere = "S",
#' rw.name = "Synthetic RW")
#'
#' # Take a look at the output data frames
#' head(n_mon_corr.out[["Correlation results"]])
#' head(n_mon_corr.out[["Correlation data"]])
#'
#' ## Method for applying n_mon_corr to each series (or tree) in a rwl file,
#' # consistent with tree-level analyses.
#'
#' # Make some synthetic stand-in tree ring data (in rwl-type format)
#' rwl.ex <- matrix(nrow = 50, ncol = 10)
#' rownames(rwl.ex) <- 1:50
#' colnames(rwl.ex) <- 1:10
#' rwl.ex <- apply(rwl.ex, MARGIN = 2, FUN = \(x) runif(length(x), 0.1, 2)) |> as.data.frame()
#'
#' # Convert rwl to long format using modendro's rwl_longer() function
#' rwl.ex.long <- modendro::rwl_longer(rwl.ex,
#' series.name = "tree",
#' dat.name = "rw.mm",
#' trim = TRUE)
#'
#' # split the data.frame into a list based on ID
#' rwl.ex.long.list <- split(rwl.ex.long, f = rwl.ex.long$tree)
#'
#' # Use mapply to run n_mon_corr for each "tree"
#' # Note that in our example, the clim data is the same for all series.
#' # You will need a different processes if you have climate and tree ring series grouped by site
#' # or plot.
#' # Also note that plots = FALSE and silent = TRUE so that these don't clog the plotting window
#' # and console.
#' n_mon_corr.out.list <- lapply(rwl.ex.long.list, FUN = \(x) {
#' n_mon_corr(rw = x,
#' rw.col = "rw.mm",
#' clim = clim,
#' clim.var = "clim.var",
#' rel.per.begin = 3,
#' hemisphere = "S",
#' plots = FALSE,
#' silent = TRUE)
#' })
#'
#' # Outputs are the same as above, but nested in one more list dimension
#' head(n_mon_corr.out.list[[1]][["Correlation results"]])
#'
#' # It might be helpful to rbind the individual tree outputs into data.frames
#'
#' data.df <- lapply(n_mon_corr.out.list, FUN = \(x, y) {
#' x[["Correlation data"]]
#' }) |>
#'  do.call(what = "rbind")
#'
#' head(data.df)
#'
#' results.df <- mapply(FUN = \(x, y) {
#'  x[["Correlation results"]]$tree <- y
#'  x[["Correlation results"]]
#' }, x = n_mon_corr.out.list, y = names(n_mon_corr.out.list), SIMPLIFY = FALSE) |>
#'  do.call(what = "rbind") |> as.data.frame()
#'
#' head(results.df)
#'
#' # Now it is possible to do things like find out how many trees have significant correlations
#' # for each month combination (these are random series so results are not meaningful):
#' results.df$sig <- ifelse(results.df$p < 0.05, "Sig.","Not sig.")
#' results.agg <- aggregate(tree ~ months + sig, data = results.df, length)
#' results.agg[order(results.agg$tree, decreasing = TRUE),] |> head()


n_mon_corr <- function(rw = NULL,
                       rw.col = "std",
                       clim = NULL,
                       clim.var = NULL,
                       rel.per.begin = NULL,
                       hemisphere = NULL,
                       agg.fun = "mean",
                       max.lag = 1,
                       prewhiten = FALSE,
                       auto.corr = FALSE,
                       corr.method = "spearman",
                       rw.name = NULL,
                       plots = TRUE,
                       silent = FALSE) {
  ## Initial error catching and interactive prompts

  stopifnot(
    "Arg rw or clim are not an object of class 'chron', 'data.frame', or 'matrix'" =
      data.class(rw) %in% "chron" |
      data.class(rw) %in% "data.frame" |
      data.class(rw) %in% "matrix" |
      data.class(clim) %in% "data.frame" |
      data.class(clim) %in% "matrix"
  )

  stopifnot(
    "clim or rw does not have a year column? (name sould start with 'y' or 'Y')" =
      any(substr(colnames(rw), 1, 1) %in% c("Y", "y")) == TRUE |
      any(substr(colnames(clim), 1, 1) %in% c("Y", "y")) == TRUE
  )

  match.test <- clim.var %in% colnames(clim)
  stopifnot("Arg clim.var must match one unique column name in clim" =
              length(match.test[match.test == TRUE]) == 1)

  match.test <- colnames(rw) %in% rw.col
  stopifnot(
    "Arg rw.col must match the name
         of the growth variable in the rw data.frame" =
      length(match.test[match.test == TRUE]) == 1
  )

  if (is.null(rel.per.begin)) {
    cat(
      "You haven't specified the beginning month of the relevant climate period -\n",
      "this is the 1st month after radial growth typically stops in a year\n",
      "(i.e., the 1st month of the fall season at your site).\n",
      "You should have an approximate idea of what month works for your study system.\n"
    )
    rel.per.begin <-
      readline(prompt = "First month of relevant climate period = ") |> as.integer()
  }

  if (is.null(hemisphere)) {
    cat(
      "You haven't specified the hemisphere from which your tree ring series comes from.\n",
      "This is important because there are different conventions for linking growth years\n",
      "to climate years for N vs. S hemispheres"
    )
    hemisphere <-
      readline(prompt = "Enter hemisphere ('N' or 'S') = ")
  }

  stopifnot(
    "Invalid climatically relevant period begin month provided (must be a single integer month)" =
      is.numeric(rel.per.begin) &
      length(rel.per.begin) == 1
  )

  stopifnot(
    "Invalid hemisphere argument provided (must be a character vector & either 'S' or 'N')" =
      is.character(hemisphere) &
      substr(hemisphere, 1, 1) %in% c("s", "S", "N", "n")
  ) # actually more permissive than the error message suggests

  stopifnot("Arg agg.fun must be either 'mean' or 'sum'" =
              agg.fun %in% "mean" |
              agg.fun %in% "sum")

  stopifnot("Arg max.lag must be a numeric vector of length = 1" =
              length(max.lag) == 1 |
              is.numeric(max.lag))

  # accept max.lag inputs that have a negative in front
  if (max.lag < 0) {
    max.lag <- as.numeric(max.lag) |> abs()
  }

  stopifnot(
    "Arg corr.method must be an exact match of one of these: c('pearson','kendall','spearman')" =
      corr.method %in% c("pearson", "kendall", "spearman")
  )


  # make sure that "year" columns are labelled as such
  colnames(clim)[which((substr(
    colnames(clim), start = 1, stop = 1
  )
  %in% c("Y", "y")) == T)] <- "year"

  # same for month
  colnames(clim)[which((substr(
    colnames(clim), start = 1, stop = 1
  )
  %in% c("M", "m")) == T)] <- "month"

  stopifnot(
    "Month variable in climate data not numeric or integer" =
      is.numeric(clim[, "month"]) |
      is.integer(clim[, "month"]) |
      is.double(clim[, "month"])
  )

  # make sure year and month are integers from here on out
  clim[, "year"] <- as.integer(clim[, "year"])
  clim[, "month"] <- as.integer(clim[, "month"])

  # Clarify the year variable
  # If there is no column named "Year" or "year"
  if (any(substr(colnames(rw), 1, 1) %in% c("Y", "y")) == FALSE) {
    rw[, "year"] <- rownames(rw) |> as.numeric() # Assume the rownames contain year
  } else {
    colnames(rw)[which((substr(
      colnames(rw), start = 1, stop = 1
    )
    %in% c("Y", "y")) == T)] <- "year"
    rw[, "year"] <- as.numeric(rw[, "year"])
  }

  # n_mon_corr assumes that all years have all 12 months! If even one month is missing somewhere,
  # this will mess up everything that follows.

  mon.count <- aggregate(month ~ year, data = clim, length)

  if (all(mon.count$month != 12)) {
    paste("Year",
          mon.count$year[mon.count$month < 12],
          "does not have all 12 months represented")
    stop("Not all years in climate data have all 12 months represented")
  } # This doesn't do what I want it too

  # n_mon_corr also assumes absolute regularity (this is true for some of the correlation
  # tests too) in both rw & clim
  rw.year.seq <- rw[, "year"]
  rw.year.seq.diff <- rw.year.seq[order(rw.year.seq)] |> diff()
  clim.year.seq <- unique(clim[, "year"])
  clim.year.seq.diff <-
    clim.year.seq[order(clim.year.seq)] |> diff()
  if (any(rw.year.seq.diff != 1) == TRUE) {
    paste("Year", rw.year.seq[which(rw.year.seq.diff > 1)], "is missing from ring width data.")
  }
  stopifnot("Ring width data does not have complete continuous years in annual steps." =
              all(rw.year.seq.diff == 1) == TRUE
            )

  if (any(clim.year.seq.diff != 1) == TRUE) {
    paste("Year", clim.year.seq[which(clim.year.seq.diff > 1)], "is missing from climate data.")
  }
  stopifnot("Climate data does not have complete continuous years in annual steps." =
              all(clim.year.seq.diff == 1) == TRUE
            )


  # Give a warning & maybe stop the function if there is autocorrelation in the tree ring series
  # There should be a prompt (verbal or otherwise)
  ac.test <- ar(x = na.omit(rw[order(rw[, "year"]), rw.col]))
  if (ac.test$order > 0 & auto.corr == FALSE & prewhiten == FALSE) {
    cat(
      "Autocorrelation detected in rw, recommend choose auto.corr = TRUE &
        corr.method = c('spearman', 'kendall') to avoid spurious correlation results.
        You can also select prewhiten = TRUE to produce ARIMA residuals of both rw
      and clim data.\n"
    )
    # auto.corr <- readline(prompt = "Enter auto.corr (TRUE or FALSE) = ")
    # corr.method <- readline(prompt = "Enter corr.method ('spearman' or 'kendall') = ")
  }

  # Clean up the hemisphere argument if needed
  hemisphere <-
    ifelse(substr(hemisphere, 1, 1) %in% c("n", "N") , "N", "S")

  # Create a chronological sequence of months starting with the numeric month
  # given as the clim.rel.per.begin argument
  mon.seq <-
    rep(1:12, 2) # A list of months representing 2 whole calendar years.
  mon.seq <-
    mon.seq[rel.per.begin:length(mon.seq)][1:12] # 12 months in a row
  # This seq of months will then be relevant to a particular current growth year, though the
  # months span 2 calendar years. According to standard practice, for N hemisphere tree rings, the
  # growth year is the later calendar year. In the S hemisphere, it is the earlier calendar year.

  # Print the relevant climate period
  if (silent == FALSE) {
    cat(
      "You have specified the following months for your relevant climate period in the\n",
      hemisphere,
      "hemisphere:",
      mon.seq,
      "\n"
    )
  }

  # The operations below assume that the climate data is arranged by month, then year.
  # Let's ensure this is the case.
  clim <- clim[order(clim$year, clim$month), ]

  # Define the "growth year" based on the hemisphere argument
  if (hemisphere %in% "S") {
    if (silent == FALSE) {
      message(
        "\nAssuming Southern hemisphere conventions for linking growth years
    and climate years (see ?n_mon_corr for details)\n"
      )
    }
    offset <- rel.per.begin - 1
    clim$growyear <-
      c(rep(min(clim[, "year"]) - 1, offset), clim[, "year"][1:(length(clim[, "year"]) - offset)])

  } else {
    if (silent == FALSE) {
      message(
        "\nAssuming Northern hemisphere conventions for linking growth years
    and climate years (see ?n_mon_corr for details)\n"
      )
    }
    offset <- 12 - rel.per.begin + 1
    clim$growyear <-
      c(clim[, "year"][(offset + 1):length(clim[, "year"])], rep(max(clim[, "year"]) + 1, offset))

  }

  # Give warnings - and stop the function - if someone uses unusual values for rel.per.begin for a
  # given hemisphere. This is to provide guardrails for users who don't understand the rel.per &
  # to minimize the chance that they will be looking at potentially spurious correlations for
  # months AFTER radial growth has ceased in a given year.

  # Find the complete years in the climate data
  clim.complete <- aggregate(month ~ growyear, data = clim, length)
  clim.complete <- clim.complete[clim.complete$month == 12, ]
  clim <- clim[clim[, "growyear"] %in% clim.complete[, "growyear"], ]

  # Find min and max complete years for the correlations, i.e., the complete overlap
  min.y <- max(min(rw[, "year"]), min(clim[, "growyear"]))
  max.y <- min(max(rw[, "year"]), max(clim[, "growyear"]))

  # Create a vector of all possible combinations of months
  # Regular calendar year first
  mos.mat <- expand.grid(mon.seq, mon.seq)

  mos <- apply(mos.mat, MARGIN = 1, FUN = \(x) {
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

  lag.list <- lapply(
    lag.seq,
    FUN = function(l) {
      # Run through all the month aggregates
      cor.res.dfs <- lapply(mos, FUN = \(x) {
        # Aggregate the variable of interest for the given month sequence
        clim.mo <-
          aggregate(formula(paste(clim.var, "growyear", sep = "~")),
                    data = clim[clim$month %in% x,],
                    FUN = \(z) {
                      ifelse(agg.fun %in% "mean", mean(z), sum(z))
                    })

        # The months vector in the desired order.
        month.vec <- ifelse(length(x) > 1,
                            paste(x[1], x[length(x)], sep = ":"),
                            paste(x))

        # attach rw to the climate data
        clim.mo.new <- clim.mo
        clim.mo.new$growyear <- clim.mo.new$growyear + l
        clim.mo.new <-
          merge(clim.mo.new, rw, by.x = "growyear", by.y = "year")
        clim.mo.new$lag <- ifelse(l == 0, paste(l), paste0("-", l))
        clim.mo.new$months <- month.vec

        # Remove any ties from the data
        clim.mo.new <-
          clim.mo.new[which(!duplicated(clim.mo.new[, clim.var])), ]
        clim.mo.new <-
          clim.mo.new[which(!duplicated(clim.mo.new[, rw.col])), ]

        # If we want to convert climate & rw to ARIMA residuals (aka, "prewhitening"), do it here
        if (prewhiten == TRUE) {
          # Save the raw data before prewhitening
          clim.mo.orig <- clim.mo.new

          # 1st the climate
          arima.mod.clim <-
            forecast::auto.arima(clim.mo.new[!is.na(clim.mo.new[, clim.var]), clim.var],
                       seasonal = FALSE)
          clim.mo.new[!is.na(clim.mo.new[, clim.var]), clim.var] <-
            residuals(arima.mod.clim) |> as.numeric()
          # Now the tree rings
          arima.mod.rw <-
            forecast::auto.arima(clim.mo.new[!is.na(clim.mo.new[, rw.col]), rw.col],
                       seasonal = FALSE)
          clim.mo.new[!is.na(clim.mo.new[, rw.col]), rw.col] <-
            residuals(arima.mod.rw) |> as.numeric()
        }

        # Run the correlation test between climate and rw
        if (auto.corr == TRUE) {
          if (corr.method %in% "pearson") {
            ct <- cor.test(clim.mo.new[, clim.var],
                           clim.mo.new[, rw.col],
                           method = corr.method,
                           alternative = "two.sided")
            # put the results together in a data.frame
            result <- data.frame(
              months = month.vec,
              coef = ct$estimate[[1]],
              p = ct$p.value[[1]],
              ci.lo = ct$conf.int[1],
              ci.hi = ct$conf.int[2]
            )
          } else {
            # if spearman or kendall
            ct <-
              corTESTsrd::corTESTsrd(
                clim.mo.new[, clim.var],
                clim.mo.new[, rw.col],
                method = corr.method,
                iid = FALSE,
                alternative = "two.sided"
              )
            # put the results together in a data.frame
            result <- data.frame(months = month.vec,
                                 coef = ct[["rho"]],
                                 p = ct[["pval"]])
          }
        } else {
          ct <-
            cor.test(clim.mo.new[, clim.var], clim.mo.new[, rw.col], method = corr.method)

          # put the results together in a data.frame
          if (corr.method %in% "pearson") {
            result <- data.frame(
              months = month.vec,
              coef = ct$estimate[[1]],
              p = ct$p.value[[1]],
              ci.lo = ct$conf.int[1],
              ci.hi = ct$conf.int[2]
            )
          } else {
            result <- data.frame(
              months = month.vec,
              coef = ct$estimate[[1]],
              p = ct$p.value[[1]]
            )
          }
        }
        # return the correlation results and the data.frame of the merged climate and rw data
        if (prewhiten == TRUE) {
          list(result, clim.mo.new, clim.mo.orig)
        } else {
          list(result, clim.mo.new)
        }
      })

      cor.results <- lapply(cor.res.dfs, FUN = \(x) {
        x[[1]]
      }) |> do.call(what = "rbind")

      cor.results$sig <-
        ifelse(cor.results$p >= 0.05, "Not sig.", "Sig.")

      cor.results$lag <- ifelse(l == 0, paste(l), paste0("-", l))

      cor.df <- lapply(cor.res.dfs, FUN = \(x) {
        x[[2]]
      }) |> do.call(what = "rbind")

      cor.df$lag <- ifelse(l == 0, paste(l), paste0("-", l))

      if (prewhiten == TRUE) {
        cor.raw.df <- lapply(cor.res.dfs, FUN = \(x) {
          x[[3]]
        }) |> do.call(what = "rbind")

        cor.raw.df$lag <- ifelse(l == 0, paste(l), paste0("-", l))

        # Return all the results
        list(cor.results, cor.df, cor.raw.df)
      } else {
        # Return all the results
        list(cor.results, cor.df)
      }
    }
  )

  lag.res <- lapply(lag.list, FUN = \(x) {
    x[[1]] # 1st list element contains the results
  }) |> do.call(what = "rbind")
  # Make the lag a factor
  lag.res$lag <- factor(lag.res$lag, levels = lag.seq * -1)
  # sort correlation results by correlation coef
  lag.res <- lag.res[order(lag.res$coef, decreasing = TRUE), ]
  # order the months as a factor
  lag.res$months <- factor(lag.res$months, levels = mos.fac)

  lag.df <- lapply(lag.list, FUN = \(x) {
    x[[2]] # 2nd list element contains the data used in the correlations
  }) |> do.call(what = "rbind")
  # Make the lag a factor
  lag.df$lag <- factor(lag.df$lag, levels = lag.seq * -1)
  # order the months as a factor
  lag.df$months <- factor(lag.df$months, levels = mos.fac)

  if (prewhiten == TRUE) {
    lag.raw.df <- lapply(lag.list, FUN = \(x) {
      x[[3]] # 2nd list element contains the data used in the correlations
    }) |> do.call(what = "rbind")
    # Make the lag a factor
    lag.raw.df$lag <- factor(lag.raw.df$lag, levels = lag.seq * -1)
    # order the months as a factor
    lag.raw.df$months <- factor(lag.raw.df$months, levels = mos.fac)
  }

  # Set up a data.frame for the plots - add some variables that are useful for plotting, but not
  # for the main output.
  if (plots == TRUE) {
    plot.df <- lag.res
    for (i in 1:nrow(plot.df)) {
      plot.df$start.mo[i] <-
        as.numeric(strsplit(as.character(plot.df$months), ":")[[i]][1])
    }
    for (i in 1:nrow(plot.df)) {
      plot.df$end.mo[i] <-
        as.numeric(ifelse(
          length(strsplit(as.character(
            plot.df$months
          ), ":")[[i]]) == 2,
          strsplit(as.character(plot.df$months), ":")[[i]][2],
          strsplit(as.character(plot.df$months), ":")[[i]][1]
        ))
    }

    # Unlist these and make sure they are numeric.
    plot.df$start.mo <- unlist(plot.df$start.mo) |> as.numeric()
    plot.df$end.mo <- unlist(plot.df$end.mo) |> as.numeric()


    # Build a nice title for the plots
    if (prewhiten == TRUE) {
      title <-
        ifelse(
          is.null(rw.name),
          paste0("Correlations with ", clim.var, " (prewhitened)"),
          paste0(rw.name, " correlations with ", clim.var, " (prewhitened)")
        )
    } else {
      title <-
        ifelse(
          is.null(rw.name),
          paste0("Correlations with ", clim.var),
          paste0(rw.name, " correlations with ", clim.var)
        )
    }

    # Make the plot
    # These 3 lines are to deal with "no visible binding" NOTEs from check()
    x_var1 <- "start.mo"
    x_var2 <- "end.mo"
    col_var <- "sig"
    out.plot <- ggplot2::ggplot(plot.df) +
      ggplot2::geom_point(aes(factor(.data[[x_var1]], levels = mon.seq), # start points
                     coef, color = .data[[col_var]]),
                 shape = 18,
                 size = 2) +
      ggplot2::geom_point(aes(factor(.data[[x_var2]], levels = mon.seq), # end points
                     coef, color = .data[[col_var]]),
                 shape = 124,
                 size = 2) +
      ggplot2::geom_segment(aes(
        x = factor(.data[[x_var1]], levels = mon.seq),
        # lines connecting
        xend = factor( .data[[col_var]], levels = mon.seq),
        y = coef,
        yend = coef,
        color = .data[[col_var]]
      )) +
      ggplot2::scale_x_discrete(breaks = mon.seq, labels = mon.seq) +
      ggplot2::scale_y_continuous(breaks = seq(-1, 1, by = 0.1)) +
      ggplot2::xlab("Month") +
      ggplot2::ylab(paste0("Correlation coefficient\n(", corr.method, ")")) +
      ggplot2::scale_color_manual(name = "", values = c("grey80", "black")) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap( ~ lag, ncol = 1, strip.position	= "right") +
      ggplot2::ggtitle(
        label = title,
        subtitle = paste0(
          "Monthly climate ",
          agg.fun,
          "s with annual lag 0:",
          max.lag,
          " (overlapping years ",
          min.y,
          "-",
          max.y,
          ")"
        )
      )

    if (prewhiten == TRUE) {
      out.list <- list(lag.res, lag.df, lag.raw.df, out.plot)
      names(out.list) <-
        c(
          "Correlation results",
          "Correlation data (prewhitened)",
          "Correlation data (raw)",
          "Results plots"
        )
    } else {
      out.list <- list(lag.res, lag.df, out.plot)
      names(out.list) <-
        c("Correlation results",
          "Correlation data",
          "Results plots")
    }

    return(out.list)
  } else {
    if (prewhiten == TRUE) {
      out.list <- list(lag.res, lag.df, lag.raw.df)
      names(out.list) <-
        c(
          "Correlation results",
          "Correlation data (prewhitened)",
          "Correlation data (raw)"
        )
    } else {
      out.list <- list(lag.res, lag.df)
      names(out.list) <-
        c("Correlation results", "Correlation data")
    }

    return(out.list)
  }

} # End of function
