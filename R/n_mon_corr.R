#' Flexible monthly aggregate growth-climate cross correlations for exploratory data analysis
#'
#' @description
#' Exploratory data analysis (EDA) function to compute correlations between tree ring data and a
#' monthly climate variable aggregated for every combination (lengths 1:12) of consecutive months
#' going back a specified number of years.
#'
#' Compared to methods that force rigidly-defined seasons of a fixed length, this approach should
#' facilitate discovery of potentially more meaningful growth-climate relationships.
#'
#' Fair warning: this is a basic function that will accept any tree ring data and climate data in
#' the proper format. It is the user's responsibility to make sure that your data is appropriate
#' to the analyses.
#'
#' @param rwl A rwl-type data.frame (e.g., read in by \code{\link[dplR]{read.rwl}}). Essentially a
#' data.frame with columns names as series IDs and years as rownames.
#' @param clim a `data.frame` with at least 3 columns: year, month (numeric), and a
#' climate variable.
#' @param clim.var character vector - the colname of the climate variable of interest in the `clim`
#'  data.frame.
#' @param gro.period.end the last month in which you expect growth to occur for your study species
#' in your study region. Not crucial in this version - only draws a line in the output plot.
#' @param agg.fun character vector specifying the function to use for aggregating monthly
#' climate combinations. Options are "mean" or "sum", e.g., for temperature or precipitation data,
#' respectively. Default is "mean".
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
#' @param make.plot logical vector indicating whether or not to produce a plot. Default is TRUE.
#' You will get a warning if you have < 10 tree-ring series. Ideally you have > 50.
#' @param group.IDs.df an optional data.frame with 2 columns: "series", representing the names of
#' the tree-ring series (and matching the colnames of rwl) and a group.var (name is your
#' specification), representing a grouping (e.g., site, plot) variable.
#' @param group.var a character vector specifying a grouping variable (e.g., site, plot).
#' Must exist and be identical in clim and group.ID.df.
#'
#' @details
#' Exploring a wide range of plausible grrowth-climate relationships can be a useful first step once
#' you have a collection of cross-dated tree ring series.
#'
#' The default correlation test method is Spearman rank correlation. This will be ±equivalent
#' to Pearson for linear relationships, but will also capture any non-linear relationships. As an
#' additional precaution, \code{\link{n_mon_corr}} uses the \code{\link[corTESTsrd]{corTESTsrd}}
#' method from Lun et al. (2022) to reduce the type I error rate  (i.e., you think there is a
#' significant cross-correlation when there isn't one), which is inflated if there is
#' autocorrelation in your tree-ring or climate series. When autocorrelation is absent, the method
#' gives similar results as the classical significance test. For these reasons,
#' \code{\link{n_mon_corr}} always uses the Lun et al. method for Spearman and Kendall rank
#' correlations. This is done here as a precaution but also as an allowance so that you can examine
#' relationships between lower-frequency (which are autocorrelated) elements of tree-ring and
#' climate series as well as between the more common high-frequency signals (e.g., "prewhitened"
#' series). You may use Pearson correlations in \code{\link{n_mon_corr}}, but be forewarned that
#' there is no adjustment for autocorrelation in the Pearson correlation tests! Therefore, it
#' is general best to use `corr.method = "pearson"` for comparison only, rather than actual results.
#'
#' A note on tree ring analyses based in the Southern hemisphere:
#' \code{\link{n_mon_corr}} is designed to work in both the Northern and Southern hemispheres.
#' Hemisphere matters for tree ring growth-climate relationships because tree ring formation in the
#' Southern hemisphere typically spans two calendar years (e.g., starting in Nov 2000 and ending
#' in Mar of 2001). It was Schulman's (1956) protocol to assign the earlier calendar year to the
#' tree rings in the Southern hemisphere, i.e., the calendar year in which growth began.
#' \code{\link{n_mon_corr}} assumes your data follows this standard as well. This has implications
#' for how the climate data is aligned with the treering data. In this version, the way this is
#' handled is that a "lag +1" year is added to suite of correlations so that the moving windows of
#' consecutive months can extend through the growing season (and into the next calendar year).
#'
#'
#' @return A 2-4 element list containing data.frames of the correlation results, the moving-window
# climate data used in the correlations (both prewhitened and raw if prewhiten = TRUE), the
#' prewhitened tree-ring data (if prewhiten = TRUE), and the default plot.
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
#' @importFrom grDevices hcl.colors
#'
#' @export
#'
#' @examples
#'## Bring in some real tree-ring and climate data
#' # Tree-ring data from Perkins & Swetnam 1996 (https://doi.org/10.1139/x26-241)
#' data(PerkinsSwetnam96)
#' # PRISM (https://www.prism.oregonstate.edu/) time series extracted for each site then averaged
#' # over all sites.
#' data(idPRISM)
#'
#' # Monthly average temperature
#' PS_gro_Tavg <- n_mon_corr(rwl = PerkinsSwetnam96,
#'                          clim = idPRISM,
#'                          clim.var = "Tavg.C",
#'                          agg.fun = "mean",
#'                          max.lag = 2,
#'                          hemisphere = "N",
#'                          prewhiten = TRUE,
#'                          corr.method = "spearman",
#'                          gro.period.end = 9,
#'                          make.plot = TRUE,
#'                          group.IDs.df = NULL,
#'                          group.var = NULL)
#' names(PS_gro_Tavg)
#' PS_gro_Tavg$`Results plot`
#' # April has the strongest signal (largest % of significant correlations), and it is a negative
#' # relationship.
#'
#'
#' # Monthly total precipitation
#' # Note that agg.fun = "sum" for precipitation data.
#' PS_gro_PPT <- n_mon_corr(rwl = PerkinsSwetnam96,
#'                         clim = idPRISM,
#'                         clim.var = "PPT.mm",
#'                         agg.fun = "sum",
#'                         max.lag = 2,
#'                         hemisphere = "N",
#'                         prewhiten = TRUE,
#'                         corr.method = "spearman",
#'                         gro.period.end = 9,
#'                         make.plot = TRUE,
#'                         group.IDs.df = NULL,
#'                         group.var = NULL)
#' names(PS_gro_PPT)
#' PS_gro_PPT$`Results plot`
#' # Strongest signal is Oct-Jan total precip.
#'
#'
#' ## Demonstrate grouped data
#' data(idPRISMgroup)
#' data(PSgroupIDs)
#'
#' PS_gro_Tavg_grouped <- n_mon_corr(rwl = PerkinsSwetnam96,
#'                                  clim = idPRISMgroup,
#'                                  clim.var = "Tavg.C",
#'                                  agg.fun = "mean",
#'                                  max.lag = 2,
#'                                  hemisphere = "N",
#'                                  prewhiten = TRUE,
#'                                  corr.method = "spearman",
#'                                  gro.period.end = 9,
#'                                  make.plot = TRUE,
#'                                  group.IDs.df = PSgroupIDs,
#'                                  group.var = "site")
#' names(PS_gro_Tavg_grouped)
#' PS_gro_Tavg_grouped$`Results plot`
#' # Similar result, but not identical - now start.month = Apr with a 2-month moving window
#' # The climate data is slightly different for each site, so some differences in results are not
#' # surprising.



n_mon_corr <- function(rwl = NULL,
                       clim = NULL,
                       clim.var = NULL,
                       agg.fun = "mean",
                       max.lag = 1,
                       hemisphere = NULL,
                       prewhiten = TRUE,
                       corr.method = "spearman",
                       gro.period.end = NULL,
                       make.plot = TRUE,
                       group.IDs.df = NULL,
                       group.var = NULL) {
  ############ Initial basic input error catching

  ###### rwl
  stopifnot(
    "Arg rwl is not an object of class 'rwl' or 'data.frame'" =
      data.class(rwl) %in% "rwl" |
      data.class(rwl) %in% "data.frame"
  )

  # Add a year variable to rwl
  # If there is no column named "Year" or "year" already
  if (any(substr(colnames(rwl), 1, 1) %in% c("Y", "y")) == FALSE) {
    rwl[, "year"] <- rownames(rwl) |> as.numeric() # Assume the rownames contain year
  } else {
    colnames(rwl)[which((substr(
      colnames(rwl), start = 1, stop = 1
    )
    %in% c("Y", "y")) == T)] <- "year"
    rwl[, "year"] <- as.numeric(rwl[, "year"])
  }

  ###### clim
  stopifnot(
    "Arg clim is not an object of class 'data.frame', or 'matrix'" =
      data.class(clim) %in% "data.frame" |
      data.class(clim) %in% "matrix"
  )

  stopifnot(
    "clim does not have a year column? (colname sould start with 'y' or 'Y')" =
      any(substr(colnames(clim), 1, 1) %in% c("Y", "y")) == TRUE
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



  ###### clim.var
  match.test <- clim.var %in% colnames(clim)
  stopifnot("Arg clim.var must match one unique column name in clim" =
              length(match.test[match.test == TRUE]) == 1)

  ###### agg.fun
  stopifnot("Arg agg.fun must be either 'mean' or 'sum'" =
              agg.fun %in% "mean" |
              agg.fun %in% "sum")

  ###### max.lag
  stopifnot("Arg max.lag must be an integer vector of length = 1" =
              length(max.lag) == 1 |
              is.integer(max.lag))

  # Accept max.lag inputs that have a negative in front
  # but change them to positive
  if (max.lag < 0) {
    max.lag <- as.integer(max.lag) |> abs()
  }

  stopifnot("Arg max.lag must have an absolute value >= 1" =
              max.lag >= 1)

  ###### hemisphere
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
    "Invalid hemisphere argument provided (must be a character vector & either 'S' or 'N')" =
      is.character(hemisphere) &
      substr(hemisphere, 1, 1) %in% c("s", "S", "N", "n")
  ) # actually more permissive than the error message suggests

  # Clean up the hemisphere argument if needed
  hemisphere <-
    ifelse(substr(hemisphere, 1, 1) %in% c("n", "N") , "N", "S")

  ###### prewhiten
  stopifnot("Arg prewhiten must be a logical vector" =
              is.logical(prewhiten))

  ###### corr.method
  stopifnot(
    "Arg corr.method must be an exact match of one of these: c('pearson','kendall','spearman')" =
      corr.method %in% c("pearson", "kendall", "spearman")
  )

  ###### gro.period.end
  stopifnot(
    "Invalid gro.period.end provided (must be a single integer month)" =
      is.numeric(gro.period.end) |
      is.integer(gro.period.end) |
      length(gro.period.end) == 1 |
      any(gro.period.end == 1:12) == TRUE
  )

  if (is.null(gro.period.end)) {
    cat(
      "You haven't specified the gro.period.end - the last month you expect growth is possible -\n",
      "You should have an approximate idea of what month works for your study system.\n"
    )
    gro.period.end <-
      readline(prompt = "Last month of radial growth = ") |> as.integer()
  }

  ###### make.plot
  stopifnot("Arg make.plot must be a logical vector" =
              is.logical(make.plot))

  ###### group.IDs.df

  # Checks for multiple climate datasets when group.IDs.df is not specified,
  # or the opposite
  multi.clim.check <- aggregate(formula(paste0(clim.var," ~ month + year")), data = clim, length)
  stopifnot(
    "Clim data has groups (multiple observations per month+year) but no group.IDs.df supplied or
    group.IDs.df was supplied and clim data has no groups" =
      # All month+year combos have 1 obs
      (all(multi.clim.check[, clim.var] == 1) &
         # group.IDs.df does not exist
         is.null(group.IDs.df)) |
      # at least some month+year combos have > 1 obs
      (!all(multi.clim.check[, clim.var] == 1) &
         # group.IDs.df does exist
         !is.null(group.IDs.df))
  )


  # If there is a group.IDs.df data.frame check that all series names are included,
  # and other things that must be true.

  if (!is.null(group.IDs.df)) {
    stopifnot(
      "Arg group.IDs.df is not an object of class 'data.frame', or 'matrix'" =
        data.class(group.IDs.df) %in% "data.frame" |
        data.class(group.IDs.df) %in% "matrix"
    )

    ##### group.var
    stopifnot(
      "group.IDs.df supplied but there is no column matching group.var in clim data" =
        any(colnames(clim) %in% group.var)
    )

    match.test2 <- group.var %in% colnames(group.IDs.df)
    stopifnot("Arg group.var must match one unique column name in group.IDs.df" =
                length(match.test2[match.test2 == TRUE]) == 1)

    stopifnot("group.IDs.df must contain a 'series' column" =
                any(colnames(group.IDs.df) %in% "series"))

    stopifnot(
      "group.IDs.df 'series' does not completely match the column names in rwl" =
        any(colnames(group.IDs.df) %in% "series")
    )

  }

  ############ Now check more complex and interactive things

  ###### rwl & clim
  # Generally must overlap in their years - this does not check each series
  year.overlap <- rwl[, "year"] %in% unique(clim[, "year"])
  stopifnot("rwl & clim have no overlap in their years" =
              any(year.overlap))

  # Give a warning if the overlap is less than 25 - general, does not check each series
  if (length(year.overlap[year.overlap == TRUE]) < 25) {
    warning(
      "Less than 25 years of overlap between rwl & clim - be cautious when interpreting
            correlations"
    )
  }

  # Check for the overlap between each series
  clim.span <- unique(clim[, "year"])
  rwl.span <- rwl[, "year"]
  check.series <- apply(rwl[, !(colnames(rwl) %in% "year")], MARGIN = 2, FUN = \(series) {
    this.series <- data.frame(series = series, year = rwl.span)
    na.omit(this.series[this.series$year %in% clim.span,])
  })

  # Identify and remove the series that have less than 4 overlap
  zero.series <- sapply(check.series, FUN = \(x) nrow(x) < 4)
  if (any(zero.series)) {
    message(paste0("The following tree-ring series have < 4 years overlap with clim data
    and will be removed from rwl:\n", paste0(names(zero.series[zero.series == TRUE]),
                                             collapse = ", ")))
    rwl <- rwl[, !(colnames(rwl) %in% names(zero.series[zero.series == TRUE]))]
  }

  # Warn the user about overlaps less than 25 years
  short.series <- sapply(check.series, FUN = \(x) nrow(x) < 25 & nrow(x) >= 4)
  if (any(short.series)) {
    message(paste0("The following tree-ring series have < 25 years overlap with clim data.
    Interpret correlations cautiously.\n", paste0(names(short.series[short.series == TRUE]),
                                                  collapse = ", ")))
  }

  # n_mon_corr & moving_win_multi assume continuity and regularity of all years and months.
  # If even one month or year is missing somewhere, this will mess up everything that follows.
  # NOTE: These tests need special accommodation for when there are groups in the data

  if (is.null(group.IDs.df)) {
    # Because of tests above, we can assume that having no group.IDs.df means we have no group
    # structure at all

    all.yrs <- clim[, "year"] |> unique()
    stopifnot(
      "One or more years are missing in the clim data (observations not regular & continuous)" =
        length(all.yrs) == length(min(all.yrs):max(all.yrs))
    )

    mo.count <- aggregate(month ~ year, data = clim, length)
    stopifnot(
      "One or more years are missing months in the clim data (observations not regular
      & continuous)" =
        all(mo.count[, "month"] == 12)
    )
  } else {
    # If we do have group structure, we must check for missing years/months in each group.
    all.yrs <- aggregate(formula(paste("year ~ ", group.var)),
                         data = clim,
                         FUN = \(x) length(unique(x)))
    all.yrs$min.yr <- aggregate(formula(paste("year ~ ", group.var)),
                                data = clim,
                                FUN = \(x) min(unique(x)))[, 2]
    all.yrs$max.yr <- aggregate(formula(paste("year ~ ", group.var)),
                                data = clim,
                                FUN = \(x) max(unique(x)))[, 2]
    all.yrs$act.len <- apply(all.yrs, MARGIN = 1, FUN = \(x) {
      x["min.yr"]:x["max.yr"] |> length()
    })

    stopifnot(
      "One of more groups is missing one or more years in the clim data (observations
              not regular & continuous)" =
        all(all.yrs$year == all.yrs$act.len)
    )

    mo.count <- aggregate(formula(paste("month ~ year +", group.var)), data = clim, length)
    stopifnot(
      "One or more groups has one or more years that are missing months in the clim data
      (observations not regular & continuous)" =
        all(mo.count[, "month"] == 12)
    )
  }


  # n_mon_corr also assumes absolute regularity (this is true for some of the correlation
  # tests too) in the rwl
  rwl.year.seq <- rwl[, "year"]
  rwl.year.seq.diff <- rwl.year.seq[order(rwl.year.seq)] |> diff()

  stopifnot(
    "Ring width data does not have complete continuous years in annual steps." =
      all(rwl.year.seq.diff == 1) == TRUE
  )


  # The operations below assume that the climate data is arranged by month, then year.
  # Let's ensure this is the case.
  # We have already specified these as integers above.
  clim <- clim[order(clim$year, clim$month), ]

  # Find the complete years in the climate data - this is a redundant safety check,
  # as we should have regulated the continuity in both climate and tree-ring data.
  # clim.complete <- aggregate(month ~ year, data = clim, length)
  # clim.complete <- clim.complete[clim.complete$month == 12, ]
  # clim <- clim[clim[, "year"] %in% clim.complete[, "year"], ]

  # I think here I have to split the data into the different groups if there are
  # is a group.IDs.df

  if (is.null(group.IDs.df)) {
    cor.res.list <- multi_clim_gro_corr(
      rwl.group = rwl,
      clim.group = clim,
      clim.var = clim.var,
      group.IDs.df = group.IDs.df,
      group.var = group.var,
      gro.period.end = gro.period.end,
      agg.fun = agg.fun,
      max.lag = max.lag,
      prewhiten = prewhiten,
      corr.method = corr.method,
      hemisphere = hemisphere
    )

  } else {
    # if there are different climate groups

    clim.group.split <- split(clim, f = clim[, group.var])

    cor.res.list.group <- lapply(clim.group.split, FUN = \(clim.group) {
      this.group <- clim.group[, group.var] |> unique()
      these.series <- group.IDs.df[group.IDs.df[, group.var] %in% this.group, "series"]

      multi_clim_gro_corr(
        rwl.group = rwl[, colnames(rwl)[colnames(rwl) %in% c(these.series, "year")]],
        clim.group = clim.group,
        clim.var = clim.var,
        group.IDs.df = group.IDs.df,
        group.var = group.var,
        gro.period.end = gro.period.end,
        agg.fun = agg.fun,
        max.lag = max.lag,
        prewhiten = prewhiten,
        corr.method = corr.method,
        hemisphere = hemisphere
      )
    })
    # This has to be reassembled after splitting

    # Get the names of the data.frames within the list within the list
    out.list.names <- names(cor.res.list.group[[1]])
    # Put the pieces back together across list elements
    cor.res.list <- lapply(out.list.names, FUN = \(element.name) {
      # The last element is an rwl, which must do cbind instead of rbind.
      if (substr(element.name, 1, 3) %in% "rwl") {
        lapply(cor.res.list.group, FUN = \(group) {
          group[[element.name]]
        }) |>
          unname() |>
          do.call(what = "cbind")
      } else {
        lapply(cor.res.list.group, FUN = \(group) {
          group[[element.name]]
        }) |> do.call(what = "rbind")
      }
    })

    names(cor.res.list) <- out.list.names
  }

  if (make.plot == TRUE) {
    sig.only <- cor.res.list[["cor.res.dat"]][cor.res.list[["cor.res.dat"]]$p <= 0.05, ]
    res.agg <- aggregate(coef ~ start.month + win.len + lag + dir,
                         data = sig.only,
                         FUN = \(x) length(x),
                         drop = FALSE) # make sure the black combinations still appear
    # Replace NAs with 0s
    res.agg$coef <- ifelse(is.na(res.agg$coef), 0, res.agg$coef)

    lag.levels <- res.agg$lag |> unique()

    res.agg$lag <- factor(res.agg$lag, levels = lag.levels[order(as.numeric(lag.levels))])
    res.agg$dir <- factor(res.agg$dir, levels = c("Pos.", "Neg."))
    # Calculate the percentage of significant correlations
    res.agg$prop.sig <- (res.agg$coef / length(unique(cor.res.list[["cor.res.dat"]][, "series"]))) *
      100

    # Make a plot.
    # These 4 lines are to deal with "no visible binding" NOTEs from check()
    x_var <- "start.month"
    y_var <- "prop.sig"
    col_var <- "win.len"
    x.intercept <- "xint"

    out.plot <- ggplot2::ggplot(res.agg,
                                ggplot2::aes(.data[[x_var]], .data[[y_var]],
                                             color = as.factor(.data[[col_var]]))) +
      ggplot2::ggtitle(paste0(ifelse(corr.method %in% "pearson", "Pearson ",
                                    ifelse(corr.method %in% "spearman", "Spearman ",
                                           "Kendall ")),
                             "correlations between ",
                             ifelse(prewhiten == TRUE, "prewhitened ", ""),
                             "tree-ring and ",
                             clim.var,
                             " series")) +
      ggplot2::scale_color_manual("Moving window\nlength\n(n months)",
                                  values = grDevices::hcl.colors(12, palette = "Spectral")) +
      ggplot2::geom_line() +
      #ggplot2::geom_point() +
      ggplot2::facet_grid(dir ~ lag, switch = "x") +
      ggplot2::geom_vline(
        data = data.frame(
          xint = gro.period.end,
          lag = factor(ifelse(hemisphere == "S", "+1", "0")),
          levels = lag.levels[order(as.numeric(lag.levels))]
        ),
        ggplot2::aes(xintercept = .data[[x.intercept]]),
        color = "white"
      ) +
      ggplot2::scale_x_continuous(breaks = c(1:12)) +
      ggplot2::ylab(
        paste(
          "Percentage of",
          length(unique(cor.res.list[["cor.res.dat"]][, "series"])),
          "total series\nrecording significant correlations"
        )
      ) +
      ggplot2::xlab("Start month") +
      ggplot2::coord_cartesian(ylim = c(0, 100)) +
      ggplot2::theme_dark() +
      ggplot2::theme(
        panel.spacing.x = ggplot2::unit(-0.1, "lines"),
        panel.background = ggplot2::element_rect(fill = "black"),
        plot.background = ggplot2::element_rect(fill = "black"),
        legend.background = ggplot2::element_rect(fill = "black"),
        panel.grid = ggplot2::element_line(color = "grey40"),
        legend.text = ggplot2::element_text(color = "white"),
        legend.title = ggplot2::element_text(color = "white"),
        axis.title = ggplot2::element_text(color = "white"),
        axis.text = ggplot2::element_text(color = "white"),
        plot.title = ggplot2::element_text(color = "white"),
        legend.position = "top"
      )

  }


  if (prewhiten == TRUE) {
    if (make.plot == TRUE) {
      out.list <- list(cor.res.list[["cor.res.dat"]],
                       cor.res.list[["clim.dat.pw"]],
                       cor.res.list[["clim.dat"]],
                       cor.res.list[["rwl.dat"]],
                       out.plot)
      names(out.list) <-
        c(
          "Correlation results",
          "Climate data (prewhitened)",
          "Climate data (raw)",
          "Ring-width series (prewhitened)",
          "Results plot"
        )
    } else {
      out.list <- list(cor.res.list[["cor.res.dat"]],
                       cor.res.list[["clim.dat.pw"]],
                       cor.res.list[["clim.dat"]],
                       cor.res.list[["rwl.dat"]])
      names(out.list) <-
        c(
          "Correlation results",
          "Climate data (prewhitened)",
          "Climate data (raw)",
          "Ring-width series (prewhitened)"
        )
    }

  } else {
    if (make.plot == TRUE) {
      out.list <- list(cor.res.list[["cor.res.dat"]], cor.res.list[["clim.dat"]], out.plot)
      names(out.list) <-
        c("Correlation results",
          "Climate data (raw)",
          "Results plot")

    } else {
      out.list <- list(cor.res.list[["cor.res.dat"]], cor.res.list[["clim.dat"]])
      names(out.list) <-
        c("Correlation results", "Correlation data")
    }
  }

  if ((ncol(rwl) - 1) < 10 & make.plot == TRUE) {
    warning("Number of tree-ring series is very low (< 10),
              results plots likely not interpretable")
  }

  # Lastly, order the Correlation results based on the highest correlation
  out.list[["Correlation results"]] <-
    out.list[["Correlation results"]][
      order(out.list[["Correlation results"]]$coef, decreasing = T), ]

  return(out.list)

} # End of function
