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
                       group.var = NULL
) {


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
  stopifnot(
    "Arg prewhiten must be a logical vector" =
      is.logical(prewhiten)
  )

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
  stopifnot(
    "Arg make.plot must be a logical vector" =
      is.logical(make.plot)
  )

  ###### group.IDs.df

  # Checks for multiple climate datasets when group.IDs.df is not specified
  multi.clim.check <- aggregate(clim.var ~ month + year, data = clim, length)
  stopifnot(
    "clim data has multiple observations per month+year yet no group.IDs.df supplied" =
      any(multi.clim.check$clim.var == 1) & # at least one month+year doesn't have 1 obs
      is.null(group.IDs.df) # group.IDs.df does not exist
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
    stopifnot("group.IDs.df supplied but there is no column matching group.var in clim data" =
                any(colnames(clim) %in% group.var)
    )

    match.test2 <- group.var %in% colnames(group.IDs.df)
    stopifnot("Arg group.var must match one unique column name in group.IDs.df" =
                length(match.test2[match.test2 == TRUE]) == 1
    )

    stopifnot("group.IDs.df must contain a 'series' column" =
                any(colnames(group.IDs.df) %in% "series")
    )

    stopifnot("group.IDs.df 'series' does not completely match the column names in rwl" =
                any(colnames(group.IDs.df) %in% "series")
    )

  }

  ############ Now check more complex and interactive things

  ###### rwl & clim
  # Generally must overlap in their years - this does not check each series
  year.overlap <- rwl[, "year"] %in% unique(clim[, "year"])
  stopifnot(
    "rwl & clim have no overlap in their years" =
      any(year.overlap)
  )

  # Give a warning if the overlap is less than 25 - general, does not check each series
  if (length(year.overlap[year.overlap == TRUE]) < 25) {
    warning("Less than 25 years of overlap between rwl & clim - be cautious when interpreting
            correlations")
  }

  # n_mon_corr & moving_win_multi assume continuity and regularity of all years and months.
  # If even one month or year is missing somewhere, this will mess up everything that follows.
  # NOTE: These tests need special accommodation for when there are groups in the data

  if (is.null(group.IDs.df)) {
    # Because of tests above, we can assume that having no group.IDs.df means we have no group
    # structure at all

    all.yrs <- clim[,"year"] |> unique()
    stopifnot(
      "One or more years are missing in the clim data (observations not regular & continuous)" =
        length(all.yrs) == length(min(all.yrs):max(all.yrs))
    )

    mo.count <- aggregate(month ~ year, data = clim, length)
    stopifnot(
      "One or more years are missing months in the clim data (observations not regular
      & continuous)" =
        all(mo.count[,"month"] == 12)
    )
  } else { # If we do have group structure, we must check for missing years/months in each group.
    all.yrs <- aggregate(formula(paste("year ~ ", group.var)),
                         data = clim,
                         FUN = \(x) length(unique(x)))
    all.yrs$min.yr <- aggregate(formula(paste("year ~ ", group.var)),
                                data = clim,
                                FUN = \(x) min(unique(x)))[,2]
    all.yrs$max.yr <- aggregate(formula(paste("year ~ ", group.var)),
                                data = clim,
                                FUN = \(x) max(unique(x)))[,2]
    all.yrs$act.len <- apply(all.yrs, MARGIN = 1, FUN = \(x) {
      x["min.yr"]:x["max.yr"] |> length()
    })

    stopifnot("One of more groups is missing one or more years in the clim data (observations
              not regular & continuous)" =
                all(all.yrs$year == all.yrs$act.len)
    )

    mo.count <- aggregate(formula(paste("month ~ year +", group.var)), data = clim, length)
    stopifnot(
      "One or more groups has one or more years that are missing months in the clim data
      (observations not regular & continuous)" =
        all(mo.count[,"month"] == 12)
    )
  }


  # n_mon_corr also assumes absolute regularity (this is true for some of the correlation
  # tests too) in the rwl
  rwl.year.seq <- rwl[, "year"]
  rwl.year.seq.diff <- rwl.year.seq[order(rwl.year.seq)] |> diff()

  stopifnot("Ring width data does not have complete continuous years in annual steps." =
              all(rwl.year.seq.diff == 1) == TRUE
  )


  # At some point I need to check that there is sufficient overlap between rwl & climate.



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


    cor.res.list <- multi_clim_gro_corr(rwl.group = rwl,
                                        clim.group = clim,
                                        clim.var = clim.var,
                                        group.IDs.df = group.IDs.df,
                                        group.var = group.var,
                                        gro.period.end = gro.period.end,
                                        agg.fun = agg.fun,
                                        max.lag = max.lag,
                                        prewhiten = prewhiten,
                                        corr.method = corr.method,
                                        hemisphere = hemisphere)

  } else { # if there are different climate groups

    clim.group.split <- split(clim, f = clim[, group.var])

    cor.res.list <- lapply(clim.group.split, FUN = \(clim.group) {
      this.group <- clim.group[, group.var] |> unique()
      these.series <- group.IDs.df[group.IDs.df[, group.var] %in% this.group, "series"]

      multi_clim_gro_corr(rwl.group = rwl[, c(these.series, "year")],
                          clim.group = clim.group,
                          clim.var = clim.var,
                          group.IDs.df = group.IDs.df,
                          group.var = group.var,
                          gro.period.end = gro.period.end,
                          agg.fun = agg.fun,
                          max.lag = max.lag,
                          prewhiten = prewhiten,
                          corr.method = corr.method,
                          hemisphere = hemisphere)
    })
    cor.res.list[["cor.res.dat"]] <- do.call(what = "rbind", cor.res.list[["cor.res.dat"]])
    cor.res.list[["cor.res.dat"]] <- merge(cor.res.list[["cor.res.dat"]],
                                           group.IDs.df, by = "series")
  }

  if (make.plot == TRUE) {

    sig.only <- cor.res.list[["cor.res.dat"]][cor.res.list[["cor.res.dat"]]$p <= 0.05,]
    res.agg <- aggregate(coef ~ start.month + win.len + lag + dir,
                         data = sig.only,
                         FUN = \(x) length(x))

    lag.levels <- res.agg$lag |> unique()

    res.agg$lag <- factor(res.agg$lag, levels = lag.levels[order(as.numeric(lag.levels))])
    res.agg$dir <- factor(res.agg$dir, levels = c("Pos.", "Neg."))
    # Calculate the percentage of significant correlations
    res.agg$prop.sig <- (res.agg$coef / length(unique(cor.res.list[["cor.res.dat"]][,"series"]))) * 100

    # Make a plot.
    out.plot <- ggplot2::ggplot(res.agg,
                                ggplot2::aes(start.month, prop.sig, color = as.factor(win.len))) +
      ggplot2::scale_color_manual("Moving window\nlength\n(n months)",
                                  values = hcl.colors(12, palette = "Spectral")) +
      ggplot2::geom_line() +
      ggplot2::facet_grid(dir ~ lag) +
      ggplot2::geom_vline(data = data.frame(xint = gro.period.end,
                                            lag = factor(ifelse(hemisphere == "S",
                                                                "+1",
                                                                "0")),
                                            levels = lag.levels[order(as.numeric(lag.levels))]),
                          ggplot2::aes(xintercept = xint),
                          color = "white") +
      ggplot2::scale_x_continuous(breaks = c(1:12)) +
      ggplot2::ylab(paste("Percentage of",
                          length(unique(cor.res.list[["cor.res.dat"]][,"series"])),
                          "total series\nrecording significant correlations")) +
      ggplot2::xlab("Start month") +
      ggplot2::coord_cartesian(ylim = c(0, 100)) +
      ggplot2::theme_dark() +
      ggplot2::theme(panel.spacing.x = ggplot2::unit(-0.1, "lines"),
                     panel.background = ggplot2::element_rect(fill = "black"),
                     plot.background = ggplot2::element_rect(fill = "black"),
                     legend.background = ggplot2::element_rect(fill = "black"),
                     panel.grid = ggplot2::element_line(color = "grey40"),
                     legend.text = ggplot2::element_text(color = "white"),
                     legend.title = ggplot2::element_text(color = "white"),
                     axis.title = ggplot2::element_text(color = "white"),
                     axis.text = ggplot2::element_text(color = "white"),
                     legend.position = "top")

  }


  if (prewhiten == TRUE) {

    if (make.plot == TRUE) {
      out.list <- list(cor.res.list[["cor.res.dat"]], cor.res.list[["clim.dat.pw"]],
                       cor.res.list[["clim.dat"]], cor.res.list[["rwl.dat"]], out.plot)
      names(out.list) <-
        c(
          "Correlation results",
          "Climate data (prewhitened)",
          "Climate data (raw)",
          "Ring-width series (prewhitened)",
          "Results plot"
        )
    } else {
      out.list <- list(cor.res.list[["cor.res.dat"]], cor.res.list[["clim.dat.pw"]],
                       cor.res.list[["clim.dat"]], cor.res.list[["rwl.dat"]])
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
      out.list <- list(cor.res.list[["cor.res.dat"]],
                       cor.res.list[["clim.dat"]],
                       out.plot)
      names(out.list) <-
        c("Correlation results",
          "Climate data (raw)",
          "Results plot"
        )

    } else {
      out.list <- list(cor.res.list[["cor.res.dat"]],
                       cor.res.list[["clim.dat"]])
      names(out.list) <-
        c("Correlation results",
          "Correlation data")
    }
  }

  if((ncol(rwl) - 1) < 10 & make.plot == TRUE) {
    warning("Number of tree-ring series is very low (< 10),
              results plots likely not interpretable")
  }

  return(out.list)

} # End of function
