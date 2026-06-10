#' Disturbance detrending of tree-ring series based on Druckenbrod et al. 2024
#'
#'
#' @description
#' \code{\link{d_detrend}} implements the disturbance detrending algorithm of Druckenbrod et al.
#' (2024). This is an updated version of the same type of analysis as \code{\link{ci_detect}} - the
#' step-by-step identification and removal of disturbance events tree-ring series from closed-canopy
#' forests. The basic assumption is that rapid changes in growth are likely driven by stand
#' dynamics, not climate. \code{\link{d_detrend}} was adapted from MatLab code provided by
#' D. Druckenbrod.
#'
#' This method combines Nowacki & Abrams's (1997) radial growth averaging to detect
#' disturbances with the intervention detection techniques of Druckenbrod et al. (2013),
#' Rydval et al. (2016), and Rydval et al. (2018) to fit curves to the disturbance event periods and
#' remove them from tree-ring series. By default, events are fit with age-dependent splines
#' (Melvin et al. 2007).
#'
#' The `modendro` implementation includes additional flexibility to try different curve fitting
#' methods for disturbance events and final age detrending, and to allow the algorithm to detect and
#' remove growth suppression events in addition to growth releases. The default settings match the
#' parameters in Druckenbrod et al. (2024). A note of caution that this algorithm was built and
#' parameterized for release events - the behavior of suppression detrending has not been evaluated.
#'
#' @param data The collection of raw tree-ring series you wish to disturbance-detrend. Can be class
#' `"rwl"` (e.g., read in by \code{\link[dplR]{read.rwl}}) or `"data.frame"` (long format). If
#' `"data.frame"` you need at least these 3 columns: `"series"`, `"year"`, & `"rw"`.
#' @param win.len The window length, in years, for calculating percent growth change from moving
#' averages of raw ring width. Default is 15 years.
#' @param pgc.thresh Percent growth change (pgc) threshold for detection of disturbance events.
#' Default is 50%.
#' @param d.detrend.method Which method to use for detrending individual disturbance events. Options
#' are `"AgeDepSpline"`, `"Spline"`, `"Friedman"`, `"Hugershoff"`, `"Linear"`, and "`Mean"`. Default
#' is `"AgeDepSpline"`. See details for more information.
#' @param detrend.method Which method to use for the final age detrending on each series. Options are
#' all the methods available in \code{\link[dplR]{detrend.series}} except for "Ar". Default is
#' "AgeDepSpline".
#' @param nyrs Two numbers (in years) setting flexibility of splines for `d.detrend.method` and
#' `detrend.method`, respectively. For "AgeDepSpline", nyrs determines the initial flexibility. For
#' regular "Spline" the flexibility is constant throughout. Defaults are c(10, 25) years.
#' @param event.type The type of disturbance events to detrend. Options are "release", "suppression",
#' or "both". Default is "release".
#'
#' @details
#' The concept of separately detrending disturbances derives from Cook's (1987) linear aggregate
#' model. Although Cook's model distinguishes between endogenous and exogenous disturbances, in
#' practice this is not possible without prior knowledge. Instead, disturbances are defined
#' quantitatively.
#'
#' In \code{\link{d_detrend}} the percentage growth change (PGC) is calculated from
#' radial growth averaging (Nowacki & Abrams 1997) - i.e., the percentage change between adjacent (1
#' year of separation) `win.len` moving averages of ring width. Periods of PGC that are greater
#' than or equal to `abs(pgc.thresh)` are identified as disturbance events. Within the event
#' periods, the maximum PGC value is selected as the start year of the event. In this `modendro`
#' implementation, events detected within the earliest 2*(win.len-1) years of the original tree-ring
#' series are ignored. This is to avoid treating juvenile growth trends as disturbance events. Also
#' unique the `modendro` version, you can control which type of event is detected & detrended with
#' the `event.type` argument.
#'
#' Once events are identified, we power- or log-transform (see \code{\link{pwr_t_rwl}}) the whole
#' series based on the variance-mean (i.e., spread-level) relationship unique to each series, then
#' detrend them individually starting with the most recent event, moving backward in time until all
#' events are detrended. Event detrending involves extracting the subset of the transformed series
#' from the estimated start year (max PGC) to the end of the series. A line is then fit to the
#' subset, subtracted from the transformed values, and the 1st value of the original
#' power-transformed event subset is added back. Each event is then determined to be "sustained" or
#' "transient". If the detrended event values stay below/above the original values for the remainder
#' of the series, the event is sustained. If the detrended event values ever cross the original
#' values, the event is transient. Transient events are only detrended for the period between event
#' start and when they cross the original values.
#'
#' You can choose the type of line used to detrend each event through `d.detrend.method`. The
#' default `"AgeDepSpline"` is the recommended method. You can adjust the initial flexibility of the
#' spline through the 1st element of `nyrs`. You may want to adjust depending on your goals/system.
#' Other methods are available as well to allow further exploration of techniques. Methods `Spline`,
#' `Friedman`, `Linear`, and `Mean` are fairly straight-forward and follow the descriptions in
#' \code{\link[dplR]{detrend.series}}. Method `Hugershoff` is a custom implementation and follows
#' the "simple increment function" described by Warren & MacWilliam (1981), with the b parameter set
#'  to 1. The curve is fit via \code{\link[stats]{nls}}. If the fit fails, the event is fit with
#' `"AgeDepSpline"`.
#'
#' After all events are detrended, the long-term (whole series) age trend lines are fit to the
#' transformed series. Under the hood, this is just a call to \code{\link[dplR]{detrend.series}}.
#' Residual series are extracted here. You can elect to skip this age-detrending step by entering
#' `"none"` to the `detrend.method` argument, returning just the disturbance detrended series.
#'
#' Finally, the power transformation is reversed for the disturbance detrended series and the fitted
#'  age trend line, putting both in the original ring width units. Note that because the power
#' used was based on the empirical variance-mean relationship for the series, the variances of the
#' detrended events are also adjusted. In this way, the disturbance detrending is not a simple
#' subtraction. The last step is to compute the standardized ring width index by dividing the
#' reverse-transformed disturbance-detrended series by the reverse-transformed age trend.
#'
#' For more details, see Druckenbrod et al. (2024).
#'
#' I highly recommend you look at the results! Convenient visualizations of the
#' \code{\link{d_detrend}} process are available by using \code{\link{plot_d_detrend}}.
#'
#' Details of output data.frames:
#' # "PGC"
#'  \describe{
#'   \item{series}{Series ID}
#'   \item{year}{Year}
#'   \item{rw}{Ring widths}
#'   \item{pt.rw}{Transformed ring widths}
#'   \item{optimal.pwr}{Estimated optimal power of transformation}
#'   \item{action}{Type of transformation applied to series}
#'   \item{win.len}{Moving average window length (constant for all series)}
#'   \item{mav.Ra.M1}{Right-aligned moving average ring width}
#'   \item{mav.La.M2}{Left-aligned moving average ring width}
#'   \item{pgc}{Percent growth change}
#'   \item{pt.rw.ddtrd}{Transformed dist.-detrended ring widths}
#'   \item{message}{Information about ... Usually NA}
#'   \item{pt.rw.ddtrd.resid}{pt.rw.ddtrd - pt.rw.ddtrd.At}
#'   \item{pt.rw.ddtrd.At}{Fitted age trend (At) to transformed dist.-detrended ring widths}
#'   \item{detrend.method}{Method used to fit At}
#'   \item{rw.ddtrd}{Back-transformed dist.-detrended ring widths}
#'   \item{rw.ddtrd.At}{Back-transformed fitted age trend (At)}
#'   \item{Dt}{Back-transformed disturbacne signal}
#'   \item{rw.ddtrd.index}{rw.ddtrd / rw.ddtrd.At}
#'   ...
#' }
#'
#' # "Events"
#'  \describe{
#'   \item{series}{Series ID}
#'   \item{max.ind}{Index position of max percentage growth change (years from innermost ring)}
#'   \item{max.val}{Value of max (release) or min (suppression) percentage growth change}
#'   \item{pgc.thresh}{User-specified percent growth change threshold}
#'   \item{year}{Estimated first year of disturbance event}
#'   \item{message}{Summary of distubance events detected for the whole series}
#'   ...
#' }
#'
#' # "Dist. detrending"
#'  \describe{
#'   \item{series}{Series ID}
#'   \item{year}{Year}
#'   \item{pt.rw.i}{Transformed ring widths for this iteration}
#'   \item{curve}{Disturbance curve fit in this iteration}
#'   \item{pt.rw.ddtrd.i}{Transformed ring widths with disturbance curve removed in this iteration}
#'   \item{iter}{Iteration ID}
#'   \item{eventID}{Estimated first year of disturbance event}
#'   \item{event.type}{release or suppression}
#'   \item{event.dur}{sustained or transient}
#'   \item{d.detrend.method}{Curve-fitting method for disturbance events}
#'   ...
#' }
#'
#' @return A list containing 3 data.frames 1) the disturbance detrended ring width results including
#'  the percent growth change series and intermediate steps ("PGC"), 2) details about the specific
#'  disturbance events detected, and 3) details on the disturbance event trend fitting and removal
#'  process.
#'
#' All data.frames are in "long-format" with `series` as a categorical variable.
#'
#'
#' @references
#' Cook ER. 1987. The Decomposition of Tree-Ring Series for Environmental Studies.
#' \emph{Tree-Ring Bulletin}. \strong{47}:37–59
#'
#' Druckenbrod DL, Cook ER, Pederson N, Martin-Benito D. 2024. Detrending tree-ring widths in
#' closed-canopy forests for climate and disturbance history reconstructions.
#' \emph{Dendrochronologia}. \strong{85}:126195. https://doi.org/10.1016/j.dendro.2024.126195
#'
#' Druckenbrod DL, Pederson N, Rentch J, Cook ER. 2013. A comparison of times series
#' approaches for dendroecological reconstructions of past canopy disturbance events.
#' \emph{Forest Ecology and Management}, \strong{302}:23–33.
#' https://doi.org/10.1016/j.foreco.2013.03.040
#'
#' Rydval M, Druckenbrod D, Anchukaitis KJ, Wilson R. 2016. Detection and removal of
#' disturbance trends in tree-ring series for dendroclimatology.
#' \emph{Canadian Journal of Forest Research}, \strong{401}:387–401.
#' https://doi.org/10.1139/cjfr-2015-0366
#'
#' Rydval M et al. 2018. Influence of sampling and disturbance history on climatic sensitivity of
#' temperature limited conifers. \emph{The Holocene}, \strong{28}(10):1574-1587.
#' https://doi.org/10.1177/0959683618782605
#'
#' Melvin TM, Briffa KR, Nicolussi K, Grabner M. 2007. Time-varying-response smoothing.
#' \emph{Dendrochronologia}. \strong{25}(1):65–69. https://doi.org/10.1016/j.dendro.2007.01.004
#'
#' Nowacki GJ, Abrams MD. 1997. Radial-Growth Averaging Criteria for Reconstruction Disturbance
#' Histories from Presettlement-Origin Oaks. \emph{Ecological Monographs}. \strong{67}(2):225.
#' https://doi.org/10.2307/2963514
#'
#' Warren WG, MacWilliam SL. 1981. Test of a new method for removing the growth trend
#' from dendrochronological data.
#' \emph{Tree Ring Bulletin} \strong{41}:55–66.
#'
#' @import dplR
#' @import stats
#' @importFrom DescTools TukeyBiweight
#'
#' @export
#'
#' @seealso \code{\link{plot_d_detrend}}, \code{\link{ci_detect}}
#'
#' @examples
#'
#' # Missouri post oak example from Druckenbrod et al. 2024
#'
#' # Load Missouri post oak ring widths
#' data(mo024)
#'
#' # Run d_detrend - use the same settings as in Druckenbrod et al. 2024
#' mo024.ddtrd <- d_detrend(data = mo024,
#'                          win.len = 15,
#'                          pgc.thresh = 50,
#'                          d.detrend.method = "AgeDepSpline",
#'                          detrend.method = "AgeDepSpline",
#'                          nyrs = c(10, 30),
#'                          event.type = "release")
#'
#' # Generate plots using plot_d_detrend (output stored in a list)
#' mo024.ddtrd.plots <- plot_d_detrend(mo024.ddtrd)
#'
#' # look at plots for one series - note that there are 3 plots for each series - in RStudio click
#' # the back arrow in the plot viewer to see all of them
#' mo024.ddtrd.plots$DEM01C
#'
#' # Run again but with event.type = "both" to detrend both releases and suppressions
#' mo024.ddtrd.both <- d_detrend(data = mo024,
#'                               win.len = 15,
#'                               pgc.thresh = 50,
#'                               d.detrend.method = "AgeDepSpline",
#'                               detrend.method = "AgeDepSpline",
#'                               nyrs = c(10, 30),
#'                               event.type = "both")
#'
#' # Summary look at the 3 output data.frames in the list
#' sapply(mo024.ddtrd, head)
#'
#' # Generate plots using plot_d_detrend (output stored in a list)
#' mo024.ddtrd.both.plots <- plot_d_detrend(mo024.ddtrd.both)
#'
#' # look at plots for one series
#' mo024.ddtrd.both.plots$DEM01C

d_detrend <- function(data = NULL,
                      win.len = 15,
                      pgc.thresh = 50,
                      d.detrend.method = "AgeDepSpline",
                      detrend.method = "AgeDepSpline",
                      nyrs = c(10, 30),
                      event.type = "release") {
  #### Initial error catching
  ## data
  stopifnot(
    "data is not an object of class 'rwl', 'data.frame', or 'matrix'" =
      data.class(data) %in% "rwl" |
      data.class(data) %in% "data.frame" |
      data.class(data) %in% "matrix"
  )

  ## If data is an rwl, do some specific checks here
  if (any(class(data) %in% "rwl")) {
    stopifnot(
      "rwl has no rownames (must be years only) or no colnames (must be series IDs only)" =
        !is.null(rownames(data)) |
        !is.null(colnames(data))
    )

    if (apply(data, MARGIN = 2, FUN = \(x) all(is.na(x))) |> any() == TRUE) {
      these_are_NA <-
        colnames(data)[which(apply(data, MARGIN = 2, FUN = \(x) all(is.na(x))) == TRUE)]
      stop("The following series have no values (all NAs): " ,
           paste(these_are_NA, collapse = ", "))
    }

    # catch any negative values here and set to 0
    if (any(data < 0, na.rm = TRUE)) {
      data[which(data < 0), ] <- 0
      warning(
        "One or more negative values detected in the rwl file. These were replaced with 0s.
              Check your data to make sure this is correct
        (e.g., sometimes -9999 is used as a place-holder value for 0 rings)."
      )
    }

    ## Convert to long format and make sure there are no internal NAs in the series
    data <- modendro::rwl_longer(
      rwl = data,
      series.name = "series",
      dat.name = "rw",
      trim = TRUE,
      new.val.internal.na = 0
    )
  }

  ##
  # If the data and references are long-format data.frames, we need certain columns
  ##
  stopifnot(
    "data does not have a year column?" =
      any(colnames(data) %in% c("Year", "year")) == TRUE
  )

  # make sure that "year" columns are labelled as such
  colnames(data)[which((colnames(data) %in% c("Year", "year")) == T)] <- "year"

  # general check
  stopifnot(
    "data needs to contain these 3 columns: 'series', 'year', & 'rw'" =
      any(colnames(data) %in% "series") ||
      any(colnames(data) %in% "year") ||
      length(grep("rw", colnames(data))) > 0
  )

  # make sure that "rw" column is labelled as such
  colnames(data)[grep("rw", colnames(data))] <- "rw"

  # Make sure these columns are the right kind of data (year & rw need to be numeric)
  stopifnot(
    "year variable in data not numeric or integer" =
      is.numeric(data[, "year"]) ||
      is.integer(data[, "year"]) ||
      is.double(data[, "year"])
  )

  stopifnot(
    "ring width ('rw') variable in data not numeric or integer" =
      is.numeric(data[, "rw"]) ||
      is.integer(data[, "rw"]) ||
      is.double(data[, "rw"])
  )

  # Also need to make sure there are no internal NAs in rw or year columns (this will not be the
  # case if data started as an rwl, but might be if user supplies data)
  noNA.data <- lapply(split(data, f = data[,"series"]), FUN = \(this.series) {
    this.series[,"rw"] <- replace_internal_NAs(this.series[,"rw"], new.val = 0)
    this.series
  }) |> do.call(what = "rbind")

  if (any(noNA.data[,"rw"] != data[,"rw"])) {
    warning("\nSome series had internal NA values. These were replaced with 0s.\n")
  }

  data <- noNA.data

  ## win.len
  stopifnot(
    "win.len must be a numeric vector >= 5 & length = 1" =
      is.numeric(win.len) &&
      length(win.len) == 1 &&
      win.len >= 5
  )

  ## pgc.thresh
  stopifnot(
    "pgc.thresh must be numeric > 5 & length = 1" =
      is.numeric(pgc.thresh) &&
      length(pgc.thresh) == 1 &&
      pgc.thresh >= 5
  )

  ## d.detrend.method
  stopifnot(
    "d.detrend.method must be a character matching options these options: c(
        'AgeDepSpline',
        'Spline',
        'Hugershoff',
        'Linear',
        'Mean'
      )" =
      is.character(d.detrend.method) &&
      d.detrend.method %in% c(
        "AgeDepSpline",
        "Spline",
        "Hugershoff",
        "Linear",
        "Mean"
      )
  )

  ## detrend.method
  stopifnot(
    "detrend.method must match options (except 'Ar') in dplR::detrend.series(), or 'none'" =
      is.character(detrend.method) &&
      detrend.method %in% c(
        "Spline",
        "ModNegExp",
        "Mean",
        "Friedman",
        "ModHugershoff",
        "AgeDepSpline",
        "none",
        "None"
      )
  )

  ## nyrs
  stopifnot(
    "nyrs must contain positive numeric values & have length <= 2" =
      is.numeric(nyrs) &&
      length(nyrs) <= 2 &&
      all(nyrs > 0) == TRUE
  )

  ## event.type
  stopifnot(
    "event.type must be a character vector matching one of 'both','release', or 'suppression'" =
      is.character(event.type) &&
      length(event.type) == 1 &&
      event.type %in% c('both', 'release', 'suppression')
  )

  #### Power transform the ring width here
  orig.IDs <- unique(data[,"series"])
  # Find the optimal power of transformation and transform the series
  pt.list <- modendro::pwr_t_rwl(longer_rwl(data, series.name = "series", dat.name = "rw"))
  # The transformed series
  pt.long <- pt.list[["Transformed ring widths"]][, orig.IDs] |>
    modendro::rwl_longer(series.name = "series", dat.name = "pt.rw")

  # Info about transformations - will need this for back transforming
  mess.df <- pt.list[["Transformation metadata"]] |> do.call(what = "rbind")

  # Merge the ring widths with the power transformed widths & the transformation info
  series.long1 <- base::merge(data,
                              pt.long,
                              by = c("series", "year"),
                              all = TRUE)

  series.long <- base::merge(series.long1, mess.df, by = "series", all = TRUE)


  ## Radial growth averaging to identify disturbance events releases and suppressions
  # Percent growth change

  mov.avgs <- base::lapply(split(series.long, f = series.long$series), FUN = \(this.series) {
    # Make sure series are in temporal order
    this.series <- this.series[order(this.series[, "year"]), ]
    # Use cumsum() as a quick vectorized way to compute moving averages
    cs <- base::cumsum(c(0, this.series[, "rw"])) # Have to attach a 0 to beginning of series
    # nth differences where n = win.len
    result <- cs[(win.len + 1):length(cs)] - cs[1:(length(cs) - win.len)]
    # Moving averages are just the nth differences divided by n. Need right-aligned and left-aligned
    mav.Ra.M1 <- c(rep(NA, win.len - 1), result) / win.len # Right align
    mav.La.M2 <- c(result[2:length(result)], rep(NA, win.len)) / win.len # Left align, shifted 1

    # Calculate the % growth change according to Nowacki & Abrams 1997
    # Say the year is 1950 and win.len = 10, then M1 = ring width mean for 1941-1950 (preceeding) and
    # M2 = ring width mean for 1951-1960 (subsequent)
    # That is, M1 is the right-aligned mav, and M2 is the left-aligned mav, offset by one year
    #mov.avg <- c(rep(NA, (win.len - 1) / 2), result, rep(NA, (win.len - 1) / 2)) / win.len # Centered
    mov.avg.df <- base::data.frame(
      series = unique(this.series[, "series"]),
      year = this.series[, "year"],
      rw = this.series[, "rw"],
      pt.rw = this.series[, "pt.rw"],
      optimal.pwr = this.series[, "optimal.pwr"],
      action = this.series[, "action"],
      win.len = win.len,
      mav.Ra.M1 = mav.Ra.M1,
      # Right align
      mav.La.M2 = mav.La.M2 # Left align, shifted 1
    )

    # Compute the percent growth change
    mov.avg.df$pgc <- ((mov.avg.df$mav.La.M2 - mov.avg.df$mav.Ra.M1) / mov.avg.df$mav.Ra.M1) * 100

    ## ID the disturbances - use the user-defined threshold (pgc.thresh)

    # Find the beginning of disturbances by the year in which pgc >= pgc.thresh, but the year before
    # is < pgc.thresh. above.thresh is a vector of index values corresponding to mov.avg.df$pgc
    # Can use the discontinuities (1st diff > 1) here to delineate the events
    above.thresh <- base::which(abs(mov.avg.df$pgc) >= pgc.thresh)

    # It is possible that the first n PGC values are above/below the threshold due to the juvenile
    # growth trend. I think the best thing to do is to not consider the first win.len-1 of
    # the PCG values. If win.len = 15, the first 28 years will be skipped (the first 14 are NA,
    # so we skip those too).
    above.thresh <- above.thresh[!(above.thresh %in% 1:(2 * (win.len - 1)))]

    # Get the 1st above.thresh value and the rest
    dist.events.start <- above.thresh[c(1, which(diff(above.thresh) > 1) + 1)]
    dist.events.end <- above.thresh[c(which(diff(above.thresh) > 1), length(above.thresh))]
    # Combine the indices and find the max
    dist.events.max <- base::mapply(
      FUN = \(start, end) {
        dist.ind <- base::seq(from = start, to = end)
        df <- base::data.frame(
          max.ind = dist.ind[which.max(abs(mov.avg.df$pgc)[dist.ind])],
          # max.val is actually min.val for suppression events
          max.val = mov.avg.df$pgc[dist.ind][which.max(abs(mov.avg.df$pgc[dist.ind]))],
          pgc.thresh = pgc.thresh
        )
        df$year <- mov.avg.df$year[df$max.ind]
        df
      },
      start = dist.events.start,
      end = dist.events.end,
      SIMPLIFY = FALSE
    ) |>
      do.call(what = "rbind")

    if (base::is.null(dist.events.max)) {
      # Return a list with the whole series pgc and the (NULL) dist.events.max
      mov.avg.list <- base::list(mov.avg.df, dist.events.max)
    } else {
      # add series to the events index df
      dist.events.max$series <- this.series$series[1]
      # Give a message about what kind of events were detected and how many
      n.release <- nrow(dist.events.max[dist.events.max$max.val > 0, ])
      n.suppression <- nrow(dist.events.max[dist.events.max$max.val < 0, ])

      dist.events.max$message <- paste0(
        n.release,
        ifelse(n.release == 1, " release & ", " releases & "),
        n.suppression,
        ifelse(
          n.suppression == 1,
          " suppression detected",
          " suppressions detected"
        )
      )
      # Return a list with the whole series pgc and the max pgc values within each event
      mov.avg.list <- base::list(mov.avg.df, dist.events.max[, c("series",
                                                                 "max.ind",
                                                                 "max.val",
                                                                 "pgc.thresh",
                                                                 "year",
                                                                 "message")])
    }

    names(mov.avg.list) <- base::c("PGC", "Events")
    mov.avg.list

  }) # End of mov.avgs lapply()


  #### Use IDs of releases and suppressions to detrend each with ADS

  d.detrended <- base::lapply(mov.avgs, FUN = \(this.series.sublist) {
    # Separate the series data.frames from the sublist just once
    pgc_s  <- this.series.sublist[["PGC"]]
    evt_s  <- this.series.sublist[["Events"]]

    # determine if there were any events detected
    if (base::is.null(evt_s)) {
      event.ind <- data.frame(
        series = pgc_s[, "series"][1],
        max.ind = as.numeric(NA),
        max.val = as.numeric(NA),
        pgc.thresh = pgc.thresh,
        year = as.numeric(NA),
        message = "No disturbance events detected"
      )

      d.detrend.iter <- list(
        data.frame(
          series = pgc_s[, "series"][1],
          year = as.numeric(NA),
          pt.rw.i = as.numeric(NA),
          curve = as.numeric(NA),
          pt.rw.ddtrd.i = as.numeric(NA),
          iter = as.numeric(NA),
          eventID = as.numeric(NA),
          event.type = as.character(NA),
          event.dur = as.character(NA),
          d.detrend.method = as.character(NA)
        )
      )

      pt.rw.detrended <- pgc_s[, "pt.rw"]

    } else {
      # detrend the latest event first, moving backward in time
      # This has to be a for loop
      event.ind <- evt_s
      event.ind <- event.ind[base::order(event.ind$max.ind, decreasing = TRUE), ]

      # Subset event types based on user input to event.type argument
      if (event.type %in% "release") {
        event.ind <- event.ind[event.ind$max.val > 0, ]
      } else {
        if (event.type %in% "suppression") {
          event.ind <- event.ind[event.ind$max.val < 0, ]
        }
      }

      # In some cases, we might erase the only detected event in the steps above
      if (nrow(event.ind) == 0) {
        event.ind <- data.frame(
          series = pgc_s[, "series"][1],
          max.ind = as.numeric(NA),
          max.val = as.numeric(NA),
          pgc.thresh = pgc.thresh,
          year = as.numeric(NA),
          message = "No disturbance events detected"
        )

        d.detrend.iter <- list(
          data.frame(
            series = pgc_s[, "series"][1],
            year = as.numeric(NA),
            pt.rw.i = as.numeric(NA),
            curve = as.numeric(NA),
            pt.rw.ddtrd.i = as.numeric(NA),
            iter = as.numeric(NA),
            eventID = as.numeric(NA),
            event.type = as.character(NA),
            event.dur = as.character(NA),
            d.detrend.method = as.character(NA)
          )
        )

        pt.rw.detrended <- pgc_s[, "pt.rw"]

      } else {

        # Separate out the power transformed ring-width
        pt.rw.detrended <- pgc_s[, "pt.rw"]
        series.length <- base::length(pt.rw.detrended)

        ### Iterate through the events
        # Need a way to save each iteration - and the curves
        d.detrend.iter <- vector("list", nrow(event.ind))
        # Create a data.frame to hold results as we go
        iter.df <- data.frame(series = pgc_s[, "series"], year = pgc_s[, "year"])
        for (i in base::seq_along(event.ind$max.ind)) {
          # Get this event's end index - ie, the start index plus the win.len
          event.end <- event.ind$max.ind[i] + win.len
          # Subset the event - max pgc until the end of the series
          this.event <- pt.rw.detrended[event.ind$max.ind[i]:series.length]

          ## Fit the disturbance lines. Conceivably appropriate options are the ADS, loess spline,
          # & Hugershoff/Warren. Perhaps I should put straight line and mean in there too.
          # Suite of options could be c("AgeDepSpline", "Spline", "Friedman", "Hugershoff", "Linear", "Mean")

          # fit the ADS curve - from dplR
          if (d.detrend.method %in% "AgeDepSpline") {
            event.fit <- dplR::ads(y = this.event,
                                   nyrs0 = nyrs[1],
                                   pos.slope = TRUE)
          }

          if (d.detrend.method %in% "Spline") {
            event.fit <- dplR::caps(y = this.event,
                                    nyrs = nyrs[1],
                                    f = 0.5)
          }

          # if (d.detrend.method %in% "Friedman") {
          #   event.fit <- stats::supsmu(y = this.event,
          #                              x = base::seq_along(this.event),
          #                              periodic = FALSE)$y
          # }

          if (d.detrend.method %in% "Hugershoff") {
            # Hugershoff - fits to the detected period and the reminder of the series too
            # The formula is modified. t is an added parameter that controls how far above/below the
            # initial fit can go beyond the asymptote. b = 1, always, to allow t to work.
            # d mainly controls the asymptote value.

            hug_form <-
              formula(rwi ~ a * ((x - t) ^ 1) * exp(-c * (x - t)) + d)

            # Set up some start values & constraints for a, based on event.type (release or suppression)
            a_start <- ifelse(event.ind$max.ind[i] > 0, 0.1, -0.1)
            if (event.ind$max.ind[i] > 0) {
              a_const <- c(0.005, 5)
            } else {
              a_const <- c(-5, -0.005)
            }

            d.range <- range(this.event, na.rm = TRUE)

            lower_const <- list(a = a_const[1],
                                c = -5,
                                t = -10,
                                d = d.range[1])
            upper_const <- list(a = a_const[2],
                                c = 5,
                                t = 10,
                                d = d.range[2])

            event.df <- data.frame(rwi = this.event, x = base::seq_along(this.event))


            hug_fit.warnings <- list()
            hug_fit.errors <- list()

            hug_fit <- tryCatch(
              withCallingHandlers(
                nls(hug_form,
                    data = event.df,
                    start = list(a = a_start,
                                 c = 0.1,
                                 t = 1.5,
                                 d = DescTools::TukeyBiweight(this.event, na.rm = TRUE)),
                    algorithm = "port",
                    lower = lower_const,
                    upper = upper_const,
                    control = nls.control(
                      maxiter = 100,
                      minFactor = 1 / 4096,
                      warnOnly = FALSE
                    )
                ),
                warning = \(w) {
                  hug_fit.warnings[[length(hug_fit.warnings) + 1]] <<- conditionMessage(w)
                  invokeRestart("muffleWarning")
                }),
              error = \(e) {
                hug_fit.errors[[length(hug_fit.errors) + 1]] <<- conditionMessage(e)
                NULL
              }
            )

            # Print the warning and ID the series that had the issue
            if (length(hug_fit.errors) > 0) {

              cat("\nd.detrend.method Hugershoff fit failed for series '",
                  pgc_s$series[1], "', ", event.ind$year[i], " event", "\nRunning AgeDepSpline instead\n", hug_fit.errors[[1]], "\n\n",
                  sep = "")

              event.fit <- dplR::ads(y = this.event,
                                     nyrs0 = nyrs[1],
                                     pos.slope = TRUE)

              d.detrend.method <- "AgeDepSpline"

            } else {
              event.fit <- predict(hug_fit, newdata = event.df)
            }
          }

          if (d.detrend.method %in% "Linear") {
            event.fit <- stats::lm(this.event ~ base::seq_along(this.event)) |>
              stats::predict() |>
              base::as.numeric()
          }

          if (d.detrend.method %in% "Mean") {
            event.fit <- base::mean(this.event,
                                    na.rm = TRUE)
          }


          # Detrend by subtracting the curve out and adding back the 1st value of the pt series
          new.detrended <- pt.rw.detrended
          new.detrended[event.ind$max.ind[i]:series.length] <- (this.event - event.fit) + this.event[1]

          # Transfer the starting series to the data.frame
          iter.df$pt.rw.i <- pt.rw.detrended


          # Druckenbrod's text on distinguishing transient vs  sustained events:
          # "...sustained release events occur when disturbance detrended ring-width values
          # remain less than initial ringwidth values for the remainder of the series after
          # the release event. Conversely, transient events occur when detrended values only
          # decrease initial ring-width values for a certain sequence of years after the release
          # event instead of the entire remaining series. Transient release events terminate when
          # the detrended value would be greater than the initial ring-width value, indicating
          # that a tree likely returned to suppressed growth conditions." (CTM - notice the implied
          # underlying condition of suppression).

          # "After the event" I think must mean anything past event.end or up to event.end if
          # event.end = series.length

          # The MatLab code looks like this:
          # find(DSti(lngthr(mav:end)) < DStT(lngthr(mav:end)),1) + mav-1

          # DSti - initial (original) ring widths of the disturbance
          # DStT - detrended ring widths of the disturbance

          # In R, this means:
          # which(DSti[lngthr[mav:length(lngthr)]] < DStT[lngthr[mav:length(lngthr)]])[1] + mav - 1
          # The 1st ring width value where the original ring width is less than the detrended value +
          # the moving window length - 1

          # I suppose if the event is a suppression, "less than" becomes "more than"

          # I interpret this (and the other code) to mean that you fit the whole series, to the end and detrend it.
          # Then decide whether or not that was justified or not based on the above description
          # of transient or sustained. If justified, you keep the detrend along the whole curve until the end of
          # the series. If not, you just keep the detrended section during the identified disturbance period.

          # Determine if event is sustained or transient
          # The sus object is an index position in the original series - the place where the detrended
          # values rise above or dip below the original value. event.end is added to translate back
          # from event index values to original series index values.
          if (event.ind$max.val[i] > 0) {
            # For releases...
            sus <- which(new.detrended[event.end:series.length] >
                           pt.rw.detrended[event.end:series.length])[1] + event.end - 1
          } else {
            # for suppressions...
            sus <- which(new.detrended[event.end:series.length] <
                           pt.rw.detrended[event.end:series.length])[1] + event.end - 1
          }


          if (is.na(sus)) {
            # if sus is NA, then the release extended to the end of the series (sustained)
            # Detrend the whole rest of series in this case
            pt.rw.detrended[event.ind$max.ind[i]:series.length] <- new.detrended[event.ind$max.ind[i]:series.length]
            event.dur <- "sustained"
          } else {
            # transient
            # Detrend just the appropriate section
            pt.rw.detrended[event.ind$max.ind[i]:sus] <- new.detrended[event.ind$max.ind[i]:sus]
            event.dur <- "transient"
          }

          iter.df$curve <- iter.df$pt.rw.i - pt.rw.detrended
          iter.df$curve <- ifelse(iter.df$curve == 0, NA, iter.df$curve + this.event[1])
          iter.df$pt.rw.ddtrd.i <- pt.rw.detrended
          iter.df$iter <- i
          iter.df$eventID <- event.ind$year[i]
          iter.df$event.type <- ifelse(event.ind$max.val[i] > 0, "release", "suppression")
          iter.df$event.dur <- event.dur
          iter.df$d.detrend.method <- d.detrend.method

          d.detrend.iter[[i]] <- iter.df
        }
      }
    } ## End of disturbance event detrending


    ## Rbind the disturbance detrending iterations to a data.frame
    ddtrd <- do.call(what = "rbind", d.detrend.iter)

    # ggplot(ddtrd) +
    #   geom_hline(yintercept = 0, color = "red") +
    #   geom_line(aes(year, pt.rw.i), color. = "grey") +
    #   geom_line(aes(year, pt.rw.ddtrd.i)) +
    #   geom_line(aes(year, curve), color = "blue") +
    #   facet_wrap(~iter, ncol = 1, strip.position = "right")

    # Transfer the disturbance detrended series to the data.frame
    pgc_s[, "pt.rw.ddtrd"] <- pt.rw.detrended

    # Make a message column to carry warnings
    pgc_s[, "message"] <- NA

    ## Option to not detrend
    if (detrend.method %in% c("none", "None")) {
      # Add the detrend.method and stop here
      pgc_s[, "detrend.method"] <- detrend.method


      ## Reverse the power transformation
      if (pgc_s[, "action"][1] %in% "log10 transformed") {
        pgc_s[, "rw.ddtrd"] <- 10^pgc_s[, "pt.rw.ddtrd"]
      } else {
        if (pgc_s[, "action"][1] %in% "Power transformed") {

          # Make sure that there are no negative values in the transformed series
          if (any(pgc_s[, "pt.rw.ddtrd"] < 0)) {
            warning(paste("Negative power transformed ring-width value(s) set to 0.001 in series",
                          pgc_s[, "series"][1]))

            pgc_s[, "pt.rw.ddtrd"][pgc_s[, "pt.rw.ddtrd"] < 0] <- 0.001
            pgc_s[, "message"] <- "Negative pt.rw.ddtrd value(s) set to 0.001"
          }

          pgc_s[, "rw.ddtrd"] <-
            pgc_s[, "pt.rw.ddtrd"]^(1 / as.numeric(pgc_s[, "optimal.pwr"][1]))
        } else {
          # No transformation
          pgc_s[, "rw.ddtrd"] <- pgc_s[, "pt.rw.ddtrd"]
        }
      }

      ## Calculate the disturbance signal/ disturbance index
      pgc_s[, "Dt"] <- pgc_s[, "rw"] - pgc_s[, "rw.ddtrd"]

    } else {
      ## Fit age trend to the whole series after removing all disturbances

      # Run detrend.series(). Note that the standardize arg passes
      # to the difference arg here. standardize = TRUE means divide out the curve, FALSE
      # means take the difference. Hence we feed the opposite to the difference arg.
      # In the case of "Ar" method, override the standardize arg - results are weird otherwise.

      ## Catch warnings and print the series name before - or edit the warning message otherwise.
      detrend.series.warnings <- list()

      dtrd.result <- withCallingHandlers(
        dplR::detrend.series(
          pgc_s[, "pt.rw.ddtrd"],
          y.name = pgc_s[, "series"][1],
          method = detrend.method,
          difference = TRUE,
          nyrs = nyrs[2],
          make.plot = FALSE,
          pos.slope = TRUE,
          return.info = TRUE
        ),
        warning = \(w) {
          detrend.series.warnings[[length(detrend.series.warnings) + 1]] <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        })

      # Print the warning and ID the series that had the issue
      if (length(detrend.series.warnings) > 0) {
        cat("\nseries '", pgc_s$series[1], "' dplR::detrend.series() warning:\n",
            detrend.series.warnings[[1]], "\n\n", sep = "")
      }

      # Extract the detrended residuals and the fitted trend curve
      pgc_s[, "pt.rw.ddtrd.resid"] <- dtrd.result$series
      pgc_s[, "pt.rw.ddtrd.At"] <- dtrd.result$curves

      # Extract the detrend method actually used.
      # This accounts for the stepped methods in ModNegExp and ModHugershoff.
      pgc_s[, "detrend.method"] <- dtrd.result[["model.info"]][[detrend.method]]$method


      ## Reverse the power transformation - for both the disturbance-detrended series and the
      # age trend curve
      if (pgc_s[, "action"][1] %in% "log10 transformed") {
        pgc_s[, "rw.ddtrd"] <- 10^pgc_s[, "pt.rw.ddtrd"]
        pgc_s[, "rw.ddtrd.At"] <- 10^pgc_s[, "pt.rw.ddtrd.At"]
      } else {
        if (pgc_s[, "action"][1] %in% "Power transformed") {

          # Make sure that there are no negative values in the transformed series
          ### This can create issues when there was no ddtrd but there were low values.
          if (any(pgc_s[, "pt.rw.ddtrd"] < 0)) {
            warning(paste("Negative power transformed ring-width value(s) set to 0.001 in series",
                          pgc_s[, "series"][1]))

            pgc_s[, "pt.rw.ddtrd"][pgc_s[, "pt.rw.ddtrd"] < 0] <- 0.001
            pgc_s[, "pt.rw.ddtrd.At"][pgc_s[, "pt.rw.ddtrd.At"] < 0] <- 0.001
            pgc_s[, "message"] <- "Negative pt.rw.ddtrd value(s) set to 0.001"
          }

          pgc_s[, "rw.ddtrd"] <-
            pgc_s[, "pt.rw.ddtrd"]^(1 / as.numeric(pgc_s[, "optimal.pwr"][1]))
          pgc_s[, "rw.ddtrd.At"] <-
            pgc_s[, "pt.rw.ddtrd.At"]^(1 / as.numeric(pgc_s[, "optimal.pwr"][1]))
        } else {
          # No transformation
          pgc_s[, "rw.ddtrd"] <- pgc_s[, "pt.rw.ddtrd"]
          pgc_s[, "rw.ddtrd.At"] <- pgc_s[, "pt.rw.ddtrd.At"]
        }
      }

      ## Calculate the disturbance signal/ disturbance index
      pgc_s[, "Dt"] <- pgc_s[, "rw"] - pgc_s[, "rw.ddtrd"]

      ## Get standardized index from dividing the At curves from the disturbance detrended series
      pgc_s[, "rw.ddtrd.index"] <- pgc_s[, "rw.ddtrd"] / pgc_s[, "rw.ddtrd.At"]
    }

    # Reassemble the series sublist
    this.series.sublist <- list(pgc_s, event.ind, ddtrd)
    names(this.series.sublist) <- c("PGC", "Events", "Dist. detrending")
    this.series.sublist

  }) # End of d.detrended lapply()

  #### Output as list of 3 data.frames (all series within)
  out.list <- list(do.call("rbind", lapply(d.detrended, FUN = \(x) {
    x[["PGC"]]
  })),
  do.call("rbind", lapply(d.detrended, FUN = \(x) {
    x[["Events"]]
  })),
  do.call("rbind", lapply(d.detrended, FUN = \(x) {
    x[["Dist. detrending"]]
  })))

  names(out.list) <- c("PGC", "Events", "Dist. detrending")

  class(out.list) <- c("list", "d_detrend")

  return(out.list)

} ## End of d_detrend
