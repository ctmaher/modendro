#' Statistical cross-dating checking using methods inspired by CooRecorder/CDendro and COFECHA.
#'
#'
#' @description
#' This function implements a custom statistical check on the dating of a collection of tree-ring
#' series. It is inspired by operations in CDendro (Larsson & Larsson 2024) and COFECHA,
#' but does not do all of what either of those programs can do. \code{\link{xd_check}} does simple
#' whole-series cross-correlations and checks for potential bad segments between leave-one-out (LOO)
#' chronologies and 0-n reference chronologies. The power of this method comes from the multiple
#' comparisons - more references = more confidence, assuming your references are appropriate (e.g.,
#' independent, properly cross-dated, same species, same region).
#'
#' You provide collections of raw tree-ring widths (not pre-made chronologies), and xd_check will
#' standardize/normalize all series and then generate the chronologies internally. If you only
#' provide the `data` collection, you get only LOO comparisons.
#'
#' The main output of the function is a summary based on how many references have best and
#' significant correlations with the series as dated and are ranked lowest, medium, and highest
#' confidence.
#'
#' One potential use for \code{\link{xd_check}} is to generate "clean" references for subsequent
#' runs of the function. You can do this by supplying the reference collection to the `data`
#' argument and not supplying any of the reference arguments (i.e., LOO-only mode). Set the other
#' arguments as you please and then take the "Highest" subset as one of the `refx` arguments in the
#' next run. Note that this assumes that you have plenty of good series in your collection!
#'
#' A note of awareness: statistical checks are correlations and can be incorrect. It is not uncommon
#' for forest interior series to have no correlation with any reference. In these cases, rely on
#' ground-truths like known sample date and presence of bark.
#'
#' Differences from CDendro: the chronology building process in \code{\link{xd_check}} is different than
#' the default in CDendro, even though the p2yrsL standardization/normalization derives from
#' CDendro. The default in CDendro is to create a mean value chronology of "Heavy" detrended series
#' (the user can select the type of detrending), then apply the p2yrsL standardization on the
#' chronology. The \code{\link{xd_check}} method is to apply p2yrsL to all the raw ring width series
#' , then take the mean of those to create the chronology. This is the simplest approach, as p2yrsL
#' is already akin to an aggressive detrending method. Becuase of these differences, the
#' correlations you see in CDendro and CooRecorder might be quantitatively different from what
#' \code{\link{xd_check}} produces.
#'
#'
#' @param data The collection of raw tree-ring series you wish to check. Can be class "rwl" or
#' "data.frame"/"data.table (long-format). This is required. If long-format you need 3 columns:
#' "series", "year", & "rw".
#' @param ref A single reference collection or a list of reference collections of raw tree-ring
#' series (NOT a chronology!). Can be rwl- or long-format. Optional, but more
#' references are better here. If long-format you need 3 columns: "series", "year", & "rw". If a
#' list, can be a named list so that you can give unique names to your references.
#' @param std.method The method of standardization/normalization to apply to the series and
#' references. Default (and currently only) option is "p2yrsL", same default as CooRecorder/CDendro.
#' @param max.offset The maximum offset (in years) you want to check in either direction from the
#' years as dated. Minimum is 5 years, default is 10 years. No true maximum limit, but overlap
#' between series and references must be at least 15 years.
#' @param p.thresh The p-value threshold to determine if a correlation is significant. Default is
#' the common 0.05.
#' @param out.format The format of subsets of data based on the dating confidence ranking. Options
#' are "rwl" or "long".
#' @param seg.penalty Sets the cost of increasingly large segment offsets
#' @param min.seg The shortest potential bad segment reported
#'
#' @details
#' \code{\link{xd_check}} does a standardization/normalization on all individual series within the
#' `data` and `refn` collections, then it computes chronologies using Tukey's Biweight Robust mean
#' (via \code{\link[DescTools]{TukeyBiweight}}). The LOO method is common in dendrochronology as
#' the method to compute interseries correlations for a tree-ring collection. In LOO, each series
#' has its own chronology made up of the remaining series.
#'
#' The bad segment checking is based on continuous (annual) contribution to the correlation
#' coefficient between a series and each reference chronology, instead of moving window correlations
#' of fixed length. In a correctly-dated (and correlated) series, the cumulative correlation
#' coefficient contributions are positive and increasing going backward in time. The function
#' attempts to find the "optimal offset path" using the Viterbi algorithm, and suggests the segments
#' and offsets that might be incorrectly dated.
#'
#' The basic ranking of dating confidence here is based on two basic criteria: 1) is the
#' whole-series correlation best as dated? and 2) is the correlation significant as dated? These
#' criteria are assessed for each reference comparison. For the highest confidence group, potential
#' bad segments are not allowed. The more references you have, the higher your overall confidence
#' should be.
#'
#' @return A list containing 1) the partitioned results based on relative confidence (itself a
#' 4-element list), 2) the raw correlation results in a nested list for easy viewing in RStudio,
#'  3) subsets of the data collection corresponding to the relative confidence (format specified
#' with `out.format`), and 4) a data.table of details on potentially bad segments.
#'
#' @references
#' Larsson & Larsson (2024) \emph{CDendro and CooRecorder programs of the CDendro package},
#'  Cybis Elektronik & Data AB. https://www.cybis.se/forfun/dendro/index.htm
#'
#' @importFrom DescTools TukeyBiweight
#'
#' @export
#'
#' @examples
#' # Examples to come


xd_check <- function(data = NULL, # the data you are checking. long format or rwl
                     ref = NULL, # a single reference data.frame or a list of references.
                     std.method = "p2yrsL", # standardization method
                     max.offset = 10, # max offset in years to check
                     p.thresh = 0.05,
                     out.format = "rwl",
                     seg.penalty = 4, # cost per year of offset change in segment path
                     min.seg = 20 # shortest segment reported (years)
) {
  . <- j_variable <- NULL
  ## Error catching
  # The function cannot do anything if data is NULL
  stopifnot(
    "You must provide tree-ring data to the data argument" =
      !is.null(data)
  )

  #### data format detection and conversion ####
  # Also does some quality checking. Leverages rwl_longer() and longer_rwl()
  data <- rwl_dt_conv(data = data,
                      series.name = "series",
                      dat.name = "rw",
                      new.val.internal.na = 0)

  ## Detailed checking of the data and references to determine if they are rwl format or
  # long format. If rwls are read in using dplR::read.rwl they will have class "rwl" already.
  # What else can signify this? Not sure of a robust test
  # Start with the simple test. For this function, I will convert rwls to long format.

  # The references should only go though these checks if they exist.
  # For the list case, it might make the most since to check through everything, then check through
  # the results once to give the entire result at once (maybe).
  if (!is.null(ref)) { # If there are any references...
    if (is.list(ref) &&
        !(is_rwl_like(ref) ||
          data.class(ref) %in% "data.table" ||
          data.class(ref) %in% "data.frame")) { # If ref is a list
      ref.names <- names(ref) # keep the names for later
      ref <- lapply(seq_along(ref), FUN = \(i) {

        this.ref <- with_context(rwl_dt_conv(data = ref[[i]],
                                             series.name = "series",
                                             dat.name = "rw",
                                             new.val.internal.na = 0), if (is.null(names(ref)[i])) {
                                               sprintf("ref element '", names(ref)[i], "'")
                                             } else {
                                               sprintf("ref element %d", i)
                                             }
        )

        this.ref
      })

      names(ref) <- ref.names

    } else { # if ref is not a list

      ref <- rwl_dt_conv(data = ref,
                         series.name = "series",
                         dat.name = "rw",
                         new.val.internal.na = 0)
    }

    ### Check that ref has multiple series - warn if there are less than 10 unique series?

  }

  ###

  ## The other args
  # std.method arg
  stopifnot(
    "std.method must be one of 'p2yrsL', 'AR', or 'ARIMA'" =
      is.character(std.method) &
      any(std.method %in% c("p2yrsL","AR","ARIMA"))
  )

  # max.offset arg
  stopifnot(
    "max.offset must be a postive integer" =
      is.numeric(max.offset) &
      max.offset > 0
  )

  # p.thresh arg
  stopifnot(
    "p.thresh must be a postive number < 1" =
      is.numeric(p.thresh) &
      p.thresh > 0 &
      p.thresh < 1
  )

  # out.format arg
  stopifnot(
    "out.format options are c('rwl', 'long')" =
      is.character(out.format) &&
      length(out.format) == 1 &&
      out.format %in% c("rwl", "long")
  )

  # segment-related args
  stopifnot(
    "seg.penalty must be a positive number" =
      is.numeric(seg.penalty) && length(seg.penalty) == 1 && seg.penalty > 0,
    "min.seg must be a positive integer" =
      is.numeric(min.seg) && length(min.seg) == 1 && min.seg > 0
  )


  ### Start the process
  # Steps: 1) convert to p2yrsL, AR, or ARIMA. 2) Build chronologies. 3) Run correlations.
  # 4) Assemble the results

  # Make a list of the data and the ref(s)
  if (is.list(ref) && # If ref is already a list
      !(is_rwl_like(ref) ||
        data.class(ref) %in% "data.table" ||
        data.class(ref) %in% "data.frame")) {

    dat.ref.list <- c(list(data), ref)

    if(is.null(names(ref))) { # NULL names
      length.ref <- ifelse(is.null(ref), 0, length(ref))
      if (length.ref == 0) { # ref is NULL altogether
        names(dat.ref.list) <- "data"
      } else {
        names(dat.ref.list) <- c("data", paste0("ref", 1:length.ref))
      }
    } else { # if ref has names
      names(dat.ref.list) <- c("data", names(ref))
    }

  } else { # ref is not a list

    dat.ref.list <- list(data, ref)
    names(dat.ref.list) <- c("data", "ref1")
  }

  # Remove the null values from the working list
  dat.ref.list <- dat.ref.list[!unlist(lapply(dat.ref.list, FUN = \(x)is.null(x)))]

  # All data is in long-format data.tables, within a list

  # Apply the standardization to each element of the list
  std.all <- lapply(dat.ref.list, FUN = \(x) {
    if (std.method %in% "p2yrsL") {
      std.df <- p2yrsL(rwl = x, limit = 2.6, out.format = "long")
      data.table::setnames(std.df, old = "p2yrsL", new = "std.series")
    }

    if (std.method %in% "AR") {
      # std.df <- longer_rwl(data = x, series.name = "series", dat.name = "rw") |>
      #   p2yrsL() |>
      #   rwl_longer()
      # data.table::setnames(std.df, old = "ar.resid", new = "std.series")
    }

    if (std.method %in% "ARIMA") {
      # std.df <- longer_rwl(data = x, series.name = "series", dat.name = "rw") |>
      #   p2yrsL() |>
      #   rwl_longer()
      # data.table::setnames(std.df, old = "arima.resid", new = "std.series")
    }
    std.df
  })


  ## Set up the sliding windows (i.e., lags) based on the provided max.offset
  lags <- seq(-max.offset, max.offset, by = 1)
  # These then have to be applied to each series + chron combination

  # set min.overlap - the way I'm handling this is causing problems I think
  min.overlap <- 15

  # Do the leave-one-out chronology method on the data (n chrons = n series)
  data.series.split <- split(std.all[["data"]], f = std.all[["data"]]$series)

  data.loo <- lapply(data.series.split,
                     FUN = \(this.series) {

                       #series$std.series <- ifelse(is.na(series$std.series), 0, series$std.series)

                       this.seriesID <- unique(this.series$series)
                       loo.chron <- std.all[["data"]][!(std.all[["data"]]$series
                                                        %in%
                                                          this.seriesID),
                                                      .(std.chron =
                                                          DescTools::TukeyBiweight(std.series, na.rm = TRUE),
                                                        samp.depth = .N),
                                                      by = year]
                       setorder(loo.chron, year)

                       loo.chron <- loo.chron[samp.depth >= 10,] # enforce minimum sample depth

                       # Get the possible lag years for the series
                       lags.df <- lapply(lags, FUN = \(n.lag) {
                         this.lag <- data.table::data.table(this.series$year + n.lag)
                         setnames(this.lag, old = "V1", new = paste0("lag", n.lag, ".year"))
                         this.lag
                       }) |> data.table::cbindlist()

                       # Merge the series with the respective loo chronology with each lag
                       merged.loo.lags <- apply(lags.df, MARGIN = 2, FUN = \(this.lag) {
                         this.lag.merged <- merge(cbind(this.series[, c("series", "std.series"), drop = FALSE],
                                                        this.lag),
                                                  loo.chron,
                                                  by.x = "this.lag", by.y = "year")

                         data.table::setnames(this.lag.merged, old = "this.lag", new = "year")

                         # Calculate the z-scores and correlation coef contribution (Pearson)
                         this.lag.merged[, cor.contr := as.numeric(scale(std.series)) *
                                           as.numeric(scale(std.chron))]
                         this.lag.merged
                       })

                       # seg_check(series.dt = this.series, chron.dt = loo.chron, lags = seq(-1, 1))
                       # plot_seg_check(this.series, loo.chron, lags = seq(-2, 2))
                       #
                       # Set min overlap
                       lag.dfs <- lapply(merged.loo.lags, FUN = \(this.lag) {
                         if (nrow(this.lag) < min.overlap) {
                           NA
                           #this.lag[0, , drop = FALSE]
                         } else {
                           this.lag
                         }
                       })

                       lag.dfs <- lag.dfs[!is.na(lag.dfs)]

                       # get lag identifiers
                       lag.dfs <- mapply(FUN = \(l, lnames) {
                         l[, lag := gsub(".*lag([^.]+)\\.year.*", "\\1", lnames)]
                         l
                       }, l = lag.dfs, lnames = names(lag.dfs), SIMPLIFY = FALSE)

                       lag.dfs

                     })


  # Make single chronologies for the rest - if they exist
  if (length(std.all) > 1) {
    # This will be a list of length 1 to 3 (or n refs)
    std.chrons <- lapply(std.all[!(names(std.all) %in% "data")],
                         FUN = \(this.ref) {
                           std.chron <- this.ref[, .(std.chron = DescTools::TukeyBiweight(std.series, na.rm = TRUE),
                                                     samp.depth = .N),
                                                 by = year]
                           setorder(std.chron, year)
                           std.chron[samp.depth >= 10,] # enforce minimum sample depth
                         })

    data.refs <- lapply(std.chrons, FUN = \(this.chron) {

      this.ref <- lapply(data.series.split,
                         FUN = \(this.series) {

                           # Get the possible lag years
                           lags.df <- lapply(lags, FUN = \(n.lag) {
                             this.lag <- data.table::data.table(this.series$year + n.lag)
                             setnames(this.lag, old = "V1", new = paste0("lag", n.lag, ".year"))
                             this.lag
                           }) |> data.table::cbindlist()

                           # Merge the series with the respective loo chronology with each lag
                           merged.lags <- apply(lags.df, MARGIN = 2, FUN = \(this.lag) {
                             this.lag.merged <- merge(cbind(this.series[, c("series", "std.series"), drop = FALSE],
                                                            this.lag),
                                                      this.chron,
                                                      by.x = "this.lag", by.y = "year")

                             data.table::setnames(this.lag.merged, old = "this.lag", new = "year")

                             # Calculate the z-scores and correlation coef contribution (Pearson)
                             this.lag.merged[, cor.contr := as.numeric(scale(std.series)) *
                                               as.numeric(scale(std.chron))]
                             this.lag.merged
                           })

                           # Set min overlap
                           lag.dfs <- lapply(merged.lags, FUN = \(this.lag) {
                             if (nrow(this.lag) < min.overlap) {
                               NA
                             } else {
                               this.lag
                             }
                           })

                           lag.dfs <- lag.dfs[!is.na(lag.dfs)]

                           # get lag identifiers
                           lag.dfs <- mapply(FUN = \(l, lnames) {
                             l[, lag := gsub(".*lag([^.]+)\\.year.*", "\\1", lnames)]
                             l
                           }, l = lag.dfs, lnames = names(lag.dfs), SIMPLIFY = FALSE)

                           lag.dfs
                         })

    })

    data.refs.list <- c(list(data.loo), data.refs)
    names(data.refs.list) <- c("LOO", names(std.chrons))

  } else {
    data.refs.list <- list(data.loo)
    names(data.refs.list) <- "LOO"
  }

  n.refs <- length(data.refs.list)   # move up from ~l.503
  seg.res <- data.table::rbindlist(
    lapply(data.refs.list, \(this.ref) {
      data.table::rbindlist(lapply(this.ref, seg_path,
                                   penalty = seg.penalty,
                                   min.seg = min.seg,
                                   p.thresh = p.thresh))
    }),
    idcol = "reference")
  seg.res <- seg_agree(seg.res, n.refs = n.refs)

  ## Run the overall correlations between each series and each reference
  # The output here is a single data.frame containing all the results - can make any other type of
  # organization from this point
  cor.res.raw <- mapply(FUN = \(this.combo, these.names) {
    mapply(FUN = \(this.series, series.name) {
      if (length(this.series) == 0) { # Handle really short series properly (no lags possible)
        all.lags <- data.table::data.table(series = series.name,
                                           reference = these.names,
                                           offset = 0,
                                           overlap = NA,
                                           cor.coef = NA,
                                           p.val = NA,
                                           T.val = NA)
      } else {
        all.lags <- lapply(this.series, FUN = \(this.lag) { # Do correlations on this layer
          overlap <- nrow(this.lag)
          if (overlap < min.overlap) {
            cor.coef <- NA
            p.val <- NA
          } else {
            suppressWarnings(
              this.cor.res <- cor.test(this.lag$std.series,
                                       this.lag$std.chron,
                                       method = "spearman")
            )
            # Round the results for easy viewing
            cor.coef <- round(this.cor.res$estimate[[1]], digits = 2)
            p.val <- round(this.cor.res$p.value[[1]], digits = 3)
          }

          res.df <- data.table::data.table(series = this.lag$series[1],
                                           reference = these.names,
                                           offset = this.lag$lag[1],
                                           overlap = overlap,
                                           cor.coef = cor.coef,
                                           p.val = p.val)
          res.df[, T.val := round(
            abs(res.df$cor.coef) *
              sqrt((res.df$overlap-2)/(1-(res.df$cor.coef)^2)),
            digits = 2)]

          res.df
        }) |> data.table::rbindlist()
      }

      all.lags

    }, this.series = this.combo, series.name = names(this.combo),
    SIMPLIFY = FALSE) |> data.table::rbindlist()

  }, this.combo = data.refs.list, these.names = names(data.refs.list),
  SIMPLIFY = FALSE) |> data.table::rbindlist()


  # Calculate the sugg.ID and sugg.OD using the original series data
  data.ID.OD <- data[, .(sugg.ID = min(year), sugg.OD = max(year)), by = series]

  cor.res.raw1 <- cor.res.raw[data.ID.OD, on = "series"]

  cor.res.raw1[, `:=`(sugg.ID = sugg.ID + as.numeric(offset),
                      sugg.OD = sugg.OD + as.numeric(offset))]

  # Reorganize the columns
  cor.res.raw1 <- cor.res.raw1[, c("series","reference","offset","sugg.ID","sugg.OD",
                                   "overlap","cor.coef","p.val","T.val")]
  # Add actual overlap to series that didn't have enough overlap for correlations
  # cor.res.raw1$overlap <- ifelse(is.na(cor.res.raw1$overlap),
  #                                cor.res.raw1$sugg.OD - cor.res.raw1$sugg.ID,
  #                                cor.res.raw1$overlap)

  ## Sort/filter/arrange this massive amount of data into a useful output
  # First step is to reorganize into a list that can be perused easily by the user
  cor.res <- lapply(split(cor.res.raw1, f = cor.res.raw1$series), FUN = \(this.series) {
    lapply(split(this.series, f = this.series$reference), FUN = \(this.combo) {
      data.table::setorder(this.combo, -cor.coef, na.last = TRUE)
      this.combo#[1:5,]
    })
  })

  # The next step is applying some decision-making to the results so that we can turn a manual
  # decision process into an automated one. Ideal output is to distill results into 1 row per
  # series.

  # Perhaps the best way (consider if all the references have best as dated and are significant)
  #n.refs <- length(data.refs.list)
  all.bad <- paste("0", n.refs, sep = "/")

  cor.res.message <- lapply(split(cor.res.raw1, f = cor.res.raw1$series), FUN = \(this.series) {

    this.series.sum <- lapply(split(this.series, f = this.series$reference),
                              FUN = \(this.combo) {
                                data.table::setorder(this.combo, -cor.coef, na.last = TRUE)
                                #rownames(this.combo) <- 1:nrow(this.combo)
                                as.dated <- this.combo[offset %in% "0",]
                                if (!nrow(as.dated)) {
                                  return(data.table::data.table(series = this.combo$series[1],
                                                                reference = this.combo$reference[1],
                                                                best.as.dated = 0, sig.as.dated = 0))
                                } else {

                                  data.table::data.table(series = as.dated$series,
                                                         reference = this.combo$reference[1],
                                                         best.as.dated = data.table::fifelse(which(this.combo$offset %in% "0")
                                                                                             %in% "1" &
                                                                                               nrow(this.combo) > 1, 1, 0),
                                                         sig.as.dated = data.table::fifelse(as.dated$p.val > p.thresh |
                                                                                              is.na(as.dated$p.val), 0, 1))
                                }
                              }) |> data.table::rbindlist()

    # Make sure it says 0/n if the overlap is inadequate to have anything

    # I need to figure out which of the references (including LOO) have correlations for
    # this.series as dated.
    these.refs <- this.series[, .(ncorrs = length(na.omit(cor.coef))), by = reference]

    agg.df <- this.series.sum[, .(best.as.dated = paste(sum(best.as.dated), n.refs, sep = "/"),
                                  sig.as.dated = paste(sum(sig.as.dated), n.refs, sep = "/")),
                              by = series]

    agg.df[, `:=`(as.dated.ID = this.series[offset %in% "0",]$sugg.ID[1],
                  as.dated.OD = this.series[offset %in% "0",]$sugg.OD[1])]

    ## Add the correlation coefficients for as dated
    # Which refs do we have correlation coefs for?
    # cor.coefs <- this.series[offset %in% "0", .(cor.coef = cor.coef), by = .(series, reference)]
    # cor.colnames <- paste0("as.dated.", names(data.refs.list), ".corr")
    # agg.df[, (cor.colnames) := as.list(cor.coefs$cor.coef)]
    cor.coefs <- this.series[offset %in% "0", .(reference, cor.coef)]
    cor.colnames <- paste0("as.dated.", names(data.refs.list), ".corr")
    agg.df[, (cor.colnames) := as.list(
      cor.coefs$cor.coef[match(names(data.refs.list), cor.coefs$reference)])]

    agg.df[, message := fifelse(nrow(these.refs[ncorrs == 0]) == 0,
                                paste("checked against",
                                      n.refs,
                                      fifelse(n.refs == 1,
                                              "reference", "references")),
                                fifelse(nrow(these.refs[ncorrs == 0]) == nrow(these.refs),
                                        paste("not checked against any",
                                              "references due to overlap <", min.overlap, "yrs"),
                                        paste("not checked against",
                                              paste(these.refs[ncorrs == 0]$reference,
                                                    sep = ifelse(nrow(these.refs) > 2, ", ", "& ")),
                                              "due to overlap <", min.overlap, "yrs")))]


    agg.df

  }) |> data.table::rbindlist()

  seg.max <- if (nrow(seg.res)) {
    seg.res[, .(n = max(as.integer(sub("/.*", "", agreement)))), by = series]
  } else {
    data.table::data.table(series = character(), n = integer())
  }
  cor.res.message[seg.max, on = "series", seg.n := i.n]
  cor.res.message[, seg.flag := paste(data.table::fcoalesce(seg.n, 0L), n.refs, sep = "/")]
  cor.res.message[, seg.n := NULL]

  # Now we can generate different subsets of cor.res.message to represent different categories
  # of caution / or confidence in the dating. c("Lowest conf.", "Medium conf.", "Highest conf.")
  all.good <- paste(n.refs, n.refs, sep = "/")

  hi <- cor.res.message[best.as.dated %in% all.good &
                          sig.as.dated %in% all.good &
                          seg.flag %in% paste("0", n.refs, sep = "/"),]

  lo <- cor.res.message[best.as.dated %in% all.bad |
                          sig.as.dated %in% all.bad,]

  mi <- cor.res.message[!series %in% c(hi$series, lo$series),]


  cor.summ.list <- list(lo, mi, hi, cor.res.message)
  names(cor.summ.list) <- c("Lowest conf.", "Medium conf.", "Highest conf.", "All")

  # Prep the data subsetting based on the results
  if (out.format %in% "rwl") {
    data.rwl <- longer_rwl(data = data,
                           series.name = "series",
                           dat.name = colnames(data)[grep("rw", colnames(data))])

    rw.out.list <- list(data.rwl[, lo$series, drop = FALSE],
                        data.rwl[, mi$series, drop = FALSE],
                        data.rwl[, hi$series, drop = FALSE])
  } else {
    rw.out.list <- list(data[series %in% lo$series,],
                        data[series %in% mi$series,],
                        data[series %in% hi$series,])
  }

  names(rw.out.list) <- c("Lowest", "Medium", "Highest")

  ## The final output is a named list
  out.list <- list(cor.summ.list, cor.res, rw.out.list, seg.res)
  names(out.list) <- c("Dating conf.", "Corr. result list", "Subset RW list",
                       "Segment flags")

  out.list

} ## End of function
