#' Statistical cross-dating checking using methods inspired by CooRecorder/CDendro and COFECHA.
#'
#'
#' @description
#' This function implements a custom statistical check on the dating of a collection of tree-ring
#' series. It is inspired by operations in CDendro (Larsson & Larsson 2024) and COFECHA,
#' but does not do all of what either of those programs can do. \code{\link{xd_check}} does simple
#' whole-series cross-correlations between leave-one-out (LOO) chronologies and 0-n reference
#' chronologies. You provide collections of raw tree-ring widths (not pre-made chronologies), and
#' xd_check will standardize/normalize all series and then generate the chronologies internally.
#' If you only provide the `data` collection, LOO comparisons are all that you'll get. Adding
#' appropriate (e.g., properly cross-dated, same species, same region) external reference
#' collections adds a lot to the confidence of your checks. The main output of the function is a
#' summary based on how many references have best and significant correlations with the series as
#' dated and are ranked lowest, medium, and highest confidence.
#'
#' One potential use for \code{\link{xd_check}} is to generate "clean" references for subsequent
#' runs of the function. You can do this by supplying the reference collection to the `data`
#' argument and not supplying any of the reference arguments (i.e., LOO only mode). Set the other
#' arguments as you please and then take the "Highest" subset as one of the `refx` arguments in the
#' next run. Note that this assumes that you have plenty of good series in your collection!
#'
#' A note of awareness: the chronology building process in \code{\link{xd_check}} is different than
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
#' "data.frame" (long format). This is required. If "data.frame" you need 3 columns: "series",
#' "year", & "rw".
#' @param ref A single reference collection or a list of reference collections of raw tree-ring
#' series (NOT a chronology!). Can be class "rwl" or "data.frame" (long format). Optional, but more
#' references are better here. If "data.frame" you need 3 columns: "series", "year", & "rw". If a
#' list, can be a named list so that you can give unique names to your references.
#' @param std.method The method of standardization/normalization to apply to the series and
#' references. Default (and currently only) option is "p2yrsL", same default as CooRecorder/CDendro.
#' @param max.offset The maximum offset (in years) you want to check in either direction from the
#' years as dated. Minimum is 5 years. No true maximum limit, but overlap between series and
#' references must be at least 20 years. Default is 10 years.
#' @param p.thresh The p-value threshold to determine if a correlation is significant. Default is
#' the common 0.05.
#' @param out.format The format of subsets of data based on the dating confidence ranking. Options
#' are "rwl" or "long".
#'
#' @details
#' \code{\link{xd_check}} does a standardization/normalization on all individual series within the
#' `data` and `refn` collections, then it computes chronologies using Tukey's Biweight Robust mean
#' (via \code{\link[DescTools]{TukeyBiweight}}). The LOO method is common in dendrochronology as
#' the method to compute interseries correlations for a tree-ring collection. In LOO, each series
#' has its own chronology made up of the remaining series.
#'
#' The basic ranking of dating confidence here is based on two basic criteria: 1) is the
#' whole-series correlation best as dated? and 2) is the correlation significant as dated? These
#' criteria are assessed for each reference comparison. The more references you have, the higher
#' your overall confidence should be.
#'
#' @return A list containing 1) the partitioned results based on relative confidence (itself a
#' 4-element list), 2) the raw correlation results in a nested list for easy viewing in RStudio,
#' and 3) subsets of the data collection corresponding to the relative confidence (format specified
#' with `out.format`).
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
                     out.format = "rwl"
) {
  ## Error catching
  # The function cannot do anything if data is NULL
  stopifnot(
    "You must provide tree-ring data to the data argument" =
      !is.null(data)
  )

  # Data must also be in the proper format
  stopifnot(
    "data is not an object of class 'data.frame', or 'matrix'" =
      is.null(data) |
      data.class(data) %in% "rwl" |
      data.class(data) %in% "data.frame" |
      data.class(data) %in% "matrix"
  )

  # Convert to long format if read in as an rwl
  if (any(class(data) %in% "rwl")) {
    data <- rwl_longer(rwl = data, series.name = "series", dat.name = "rw", trim = TRUE,
                       new.val.internal.na = 0)
  }

  # If the data and references are long-format data.frames, we need certain columns
  ##
  stopifnot(
    "data does not have a year column?" =
      any(colnames(data) %in% c("Year", "year")) == TRUE
  )

  # make sure that "year" columns are labelled as such
  colnames(data)[which((colnames(data) %in% c("Y", "y")) == T)] <- "year"

  # general check
  stopifnot(
    "data needs to contain these 3 columns: 'series', 'year', & 'rw'" =
      any(colnames(data) %in% "series") |
      any(colnames(data) %in% "year") |
      length(grep("rw", colnames(data))) > 0
  )

  # make sure that "rw" column is labelled as such
  colnames(data)[grep("rw", colnames(data))] <- "rw"

  # Make sure these columns are the right kind of data (year & rw need to be numeric)
  stopifnot(
    "year variable in data not numeric or integer" =
      is.numeric(data[, "year"]) |
      is.integer(data[, "year"]) |
      is.double(data[, "year"])
  )

  stopifnot(
    "ring width ('rw') variable in data not numeric or integer" =
      is.numeric(data[, "rw"]) |
      is.integer(data[, "rw"]) |
      is.double(data[, "rw"])
  )

  ## Detailed checking of the data and references to determine if they are rwl format or
  # long format. If rwls are read in using dplR::read.rwl they will have class "rwl" already.
  # What else can signify this? Not sure of a robust test
  # Start with the simple test. For this function, I will convert rwls to long format.

  # The references should only go though these checks if they exist.
  # For the list case, it might make the most since to check through everything, then check through
  # the results once to give the entire result at once (maybe).
  if (!is.null(ref)) {
    if (is.list(ref) &&
        !(data.class(ref) %in% "rwl" ||
          data.class(ref) %in% "data.frame" ||
          data.class(ref) %in% "matrix")) { # If ref is a list
      ref.names <- names(ref) # keep the names for later
      ref <- lapply(seq_along(ref), FUN = \(i) {
        this.ref <- ref[[i]]

        if (!(data.class(this.ref) %in% "rwl" ||
              data.class(this.ref) %in% "data.frame" ||
              data.class(this.ref) %in% "matrix")) {

          stop(sprintf("ref element %d is not an object of class 'data.frame', or 'matrix'", i))
        }

        # Convert to long format if a rwl
        if (any(class(this.ref) %in% "rwl")) {
          this.ref <- rwl_longer(rwl = this.ref,
                                 series.name = "series",
                                 dat.name = "rw",
                                 trim = TRUE,
                                 new.val.internal.na = 0)
        }

        ## Columns check
        if (any(colnames(this.ref) %in% c("Year", "year")) == FALSE) {

          stop(sprintf("ref element %d does not have a year column",
                       i))
        }

        # make sure that "year" columns are labelled as such
        colnames(this.ref)[which((colnames(this.ref) %in%
                                    c("Year", "year")) == T)] <- "year"

        # general check
        if (!(any(colnames(this.ref) %in% "series") ||
              any(colnames(this.ref) %in% "year") ||
              length(grep("rw", colnames(this.ref))) > 0)) {

          stop(sprintf("ref element %d needs to contain these 3 columns: 'series', 'year', & 'rw'",
                       i))
        }

        # make sure that "rw" column is labelled as such
        colnames(this.ref)[grep("rw", colnames(this.ref))] <- "rw"

        # Make sure these columns are the right kind of data (year & rw need to be numeric)
        if (!(is.numeric(this.ref[, "year"]) ||
              is.integer(this.ref[, "year"]) ||
              is.double(this.ref[, "year"]))) {

          stop(sprintf("year variable in ref element %d not numeric or integer",
                       i))
        }

        if (!(is.numeric(this.ref[, "rw"]) ||
              is.double(this.ref[, "rw"]))) {

          stop(sprintf("ring width ('rw') in ref element %d not numeric",
                       i))
        }
        this.ref
      })

      names(ref) <- ref.names

    } else { # if ref is not a list

      stopifnot(
        "ref is not an object of class 'data.frame', 'matrix', or 'list'" =
          data.class(ref) %in% "rwl" |
          data.class(ref) %in% "data.frame" |
          data.class(ref) %in% "matrix"
      )

      # Convert to long format if a rwl
      if (any(class(ref) %in% "rwl")) {
        ref <- rwl_longer(rwl = ref, series.name = "series", dat.name = "rw", trim = TRUE,
                          new.val.internal.na = 0)
      }

      ## Columns check
      stopifnot(
        "ref does not have a year column?" =
          any(colnames(ref) %in% c("Year", "year")) == TRUE
      )

      # make sure that "year" columns are labelled as such
      colnames(ref)[which((colnames(ref) %in% c("Year", "year")) == T)] <- "year"

      # general check
      stopifnot(
        "ref needs to contain these 3 columns: 'series', 'year', & 'rw'" =
          any(colnames(ref) %in% "series") |
          any(colnames(ref) %in% "year") |
          length(grep("rw", colnames(ref))) > 0
      )

      # make sure that "rw" column is labelled as such
      colnames(ref)[grep("rw", colnames(ref))] <- "rw"

      # Make sure these columns are the right kind of data (year & rw need to be numeric)
      stopifnot(
        "year variable in ref not numeric or integer" =
          is.numeric(ref[, "year"]) |
          is.integer(ref[, "year"]) |
          is.double(ref[, "year"])
      )

      stopifnot(
        "ring width ('rw') variable in ref not numeric or integer" =
          is.numeric(ref[, "rw"]) |
          is.integer(ref[, "rw"]) |
          is.double(ref[, "rw"])
      )
    }
  }

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
      is.character(out.format) &
      length(out.format) == 1 &
      out.format %in% c("rwl", "long")
  )


  ### Start the process
  # Steps: 1) convert to p2yrsL, AR, or ARIMA. 2) Build chronologies. 3) Run correlations.
  # 4) Assemble the results

  # Make a list of the data and the ref(s)
  if (is.list(ref) &&
      !(data.class(ref) %in% "rwl" ||
        data.class(ref) %in% "data.frame" ||
        data.class(ref) %in% "matrix")) {

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

  } else { # Ref is not a list

    dat.ref.list <- list(data, ref)
    names(dat.ref.list) <- c("data", "ref1")
  }

  # Remove the null values
  dat.ref.list <- dat.ref.list[!unlist(lapply(dat.ref.list, FUN = \(x)is.null(x)))]


  # Apply the standardization to each element of the list
  std.all <- lapply(dat.ref.list, FUN = \(x) {
    if (std.method %in% "p2yrsL") {
      std.df <- longer_rwl(df = x, series.name = "series", dat.name = "rw") |>
        p2yrsL() |>
        rwl_longer()
    }

    if (std.method %in% "AR") {
      # std.df <- longer_rwl(df = x, series.name = "series", dat.name = "rw") |>
      #   p2yrsL() |>
      #   rwl_longer()
    }

    if (std.method %in% "ARIMA") {
      # std.df <- longer_rwl(df = x, series.name = "series", dat.name = "rw") |>
      #   p2yrsL() |>
      #   rwl_longer()
    }

    colnames(std.df)[which(colnames(std.df) %in% "rw")] <- "std.series"
    std.df
  })


  ## Set up the sliding windows (i.e., lags) based on the provided max.offset
  lags <- seq(-max.offset, max.offset, by = 1)
  # These then have to be applied to each series + chron combination

  # set min.overlap - the way I'm handling this is causing problems I think
  min.overlap <- 20

  # Do the leave-one-out chronology method on the data (n chrons = n series)
  data.series.split <- split(std.all[["data"]], f = std.all[["data"]]$series)

  data.loo <- lapply(data.series.split,
                     FUN = \(this.series) {

                       #series$std.series <- ifelse(is.na(series$std.series), 0, series$std.series)

                       this.seriesID <- unique(this.series[,"series"])

                       loo.chron1 <- aggregate(std.series ~ year,
                                               data = std.all[["data"]][!(std.all[["data"]]$series
                                                                          %in%
                                                                            this.seriesID),],
                                               FUN = \(x) DescTools::TukeyBiweight(x, na.rm = TRUE))
                       colnames(loo.chron1)[colnames(loo.chron1) %in% "std.series"] <- "std.chron"

                       loo.depth <- aggregate(std.series ~ year,
                                              data = std.all[["data"]][!(std.all[["data"]]$series
                                                                         %in%
                                                                           this.seriesID),],
                                              length)
                       colnames(loo.depth)[colnames(loo.depth) %in% "std.series"] <- "samp.depth"

                       loo.chron <- merge(loo.chron1, loo.depth, by = "year")

                       # Get the possible lag years for the series
                       lags.df <- sapply(lags, FUN = \(n.lag) {
                         this.lag <- data.frame(this.series[, "year"] + n.lag)
                         colnames(this.lag) <- paste0("lag", n.lag, ".year")
                         this.lag
                       }) |> do.call(what = "cbind") |> as.data.frame()

                       # Merge the series with the respective loo chronology with each lag
                       merged.loo.lags <- apply(lags.df, MARGIN = 2, FUN = \(this.lag) {
                         merge(cbind(this.series[, c("series", "std.series"), drop = FALSE],
                                     this.lag),
                               loo.chron,
                               by.x = "this.lag", by.y = "year")
                       })


                       # Set min overlap
                       lag.dfs <- lapply(merged.loo.lags, FUN = \(this.lag) {
                         if (nrow(this.lag) < min.overlap) {
                           NA
                         } else {
                           this.lag
                         }
                       })

                       lag.dfs <- lag.dfs[!is.na(lag.dfs)]

                       # get lag identifiers
                       lag.dfs <- mapply(FUN = \(l, lnames) {
                         l[,"lag"] <- gsub(".*lag([^.]+)\\.year.*", "\\1", lnames)
                         l
                       }, l = lag.dfs, lnames = names(lag.dfs), SIMPLIFY = FALSE)

                       lag.dfs

                     })


  # Make single chronologies for the rest - if they exist
  if (length(std.all) > 1) {
    # This will be a list of length 1 to 3 (or n refs)
    std.chrons <- lapply(std.all[!(names(std.all) %in% "data")],
                         FUN = \(this.ref) {
                           std.chron <- aggregate(std.series ~ year, data = this.ref,
                                                  FUN = \(x) {
                                                    DescTools::TukeyBiweight(x, na.rm = TRUE)
                                                  })
                           colnames(std.chron)[which(colnames(std.chron) %in% "std.series")] <-
                             "std.chron"
                           std.chron
                         })

    data.refs <- lapply(std.chrons, FUN = \(this.chron) {

      this.ref <- lapply(data.series.split,
                         FUN = \(this.series) {

                           # Get the possible lag years
                           lags.df <- sapply(lags, FUN = \(n.lag) {
                             this.lag <- data.frame(this.series[, "year"] + n.lag)
                             colnames(this.lag) <- paste0("lag", n.lag, ".year")
                             this.lag
                           }) |> do.call(what = "cbind") |> as.data.frame()

                           # Merge the series with the respective loo chronology with each lag
                           merged.lags <- apply(lags.df, MARGIN = 2, FUN = \(this.lag) {
                             merge(cbind(this.series[, c("series", "std.series"), drop = FALSE],
                                         this.lag),
                                   this.chron,
                                   by.x = "this.lag", by.y = "year")
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
                             l[,"lag"] <- gsub(".*lag([^.]+)\\.year.*", "\\1", lnames)
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


  ## Run the correlations between each series and each reference
  # The output here is a single data.frame containing all the results - can make any other type of
  # organization from this point
  cor.res.raw <- mapply(FUN = \(this.combo, these.names) {
    lapply(this.combo, FUN = \(this.series) {
      all.lags <- lapply(this.series, FUN = \(this.lag) { # Do correlations on this layer
        overlap <- nrow(this.lag)
        if (overlap < min.overlap) {
          cor.coef <- NA
          p.val <- NA
        } else {
          suppressWarnings(
            this.cor.res <- cor.test(this.lag[, "std.series"],
                                     this.lag[, "std.chron"],
                                     method = "spearman")
          )
          # Round the results for easy viewing
          cor.coef <- round(this.cor.res$estimate[[1]], digits = 2)
          p.val <- round(this.cor.res$p.value[[1]], digits = 3)
        }

        res.df <- data.frame(series = unique(this.lag[, "series"]),
                             reference = these.names,
                             offset = unique(this.lag[, "lag"]),
                             overlap = overlap,
                             cor.coef = cor.coef,
                             p.val = p.val)
        res.df$T.val <- round(
          abs(res.df$cor.coef) *
            sqrt((res.df$overlap-2)/1-(res.df$cor.coef)^2),
          digits = 2)

        res.df
      }) |> do.call(what = "rbind")

      # Calculate correct sugg.ID and sugg.OD that can go beyond the references
      # This will get the LOO and the reference values - only the LOO is guaran
      # all.lags$sugg.ID <- all.lags[all.lags$offset %in% "0", "sugg.ID"] +
      #   as.numeric(all.lags$offset)
      #
      # all.lags$sugg.OD <- all.lags[all.lags$offset %in% "0", "sugg.OD"] +
      #   as.numeric(all.lags$offset)
      #
      all.lags

    }) |> do.call(what = "rbind")

  }, this.combo = data.refs.list, these.names = names(data.refs.list),
  SIMPLIFY = FALSE) |> do.call(what = "rbind")


  # Calculate the sugg.ID and sugg.OD using the original series data
  data.ID <- aggregate(year ~ series, data = dat.ref.list$data, FUN = \(x) min(x))
  data.OD <- aggregate(year ~ series, data = dat.ref.list$data, FUN = \(x) max(x))
  colnames(data.ID)[colnames(data.ID) %in% "year"] <- "sugg.ID"
  colnames(data.OD)[colnames(data.OD) %in% "year"] <- "sugg.OD"

  data.ID.OD <- merge(data.ID, data.OD, by = "series")
  cor.res.raw1 <- merge(cor.res.raw, data.ID.OD, by = "series", no.dups = FALSE)

  cor.res.raw1$sugg.ID <- cor.res.raw1$sugg.ID + as.numeric(cor.res.raw1$offset)
  cor.res.raw1$sugg.OD <- cor.res.raw1$sugg.OD + as.numeric(cor.res.raw1$offset)
  # Reorganize the columns
  cor.res.raw1 <- cor.res.raw1[, c("series","reference","offset","sugg.ID","sugg.OD",
                                   "overlap","cor.coef","p.val","T.val")]

  ## Sort/filter/arrange this massive amount of data into a useful output
  # First step is to reorganize into a list that can be perused easily by the user
  cor.res <- lapply(split(cor.res.raw1, f = cor.res.raw1$series), FUN = \(this.series) {
    lapply(split(this.series, f = this.series[, "reference"]), FUN = \(this.combo) {
      this.combo <- this.combo[order(this.combo$cor.coef, decreasing = TRUE),]
      this.combo#[1:5,]
    })
  })

  # The next step is applying some decision-making to the results so that we can turn a manual
  # decision process into an automated one. Ideal output is to distill results into 1 row per
  # series.

  # Perhaps the best way (consider if all the references have best as dated and are significant)
  n.refs <- length(data.refs.list)
  all.bad <- paste("0", n.refs, sep = "/")

  cor.res.message <- lapply(split(cor.res.raw1, f = cor.res.raw1$series), FUN = \(this.series) {

    this.series.sum <- lapply(split(this.series, f = this.series[, "reference"]),
                              FUN = \(this.combo) {
                                this.combo <- this.combo[order(this.combo$cor.coef,
                                                               decreasing = TRUE),]
                                rownames(this.combo) <- 1:nrow(this.combo)
                                as.dated <- this.combo[this.combo$offset %in% "0",]

                                data.frame(series = as.dated$series,
                                           best.as.dated = ifelse(which(this.combo$offset %in% "0")
                                                                  %in% "1", 1, 0),
                                           sig.as.dated = ifelse(as.dated$p.val < p.thresh, 1, 0))

                              }) |> do.call(what = "rbind")

    # Make sure it says 0/n if the overlap is inadequate to have anything

    # I need to figure out which of the references (including LOO) have correlations for
    # this.series as dated.
    these.refs <- rownames(this.series.sum)[which(names(data.refs.list) %in%
                                                    rownames(this.series.sum))]


    if (nrow(this.series.sum) < 1) { # No refs had adequate overlap
      agg.df <- data.frame(series = unique(this.series[, "series"]),
                           best.as.dated = all.bad,
                           sig.as.dated = all.bad,
                           as.dated.ID = unique(this.series[this.series$offset %in% "0",
                                                            "sugg.ID"]),
                           as.dated.OD = unique(this.series[this.series$offset %in% "0",
                                                            "sugg.OD"]))
      agg.df[, paste0("as.dated.", names(data.refs.list), ".cor")] <- NA

      agg.df$message <-  paste("not checked against any",
                               "references due to overlap < 20yrs")

    } else {
      agg.df <- aggregate(cbind(best.as.dated, sig.as.dated) ~ series, data = this.series.sum,
                          FUN = \(x) paste(sum(x), n.refs, sep = "/"), drop = FALSE)

      agg.df$as.dated.ID <- unique(this.series[this.series$offset %in% "0",
                                               "sugg.ID"])

      agg.df$as.dated.OD <- unique(this.series[this.series$offset %in% "0",
                                               "sugg.OD"])

      ## Add the correlation coefficients for as dated
      # Which refs do we have correlation coefs for?
      cor.coefs <- lapply(1:n.refs, FUN = \(i) {
        this.coef <- this.series[this.series$offset %in% "0" &
                                   this.series$reference %in% names(data.refs.list)[i], "cor.coef"]
        ifelse(length(this.coef) == 0, NA, this.coef)
      }) |> do.call(what = "c")

      agg.df[, paste0("as.dated.", names(data.refs.list), ".corr")] <- cor.coefs

      agg.df$message <- ifelse(nrow(this.series.sum) < n.refs,
                               paste("not checked against",
                                     n.refs - nrow(this.series.sum),
                                     "reference due to overlap < 20yrs"),
                               paste("checked against",
                                     n.refs,
                                     ifelse(n.refs == 1,
                                            "reference", "references")))
    }

    agg.df

  }) |> do.call(what = "rbind")

  # Now we can generate different subsets of cor.res.message to represent different categories
  # of caution / or confidence in the dating. c("Lowest conf.", "Medium conf.", "Highest conf.")
  all.good <- paste(n.refs, n.refs, sep = "/")
  lo <- cor.res.message[cor.res.message$best.as.dated %in% all.bad |
                          cor.res.message$sig.as.dated %in% all.bad,]

  mi <- cor.res.message[!((cor.res.message$best.as.dated %in% all.bad |
                             cor.res.message$sig.as.dated %in% all.bad) |
                            (cor.res.message$best.as.dated %in% all.good &
                               cor.res.message$sig.as.dated %in% all.good)),]

  hi <- cor.res.message[cor.res.message$best.as.dated %in% all.good &
                          cor.res.message$sig.as.dated %in% all.good,]

  cor.summ.list <- list(lo, mi, hi, cor.res.message)
  names(cor.summ.list) <- c("Lowest conf.", "Medium conf.", "Highest conf.", "All")

  # Prep the data subsetting based on the results
  if (out.format %in% "rwl") {
    data.rwl <- longer_rwl(data,
                           series.name = "series",
                           dat.name = colnames(data)[grep("rw", colnames(data))])

    rw.out.list <- list(data.rwl[, lo$series, drop = FALSE],
                        data.rwl[, mi$series, drop = FALSE],
                        data.rwl[, hi$series, drop = FALSE])
  } else {
    rw.out.list <- list(data[data$series %in% lo$series,],
                        data[data$series %in% mi$series,],
                        data[data$series %in% hi$series,])
  }

  names(rw.out.list) <- c("Lowest","Medium","Highest")

  ## The final output is a named list
  out.list <- list(cor.summ.list, cor.res, rw.out.list)
  names(out.list) <- c("Dating conf.", "Corr. result list", "Subset RW list")

  out.list

} # End of function
