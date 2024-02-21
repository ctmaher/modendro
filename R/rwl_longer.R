#' Convert rwl-format dataframes/matrices to long format
#'
#' @description
#' Simple function to convert rwl-format data (columns are series, rows are yeqrs) to long format
#' (3 columns: year, series, rw). This is useful for getting your data ready for plotting in ggplot or
#' for correlation analyses.
#'
#'
#' @param rwl A rwl-type data.frame (e.g., read in by \code{\link[dplR]{read.rwl}}). Essentially a data.frame with columns names as series IDs and years as rownames.
#' @param series.name Character vector of length 1 for the column name you want for the series IDs. Default is "series", but perhaps "tree" or "sample" is more appropriate for your work.
#' @param dat.name Character vector of length 1 for the column name you want for the the tree ring data. Default is "rw".
#' @param trim Logical vector indicating whether to trim off NA sequences at the beginning or end of individual series. Default is TRUE. Will not remove missing rings that are represented with NA (i.e., NA values within the series).
#'
#' @return A data.frame with 3 columns: 1) "year", 2) "series", 3) dat.name ("rw", "rw.mm", "bai.mm" or whatever you name this)
#'
#' @import stats
#' @import dplR
#' @export
#'
#' @examples
#' library(dplR)
#' # Bristlecone pine tree ring collection from Campito mountain, White Mountains, CA
#' # Many 0 value rings in this collection
#' data("ca533")
#'
#' ca533.long <- rwl_longer(rwl = ca533, dat.name = "rw.mm", trim = TRUE)
#' head(ca533.long)


rwl_longer <- function(rwl = NULL,
                       series.name = "series",
                       dat.name = "rw",
                       trim = TRUE) {
  ## Error catching & warnings
  #
  stopifnot(
    "rwl is not an object of class 'rwl', 'data.frame', or 'matrix'" =
      data.class(rwl) %in% "rwl" |
      data.class(rwl) %in% "data.frame" |
      data.class(rwl) %in% "matrix"
  )

  #
  stopifnot(
    "rwl has no rownames (must be years only) or no colnames (must be series IDs only)" =
      !is.null(rownames(rwl)) |
      !is.null(colnames(rwl))
  )

  #
  stopifnot(
    "dat.name can only be a character vector of length 1" =
      length(dat.name) == 1 |
      is.character(dat.name)
  )

  stopifnot("trim must be a logical vector (TRUE or FALSE)" =
              is.logical(trim))

  # Get the series names before we add a year column
  series.cols <- colnames(rwl)

  # Make a year column based on the rownames
  rwl[, "year"] <- rownames(rwl) |> as.numeric()


  long.rwl <- stats::reshape(
    rwl,
    idvar = "year",
    ids = x$year,
    times = series.cols,
    timevar = series.name,
    v.names = dat.name,
    varying = list(series.cols),
    direction = "long",
    new.row.names = 1:(length(series.cols) * nrow(rwl))
  )
  if (trim == TRUE) {
    # Find the earliest & latest years that have rw values
    long.rwl.nonas <-  na.omit(long.rwl)
    long.rwl.nonas.split <- split(long.rwl.nonas, f = long.rwl.nonas[, series.name])
    year.spans <- lapply(long.rwl.nonas.split, FUN = \(series) {
      data.frame(series = unique(series[, series.name]),
                 min.year = series[, "year"] |> min(),
                 max.year = series[, "year"] |> max())
    })
    long.rwl.split <- split(long.rwl, f = long.rwl[, series.name])
    long.rwl.trim <- mapply(FUN = \(series.dfs, end.years) {
      series.dfs[series.dfs[, "year"] %in% end.years[, "min.year"]:end.years[, "max.year"], ]
    }, series.dfs = long.rwl.split, end.years = year.spans,
    SIMPLIFY = FALSE) |> do.call(what = "rbind")
    long.rwl.trim

    if (any(is.na(long.rwl.trim[,dat.name])) == TRUE) {
      warning("NA values found within tree ring series. Are NAs used to represent missing rings? 0 would be better.")
    }

    } else { # just return everything, NAs & all.
    long.rwl
  }
}
