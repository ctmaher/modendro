#' Convert long format ring-width data to rwl-format
#'
#' @description
#' Simple function to convert long format tree ring data (3 columns: year, series, rw) to
#' rwl-format data (columns are series, rows are years). This is useful for grouping and
#' reassembling data. Note that \code{\link{longer_rwl}} does not enforce series name length limits
#' like a true rwl does.
#'
#'
#' @param df A long-format data.frame with "year" & "series" columns and a data column, specified
#' by the dat.name argument
#' @param series.name Character vector of length 1 for the column name representing series IDs.
#' E.g., "series", "tree", or "sample".
#' @param dat.name A character vector specifying the name of the tree ring data
#' (e.g., "rw", "bai.cm2", "d13C", etc.)
#'
#' @return A data.frame with series names as columns and years as rownames.
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
#' ca533.long <- rwl_longer(rwl = ca533, series.name = "series", dat.name = "rw.mm", trim = TRUE)
#' head(ca533.long)
#'
#' # Convert it back to rwl-like format
#'
#' ca533.rwl <- longer_rwl(df = ca533.long, series.name = "series", dat.name = "rw.mm")
#' head(ca533.rwl)
#' all(ca533.rwl %in% ca533)
#'


longer_rwl <- function(df = NULL,
                       series.name = NULL,
                       dat.name = NULL) {
  ## Error catching & warnings
  stopifnot(
    "df is not an object of class 'data.frame', or 'matrix'" =
      data.class(df) %in% "data.frame" |
      data.class(df) %in% "matrix"
  )

  stopifnot(
    "series.name must be a character vector of length 1" =
      is.character(series.name) |
      length(series.name) == 1
  )

  stopifnot(
    "dat.name must be a character vector of length 1" =
      is.character(dat.name) |
      length(dat.name) == 1
  )

  stopifnot(
    "df does not contain column 'year'" =
      any(colnames(df) %in% "year")
  )

  stopifnot(
    "df does not contain column matching series.name argument" =
      any(colnames(df) %in% series.name)
  )

  stopifnot(
    "df does not contain column matching dat.name argument" =
      any(colnames(df) %in% dat.name)
  )


  # Get the series names before we add a year column
  series <- unique(df[, series.name])


  wide.rwl <- stats::reshape(
    data = df[, c("year", series.name, dat.name)],
    idvar = "year",
    ids = .data[["year"]], # This was causing a "no visible binding for global variable" note
    times = series,
    timevar = series.name,
    #varying = as.list(series),
    direction = "wide"#,
    #new.row.names = min(df[, "year"]):max(df[, "year"])
  )

  rownames(wide.rwl) <- wide.rwl[, "year"]
  wide.rwl <- wide.rwl[order(wide.rwl[, "year"]),]
  wide.rwl <- wide.rwl[,!(colnames(wide.rwl) %in% "year")]
  colnames(wide.rwl) <- series
  class(wide.rwl) <- c("rwl","data.frame")
  wide.rwl
  }
