#' Convert long format ring-width data to rwl-format
#'
#' @description
#' Simple function to convert long format tree ring data (3 columns: year, series, rw) to rwl-format
#' data (columns are series, rows are years). This is useful for grouping and reassubling data.
#'
#'
#' @param df A long-format data.frame with "year" & "series" columns and a data column, specified
#' by the dat.name argument
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
#' ca533.long <- rwl_longer(rwl = ca533, dat.name = "rw.mm", trim = TRUE)
#' head(ca533.long)
#'
#' # Convert it back to rwl-like format
#'
#' ca533.rwl <- longer_rwl(df = ca533.long, dat.name = "rw.mm")
#' head(ca533.rwl)
#' all(ca533.rwl %in% ca533)
#'


longer_rwl <- function(df = NULL, dat.name = NULL) {
  ## Error catching & warnings
  stopifnot(
    "df is not an object of class 'data.frame', or 'matrix'" =
      data.class(df) %in% "rwl" |
      data.class(df) %in% "data.frame" |
      data.class(df) %in% "matrix"
  )

  stopifnot(
    "dat.name must be a character vector of length 1" =
      is.character(dat.name) |
      length(dat.name) == 1
  )


  stopifnot(
    "df does not contain columns 'year' or 'series'" =
      any(colnames(df) %in% "year") |
      any(colnames(df) %in% "series")
  )

  stopifnot(
    "df does not contain column matching dat.name argument" =
      any(colnames(df) %in% dat.name)
  )


  # Get the series names before we add a year column
  series <- unique(df[, "series"])


  wide.rwl <- stats::reshape(
    df[, c("year", "series", dat.name)],
    idvar = "year",
    ids = x$year,
    times = series,
    timevar = "series",
    varying = list(series),
    direction = "wide"#,
    #new.row.names = min(df[, "year"]):max(df[, "year"])
)

  rownames(wide.rwl) <- wide.rwl[, "year"]
  wide.rwl <- wide.rwl[order(wide.rwl[, "year"]),]
  wide.rwl[,!(colnames(wide.rwl) %in% "year")]

  }
