#' Convert long format ring-width data to rwl-format
#'
#' @description
#' Simple function to convert long format tree ring data (3 columns: year, series, rw) to rwl-format
#' data (columns are series, rows are years). This is useful for grouping and reassubling data.
#'
#'
#' @param df A long-format data.frame with 3 columns ("year", "series", "rw")
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


longer_rwl <- function(df = NULL) {
  ## Error catching & warnings
  #
  stopifnot(
    "df is not an object of class 'data.frame', or 'matrix'" =
      data.class(df) %in% "rwl" |
      data.class(df) %in% "data.frame" |
      data.class(df) %in% "matrix"
  )

  # Get the series names before we add a year column
  series <- unique(df[, "series"])


  wide.rwl <- stats::reshape(
    df,
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
