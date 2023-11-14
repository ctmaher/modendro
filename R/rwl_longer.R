#' Convert rwl-format dataframes/matrices to long format
#'
#' @description
#' Simple function to convert rwl-format data (columns are series, rows are yeqrs) to long format
#' (3 columns: year, series, rw). This is useful for getting your data ready for plotting in ggplot or
#' for correlation analyses.
#'
#' @param rwl A rwl-type data.frame (e.g., read in by \code{\link[dplR]{read.rwl}}). Essentially a data.frame with columns names as series IDs and years as rownames.
#' @param dat.name Character vector of length 1 representing the column name you want for the the tree ring data. Default is "rw".
#' @param omit.NAs Logical vector indicating whether to keep NAs in the output or not. Default is TRUE. This can cause issues if missing rings are represented with NA.
#'
#' @return A data.frame with 3 columns
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
#' ca533.long <- rwl_longer(ca533, dat.name = "rw.mm", omit.NAs = TRUE)
#' head(ca533.long)


rwl_longer <- function(rwl = NULL,
                       dat.name = "rw",
                       omit.NAs = TRUE) {
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

  stopifnot("omit.NAs must be a logical vector (TRUE or FALSE)" =
              is.logical(omit.NAs))

  # Get the series names before we add a year column
  series.cols <- colnames(rwl)

  # Make a year column based on the rownames
  rwl$year <- rownames(rwl) |> as.numeric()


  long.rwl <- stats::reshape(
    rwl,
    idvar = "year",
    ids = x$year,
    times = series.cols,
    timevar = "series",
    v.names = "rw",
    varying = list(series.cols),
    direction = "long",
    new.row.names = 1:(length(series.cols) * nrow(rwl))
  )
  if (omit.NAs == TRUE) {
    long.rwl |>
      na.omit()
  } else {
    long.rwl
  }
}
