#' Convert long format ring-width data to rwl-format
#'
#' @description
#' Simple function to convert long format tree-ring data to wl-format data (columns are series, rows
#'  are years). This is useful for export (e.g., with \code{\link[dplR]{write.rwl}}) or for grouping
#' and reassembling data. Note that \code{\link{longer_rwl}} does not enforce series name length
#' limits like a true rwl does.
#'
#'
#' @param data A long-format data.table or data.frame with a "year" column, an ID column
#' (specified by series.name), and a data column (specified by the dat.name argument)
#' @param series.name Character vector of length 1 for the column name representing series IDs.
#' E.g., "series", "tree", or "sample".
#' @param dat.name A character vector specifying the name of the tree-ring data
#' (e.g., "rw", "bai.cm2", "d13C", etc.)
#'
#' @return An rwl-data.frame with series names as columns and years as rownames.
#'
#' @import data.table
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
#' # Convert it back to rwl format
#'
#' ca533.rwl <- longer_rwl(data = ca533.long, series.name = "series", dat.name = "rw.mm")
#' head(ca533.rwl)
#' all(ca533.rwl %in% ca533)
#'
#' data.class(ca533.rwl)
#'


longer_rwl <- function(data = NULL,
                       series.name = NULL,
                       dat.name = NULL) {
  . <- j_variable <- IDs <- series <- year <- dup.years <- cont.check <- NULL #

  #### Error catching ####
  stopifnot(
    "data cannot already be an 'rwl'" =
      !(data.class(data) %in% "rwl")
  )

  stopifnot(
    "data must be an object of class 'data.table' or 'data.frame'" =
      is.data.frame(data)
  )

  stopifnot(
    "series.name must be a character vector of length 1" =
      is.character(series.name) &&
      length(series.name) == 1
  )

  stopifnot(
    "dat.name must be a character vector of length 1" =
      is.character(dat.name) &&
      length(dat.name) == 1
  )

  stopifnot(
    "data does not contain column 'year'" =
      any(colnames(data) %in% "year")
  )

  stopifnot(
    "data does not contain column matching series.name argument" =
      any(colnames(data) %in% series.name)
  )

  stopifnot(
    "data does not contain column matching dat.name argument" =
      any(colnames(data) %in% dat.name)
  )

  # Prevent user from specifying series name & dat.name the same and neither can be "year"
  stopifnot(
    "series.name and dat.name must differ and neither can be 'year' or 'Year'" =
      !anyDuplicated(c(series.name, dat.name, "year", "Year"))
  )

  # Check for NAs in year and series.name columns
  stopifnot(
    "'year' variable cannot contain NAs" =
      !anyNA(data$year)
  )

  stopifnot(
    "series.name variable cannot contain NAs" =
      !anyNA(data[[series.name]])
  )

  # Make sure year column contains integer values
  stopifnot(
    "'year' must be whole numbers" =
      is.numeric(data$year) &&
      all(data$year %% 1 == 0)
  )

  # If not a data.table already, convert it. If it is, copy it
  data <- if (data.table::is.data.table(data)) {
    data.table::copy(data)
  } else {
    data.table::as.data.table(data)
  }

  data[, year := as.integer(year)]

  # Ensure proper sorting of years within series
  data.table::setorderv(data, c(series.name, "year"))

  # Make sure dat.name is numeric
  stopifnot(
    "dat.name variable must be numeric" =
      is.numeric(data[[dat.name]])
  )

  # Make sure that years are not duplicated and continuous within each series
  dup.year.test <- data[, .(dup.years = anyDuplicated(year)), by = series,
                        env = list(series = series.name)]
  if (any(dup.year.test$dup.years > 0)) {
    stop("These series have duplicated years: ",
         paste(dup.year.test[dup.years > 0, ][[series.name]], collapse = ", "), call. = FALSE)
  }

  cont.year.test <- data[, .(cont.check = any(diff(year) != 1)), by = series,
                         env = list(series = series.name)]
  if (any(cont.year.test$cont.check)) {
    stop("These series have non-continuous years: ",
         paste(cont.year.test[(cont.check),][[series.name]], collapse = ", "), call. = FALSE)
  }


  #### Process the data ####
  # Get the series names
  ser.names <- unique(data[, IDs, env = list(IDs = series.name)]) |> as.character()

  # long-to-wide with dcast()
  wide.rwl <- data.table::dcast(
    data = data,
    formula = year ~ series,
    value.var = dat.name,
    env = list(series = series.name)
  )
  # Get full sequence of years to ensure no discontinuities in year/rownames in rwl output
  yrs <- seq(min(wide.rwl$year), max(wide.rwl$year))
  wide.rwl <- wide.rwl[data.table::data.table(year = yrs), on = "year"]
  # Convert to data.frame
  wide.rwl <- as.data.frame(wide.rwl)
  # Transfer years to rownames
  rownames(wide.rwl) <- wide.rwl[, "year"]
  # Arrange by series name order and remove year column
  wide.rwl <- wide.rwl[, ser.names, drop = FALSE]
  # Assign rwl class
  class(wide.rwl) <- c("rwl","data.frame")
  # Return
  wide.rwl

} ## End of function
