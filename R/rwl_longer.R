#' Convert rwl-format data.frames/matrices to long format
#'
#' @description
#' Simple function to convert rwl-format data (columns are series, rows are years) to long format
#' (3 columns: year, series, rw). This is useful for getting your data ready for plotting in ggplot
#' or for other analyses.
#'
#' \code{\link{rwl_longer}} doesn't tolerate replicate column names. Column names are series names,
#' after all. If you have replicate column names in your rwl, you will get an error message
#' that lists the replicated columns. Fix these and run \code{\link{rwl_longer}} again.
#'
#' @param rwl A rwl-type data.frame (e.g., read in by \code{\link[dplR]{read.rwl}}). Essentially a
#' data.frame with columns names as series IDs and years as rownames.
#' @param series.name Character vector of length 1 for the column name you want for the series IDs.
#' Default is "series", but perhaps "tree" or "sample" is more appropriate for your work.
#' @param dat.name Character vector of length 1 for the column name you want for the the tree ring
#' data. Default is "rw".
#' @param new.val.internal.na Numeric vector giving a new value to replace NAs in the
#' middle of measurement series (i.e., some measurement protocols leave NAs for missing/locally
#' absent rings). E.g., `0` makes sense. If `NULL`, these values are left as NAs.
#' @param trim Logical vector indicating whether to trim off NA sequences at the beginning or end
#' of individual series. Default is TRUE. Will not remove missing rings that are represented with
#' NA (i.e., NA values within the series).
#'
#' @return A data.table/data.frame with 3 columns: 1) "year", 2) series.name (e.g., "series", "tree"
#' , etc.),3) dat.name (e.g., "rw", "rw.mm", "bai.mm", etc.)
#'
#' @details
#' \code{\link{rwl_longer}} does a bit more than convert to long format. It also gives the output
#' as a `data.table` (see \code{\link[data.table]}), which is essentially a special kind of
#' data.frame with some different characteristics. The `data.table` package enables streamlined and
#' very memory efficient code. If you're not ready to dive into `data.table`, use
#' \code{\link[base]{as.data.frame()}} to convert to a regular data.frame.
#'
#' \code{\link{rwl_longer}} also has built in data-continuity checking. For example, if your data
#' skips a year (this should almost never happen), like this: ...1999, 2000, 2002, 2003, ..., then
#' \code{\link{rwl_longer}} will add a row for the missing year and give an NA ring-width value.
#' A note of caution: any NA values added in this way will be replaced with `new.val.internal.na` if
#' you've given a value for that argument.
#'
#' @seealso [longer_rwl()]
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
#' ca533.long <- rwl_longer(rwl = ca533, dat.name = "rw.mm", trim = TRUE)
#' head(ca533.long)


rwl_longer <- function(rwl = NULL,
                       series.name = "series",
                       dat.name = "rw",
                       new.val.internal.na = NULL,
                       trim = TRUE) {
  . <- j_variable <- year <- NULL # 1st line of function - data.table

  #### Error catching ####
  # Basic checking 1st. Does not strictly require an "rwl" class.
  stopifnot("rwl is not an object of class 'rwl', 'data.frame', or 'matrix'" =
              data.class(rwl) %in% "rwl" ||
              is.data.frame(rwl) ||
              is.matrix(rwl)
  )

  if (is.matrix(rwl)) {
    rwl <- as.data.frame(rwl)
  }

  # More tests for rwl format beyond "rwl" class
  yrs <- suppressWarnings(as.integer(rownames(rwl))) # Coerce rownames to integer - should be years

  stopifnot(
    "rwl rownames are not years - are data in rwl-format (ie, wide-format)?" =
      # .row_names_info() gives negative values when rownames are auto-assigned, which means not rwl
      .row_names_info(rwl) > 0 &&
      !anyNA(yrs) # there cannot be NA yrs
  )

  stopifnot(
    "rwl rownames are not integer years - are data in rwl-format (ie, wide-format)?" =
      all(yrs == suppressWarnings(as.numeric(rownames(rwl))))
  )

  stopifnot(
    "rwl rownames are not consecutive years - are data in rwl-format (ie, wide-format)?" =
      # All years must be consecutive
      all(diff(yrs) == 1)
  )

  stopifnot(
    "rwl already has colname 'year' - must only have rownames as years and colnames as series IDs" =
      !any(tolower(colnames(rwl)) == "year")
  )

  # All the columns in an rwl must be numeric or double
  stopifnot(
    "rwl data must be all numeric" =
      all(vapply(rwl, is.numeric, logical(1)))
  )

  # Don't allow all NA series to pass through
  all_na <- colnames(rwl)[colSums(!is.na(rwl)) == 0]
  if (length(all_na) > 0) {
    stop("These series have no values (all NAs): ", paste(all_na, collapse = ", "), call. = FALSE)
  }

  # series.name arg
  stopifnot(
    "series.name can only be a character vector of length 1" =
      length(series.name) == 1 &&
      is.character(series.name)
  )

  # dat.name arg
  stopifnot(
    "dat.name can only be a character vector of length 1" =
      length(dat.name) == 1 &&
      is.character(dat.name)
  )

  # Prevent user from specifying series name & dat.name the same and neither can be "year"
  stopifnot("series.name and dat.name must differ and neither can be 'year' or 'Year'" =
              !anyDuplicated(c(series.name, dat.name, "year", "Year")))

  # new.val.internal.na arg
  stopifnot("new.val.internal.na must be NULL or a numeric vector of length = 1" =
              is.null(new.val.internal.na) ||
              (length(new.val.internal.na) == 1 &&
                 is.numeric(new.val.internal.na))
  )

  # trim arg
  stopifnot("trim must be a logical vector of length = 1 (TRUE or FALSE)" =
              isTRUE(trim) ||
              isFALSE(trim)
  )

  #### Process the data ####
  # Get the series names
  series.cols <- colnames(rwl)

  # Check for any replicated series names
  dups <- unique(series.cols[duplicated(series.cols)])

  if (length(dups) > 0) {
    stop("rwl has multiple series (columns) with the same name:\n",
         paste(dups, collapse = ", "), call. = FALSE)
  }

  # Convert to data.table - we can convert the rownames to a year column in the same step
  rwl.dt <- data.table::as.data.table(rwl, keep.rownames = "year")
  rwl.dt[, year := base::as.integer(year)] # Have to follow up with conversion to integer

  # Reshape to long format
  long.rwl <- data.table::melt(
    data = rwl.dt,
    id.vars = "year",
    variable.name = series.name,
    value.name = dat.name,
    variable.factor = FALSE
  )

  # Check that data are continuous - redundant and probably unnecessary, but good to have
  # Full year sequence per series, from each series' first year to its last
  full <- long.rwl[!is.na(dat), .(year = seq(min(year), max(year))),
                   by = series.name,
                   env = list(dat = dat.name)]

  # Which years are missing from which series
  trimmed <- long.rwl[full, on = c(series.name, "year")]
  missing <- trimmed[is.na(dat), env = list(dat = dat.name)]

  # Give a warning if there are any missing values
  if (nrow(missing) > 0) {
    msum <- missing[, .(years = paste(year, collapse = ",")), by = series.name]
    warning("Some series have internal NAs:\n",
            paste0(msum[[series.name]], ": ", msum$years, collapse = "\n"),
            if (!is.null(new.val.internal.na))
              paste0("\nThese will be replaced with ", new.val.internal.na),
            call. = FALSE)
  }

  # Fill the gaps: missing years get rows with NA rw
  long.rwl_full <- if (trim) {
    trimmed
  } else {
    long.rwl
  }
  data.table::setorderv(long.rwl_full, c(series.name, "year"))

  # Replace the internal NAs if new.val.internal.na is specified
  if (is.numeric(new.val.internal.na)) {
    long.rwl_full[, dat := modendro::replace_internal_NAs(dat, new.val = new.val.internal.na),
                  env = list(dat = dat.name), by = series.name]
  }

  # Return the resulting data.table
  long.rwl_full

} ## End of function
