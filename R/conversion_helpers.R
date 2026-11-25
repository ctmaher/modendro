### Tree-ring data format conversion helpers

## Internal helper function to ID rwl-format data
is_rwl_like <- function(x) {
  if (inherits(x, "rwl")) {
    return(TRUE) # one could give any data.frame "rwl" class, but why would you do this?
  }
  if (is.matrix(x)) { # matrices don't have row names
    x <- as.data.frame(x)
  }
  if (!is.data.frame(x) || .row_names_info(x) <= 0) {
    return(FALSE)  # data.tables always fail here because they don't have row names
  }
  !anyNA(suppressWarnings(as.numeric(rownames(x)))) && # all numeric rownames
    all(vapply(x, is.numeric, logical(1))) # all numeric values in columns
} ## End of function


## Internal helper function to detect format and to convert to long format data.tables
# Leverages rwl_longer() and longer_rwl()
rwl_dt_conv <- function(data,
                        series.name = "series",
                        dat.name = "rw",
                        new.val.internal.na = NULL) {
  if (is_rwl_like(data)) { # if data is an rwl
    return(rwl_longer(data, series.name = series.name, dat.name = dat.name,
                      new.val.internal.na = new.val.internal.na, trim = TRUE))
  }
  stopifnot("data must be an rwl, data.frame, or data.table" =
              is.data.frame(data)
  )
  data <- data.table::as.data.table(data)  # copies data.frames; longer_rwl() copies data.tables
  if (!"year" %in% names(data) &&
      "Year" %in% names(data)) {
    data.table::setnames(data, "Year", "year")
  }

  # Handle the read_pos convention of rw.mm
  if (!"rw" %in% names(data) &&
      "rw.mm" %in% names(data)) {
    data.table::setnames(data, "rw.mm", "rw")
  }

  # Check that crucial 3 columns are present
  miss <- setdiff(c(series.name, "year", dat.name), names(data))
  if (length(miss) > 0) {
    stop("data is missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  }

  # Round trip through longer_rwl and rwl_longer reuses all existing validation:
  # NAs, duplicate years, gap filling, internal-NA warning
  rwl_longer(longer_rwl(data, series.name = series.name, dat.name = dat.name),
             series.name = series.name, dat.name = dat.name,
             new.val.internal.na = new.val.internal.na, trim = TRUE)
  # Note that longer_rwl() simplifies data to just series names, years, and ring widths,
  # so any additional data is lost here
} ## End of function


# Wrapper function for rwl_dt_conv to call out where errors occur in a list. E.g., for when
# we use rwl_dt_conv in lapply, as in xd_check
with_context <- function(expr, label) {
  withCallingHandlers(
    expr,
    error = \(e) stop(label, ": ", conditionMessage(e), call. = FALSE),
    warning = \(w) {
      warning(label, ": ", conditionMessage(w), call. = FALSE)
      invokeRestart("muffleWarning")
    }
  )
} ## End of function

# E.g.:
# ref <- lapply(seq_along(ref), \(i) {
#   with_context(as_long_rw(ref[[i]]), sprintf("ref element %d", i))
# })
