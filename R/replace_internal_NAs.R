#' Replace NA values when they occur within a numeric series
#'
#' @description
#' Simple utility function to replace NA values when they occur within a numeric series and the
#' series has NA values on either end (or both ends). Chains of NA values are common at the
#' beginnings and ends of tree-ring series within rwl files, and by some conventions missing or
#' locally absent rings are marked with NA. This can cause problems for many analyses.
#'
#' This function is tailored to work with tree-ring series, but could do simple corrections to
#' climate data as well. Instead of 0, which would be non-sensical for climate data, you could use
#' mean of the whole series.
#'
#' Warning: for the default of new.val = 0, make sure that the NAs do indeed represent missing or
#' locally absent rings!
#'
#' @param x A numeric vector, e.g., a single tree-ring series from a rwl
#' @param new.val a numeric vector of length 1. The default is zero
#'
#' @export
#'
#' @examples
#' a.series <- c(NA,NA,NA,NA,1,2,1.5,1,0.5,NA,1,0.5,1,2,1.5,NA,NA,NA)
#' replace_internal_NAs(x = a.series, new.val = 0)


replace_internal_NAs <- function(x = NULL,
                                 new.val = 0) {

  stopifnot("x must be a numeric vector of length > 1" =
              is.numeric(x) &
              length(x) > 1
  )

  stopifnot("new.val must be a numeric vector of length == 1" =
              is.numeric(new.val) &
              length(new.val) == 1
  )

  # Get the NAs
  notNAs <- which(!is.na(x))
  # Get the "body" of x, which may include NAs within
  x.body <- x[min(notNAs):max(notNAs)]
  # Replace any NAs in the body
  x.body[is.na(x.body)] <- new.val
  # Put the fixed body back into the whole x
  x[min(notNAs):max(notNAs)] <- x.body

  return(x)
}
