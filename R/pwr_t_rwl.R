#' Find optimal power and apply Cook & Peters (1997) transformation on tree ring series
#'
#' @description
#' Function to apply power transformation to a tree ring series.
#'
#' @details
#' This function uses the estimated optimal power of transformation and a simple set of selection
#' criteria to transform a heteroscastic (variance changes with mean) tree ring series into a ±
#' homoscedastic one. The selection criteria are as follows: If the estimated power is ≤ 0.1, then
#' log10 transform. If greater than 1, don't transform the series. The set of estimated optimal
#' powers that is actually used (applied as: ring width ^ opt pwr) lies between 0.1 < opt pwr ≤ 1.
#'
#' See \code{\link{find_opt_pwr}}, \code{\link{cp_detrend}},and Cook and Peters (1997) for more
#' details.
#'
#' Some cleaning of the rwl data is performed before transformation: because log(0) is undefined,
#' 0 ring width values are replaced with the minimum value possible given the resolution of the
#' data, which for tree ring data is typically 0.01 or 0.001 mm.
#'
#'
#' @param rwl A rwl object (read in by dplR's \code{\link[dplR]{read.rwl}}). Essentially a
#' data.frame with columns names as series IDs and years as rownames.
#' @param universal A logical vector indicating whether to compute a single "universal" optimal
#' power for all series within a group (if `TRUE` and arg `ID.group.substr` is unspecified, will
#' assume all series in the rwl). If `FALSE`, the function estimates a separate optimal power for
#' each individual series.
#' @param ID.group.substr A numeric vector of length 2 that determines a sub grouping for
#' calculating "universal" optimal power. Ulitmately passes to \code{\link[base]{substr}} as the
#' start (1st number) and stop (2nd number) args to split the series IDs in your rwl into groups
#' (performed in \code{\link{find_opt_pwr}}).
#'
#' @return A two-element list, 1 is the transformed series and 2 contains the messages about the
#' transformations
#'
#' @references
#' Cook, E. R., and Peters, K. (1997) Calculating unbiased tree-ring indices for the study of
#' climatic and environmental change.
#' \emph{The Holocene}, \strong{7}(3), 361-370.
#'
#' @seealso \code{\link{find_opt_pwr}}, \code{\link{cp_detrend}}
#'
#' @export
#'
#' @examples
#' library(dplR)
#' data("ca533")
#' pwr_t_rwl(rwl = ca533)

pwr_t_rwl <- function(rwl = NULL, universal = FALSE, ID.group.substr = NULL) {

  # Error catching
  stopifnot("rwl is not an object of class 'rwl', 'data.frame', or 'matrix'" =
              data.class(rwl) %in% "rwl" |
              data.class(rwl) %in% "data.frame" |
              data.class(rwl) %in% "matrix"
  )

  stopifnot("rwl has no rownames (must be years only) or no colnames (must be series IDs only)" =
              !is.null(rownames(rwl)) |
              !is.null(colnames(rwl))
  )

  stopifnot(
    "'universal' argument must be a logical vector (TRUE or FALSE)" =
      is.logical(universal)
  )

  stopifnot(
    "'ID.group.substr' argument must be NULL or a 2-element vector" =
      length(ID.group.substr) == 2 |
      is.null(ID.group.substr)
  )

  if (apply(rwl, MARGIN = 2, FUN = \(x) all(is.na(x))) |> any() == TRUE) {
    these_are_NA <- colnames(rwl)[which(apply(rwl, MARGIN = 2, FUN = \(x) all(is.na(x))) == TRUE)]
    stop("The following series have no values (all NAs): " , paste(these_are_NA, collapse = ", "))
  }

  orig.IDs <- colnames(rwl) # capture the original names in their original order

  # Replace ≤0 values in both rwls with the minimum possible non-zero value given the resolution of
  # the data
  # & make sure they match
  min.value <- ifelse(sapply(na.omit(unlist(rwl)),
                             FUN = \(x) nchar(sub(".", "", x, fixed=TRUE))) |>
                        max() <= 3,
                      0.01, 0.001)

  # Create an indexed array for each rwl that indicates very small values
  # These values (or 0s) don't play well with log10().
  rwl.ind <- which(rwl <= 0.001, arr.ind = TRUE)
  rwl0 <- rwl
  rwl0[rwl.ind] <- min.value

  # Get the estimated optimal power
  # if (add1 == TRUE) {
  #   optimal.pwr.t <- find_opt_pwr(rwl0, add1 = TRUE)
  # } else {
  optimal.pwr.t <- find_opt_pwr(rwl0, universal = universal, ID.group.substr = ID.group.substr)
  # }

  optimal.pwr.t <- optimal.pwr.t[orig.IDs]

  # if power is very small, then just log10 transform
  # If greater than 1, power transform with power = 1, same as the untransformed series
  # Pwr transform w/ optimal power. Cook & Peters 1997
  to.log <- which(optimal.pwr.t <= 0.1)
  no.trans <- which(optimal.pwr.t > 1)
  pwr.trans <- which(optimal.pwr.t > 0.1 & optimal.pwr.t <= 1)

  # Do the transformations - the no trans option is implicit
  rwl0[, to.log] <- log10(rwl0[,to.log])
  rwl0[, pwr.trans] <- mapply(FUN = function(x, y) {x^y},
                              x = rwl0[, pwr.trans], y = optimal.pwr.t[pwr.trans])

  pwr.t.df <- data.frame(series = names(optimal.pwr.t),
                         optimal.pwr = optimal.pwr.t)
  pwr.t.df$action <- ifelse(optimal.pwr.t <= 0.1, "log10 transformed",
                            ifelse(optimal.pwr.t > 1, "No transformation",
                                   "Power transformed"))

  # split back into a list
  pwr.t.list <- split(pwr.t.df, f = pwr.t.df$series)[orig.IDs]
  rwl0.out <- as.data.frame(rwl0)
  # Assign data class
  class(rwl0.out) <- c("rwl","data.frame")
  out.list <- list(rwl0.out, pwr.t.list, rwl)
  names(out.list) <- c("Transformed ring widths", "Transformation metadata", "Raw ring widths")
  out.list

} # end of function
