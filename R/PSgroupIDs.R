#' Data.frame containing grouping data for tree-ring series (tree IDs) from Perkins & Swetnam 1996
#' whitebark pine (Pinus albicaulis) data from central Idaho, USA. This is an example data.frame
#' for demonstrating the use of grouped data in \code{\link{n_mon_corr}}. PSgroupIDs is fed to the
#' `group.IDs.df` argument, and `"site"` is fed to the `group.var` argument. To be used when
#' `idPRISMgroup` is used as the climate data.
#'
#' @format ## `PSgroupIDs`
#' A data.frame with 68 rows and 2 columns:
#' \describe{
#'   \item{series}{The series ID code - in this case the tree ID code}
#'   \item{site}{Site code corresponding to the Perkins & Swetnam study. SDP & UPS are combined as
#'   coordinates were the same in the ITRDB}
#'   ...
#' }
#'
#' @source <https://doi.org/10.1139/x26-241>
"PSgroupIDs"
