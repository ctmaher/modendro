#' Post oak (QUST) tree-ring width collection from Democrat Ridge, Missouri 1992 update
#' D. Stahle, G. Hawks, S. Sierzchula, N. Montagu
#'
#' @format ## `mo024`
#' An rwl data.frame with 373 rows and 61 columns (i.e., 61 tree-ring series spanning 373 years)
#'
#' @source <https://www.ncei.noaa.gov/access/paleo-search/study/4836>
"mo024"

#' Perkins & Swetnam whitebark pine ring width collection, four Idaho sites, 1996
#' D. Perkins, T. Swetnam
#' Four ITRDB collections are combined here and tree-level mean ring widths calculated.
#'
#' @format ## `ps96`
#' An rwl data.frame with 1267 rows and 68 columns (68 tree-ring series spanning 1267 years)
#'
#' @source Publication: <https://doi.org/10.1139/x26-241>
#' ITRDB collections:
#'  <https://www.ncei.noaa.gov/access/paleo-search/study/4108>
#'  <https://www.ncei.noaa.gov/access/paleo-search/study/4108>
#'  <https://www.ncei.noaa.gov/access/paleo-search/study/4110>
#'  <https://www.ncei.noaa.gov/access/paleo-search/study/4111>
"ps96"

#' Site grouping data for Perkins & Swetnam whitebark pine collection, four Idaho sites, 1996
#' D. Perkins, T. Swetnam
#' Four sites (TWP, SDP, UPS, & RRR) from the original study are condensed into three
#' (TWP, SDP_UPS, & RRR). The SDP & UPS sites colocated - gridded climate estimates would not be
#' different.
#'
#' @format ## `ps.groupIDs`
#' A data.frame with 68 rows and 2 columns
#' \describe{
#'   \item{tree}{Tree identifier}
#'   \item{site}{Site identifier}
#'   ...
#' }
#'
#' @source Publication: <https://doi.org/10.1139/x26-241>
"ps.groupIDs"

#' PRISM climate data to accompany Perkins & Swetnam whitebark pine collection
#' 1895-1991
#'
#' @format ## `idPRISM`
#' An rwl data.frame with 1164 rows and 5 columns:
#' \describe{
#'   \item{Date}{Date in %Y-%m-%d format}
#'   \item{year}{Year}
#'   \item{month}{Month}
#'   \item{PPT.mm}{Monthly total precipitation, in mm}
#'   \item{Tavg.C}{Monthly average daily temperature, in °C}
#'   ...
#' }
#'
#' @source <https://www.prism.oregonstate.edu/explorer/>
"idPRISM"

#' Site-specific PRISM climate data to accompany Perkins & Swetnam whitebark pine collection
#' 1895-1991
#'
#' @format ## `idPRISMgroup`
#' An rwl data.frame with 3492 rows and 6 columns:
#' \describe{
#'   \item{Date}{Date in %Y-%m-%d format}
#'   \item{year}{Year}
#'   \item{month}{Month}
#'   \item{PPT.mm}{Monthly total precipitation, in mm}
#'   \item{Tavg.C}{Monthly average daily temperature, in °C}
#'   \item{site}{Whitebark pine sampling location}
#'   ...
#' }
#'
#' @source <https://www.prism.oregonstate.edu/explorer/>
"idPRISMgroup"



