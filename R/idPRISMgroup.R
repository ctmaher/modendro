#' PRISM climate time series 1895-1991 for the 4 sites (grouped data) where Perkins & Swetnam (1996)
#' collected whitebark pine (Pinus albicaulis) cores from in central Idaho, USA.
#'
#' Tree-level mean ring widths (mm) calculated from Perkins & Swetnam's (1996) original raw data.
#' Downloaded from the International Tree-Ring Data Bank (ITRDB).
#'
#' The tree series derive from 4 study sites in the mountains of central Idaho.
#'
#' @format ## `idPRISMgroup`
#' A data.frame with 3492 rows and 6 columns:
#' \describe{
#'   \item{Date}{Year-Month-Day in date format. This is monthly data, day is assumed to be 15th}
#'   \item{year}{Year}
#'   \item{month}{Month}
#'   \item{PPT.mm}{Total monthly precipitation in mm}
#'   \item{Tavg.C}{Average monthly temperature in degrees C}
#'   \item{site}{Site code corresponding to the Perkins & Swetnam study. SDP & UPS are combined as
#'   coordinates were the same in the ITRDB}
#'   ...
#' }
#'
#'
#' @source <https://www.prism.oregonstate.edu/explorer/>
"idPRISMgroup"
