#' PRISM climate time series 1895-1991 averaged over the 4 sites where Perkins & Swetnam (1996)
#' collected whitebark pine (Pinus albicaulis) cores from in central Idaho, USA
#'
#' Tree-level mean ring widths (mm) calculated from Perkins & Swetnam's (1996) original raw data.
#' Downloaded from the International Tree-Ring Data Bank (ITRDB).
#'
#' The tree series derive from 4 study sites in the mountains of central Idaho.
#'
#' @format ## `idPRISM`
#' A data.frame with 1164 rows and 4 columns:
#' \describe{
#'   \item{Date}{Year-Month-Day in date format. This is monthly data, day is assumed to be 15th}
#'   \item{year}{Year}
#'   \item{month}{Month}
#'   \item{PPT.mm}{Total monthly precipitation in mm}
#'   \item{Tavg.C}{Average monthly temperature in degrees C}
#'   ...
#' }
#'
#'
#' @source <https://www.prism.oregonstate.edu/explorer/>
"idPRISM"
