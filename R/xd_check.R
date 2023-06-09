#' Cross-dating checking using visual and statistical methods inspired by CooRecorder/CDendro and COFECHA.
#'
#' @description
#' This function implements Lars-Ake Larsson's (author of the popular CooRecorder and CDendro programs) proportion of the last two years
#' (limited) P2YrsL standardizing. This is calculated as i = rw(c)/rw(c) + rw(c-1),
#' with a cut off of 2.6x the SD for the whole series (for details see: https://cdendro.se/wiki/index.php/Normalization).
#' The cut off is imposed to reduce (not eliminate!) the chances that you will get spurious correlations because of the alignment of relatively tiny rings in a reference and the series you are evaluating.
#'
#' P2YrsL is similar in concept to computing AR residuals (aka "prewhitening") - the default in dplR. P2YrsL is the default in CDendro and CooRecorder.
#' The essential element of both approaches is that we remove the mid- to low-frequency variation
#' and highlight the high-frequency for cross-dating purposes only. The high-frequency variation facilitates checking both statistical (correlations) and visual (comparing the squiggly line plots of two series) correspondence of two series.
#'
#' @param rwl A rwl object (read in by dplR's `read.rwl()`). Essentially a data.frame with columns names as series IDs and years as rownames.
#' @param std.method A character vector specifying the standardization method to use. Choices are c("Ar","p2yrsL").
#'
#'
#' @details
#' See https://cdendro.se/wiki/index.php/Normalization
#'
#' @return A series of plots and text outputs?
#'
#' @references
#' Larsson & Larsson (2023) \emph{CDendro and CooRecorder programs of the CDendro package},
#'  Cybis Elektronik & Data AB. https://www.cybis.se/forfun/dendro/index.htm
#'
#' @import dplR
#'
#' @export
#'
#' @examples
#' library(dplR)
#' data(ca533)
#' # before
#' ca533[1000:1358,] |> spag.plot()
#' # after
#' ca533[1000:1358,] |> p2yrsL() |> spag.plot()
#'
