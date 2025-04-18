% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p2yrsL.R
\name{p2yrsL}
\alias{p2yrsL}
\title{Highlight interannual variability in tree ring series by computing ring width as the proportion
of the last 2 years, limited by a variance cut-off}
\usage{
p2yrsL(rwl = NULL, limit = 2.6)
}
\arguments{
\item{rwl}{A rwl object (read in by dplR's \code{\link[dplR]{read.rwl}}). Essentially a
data.frame with columns names as series IDs and years as rownames.}

\item{limit}{A numeric vector specifying the number of standard deviations beyond the mean to
limit (cut off) the resulting series. Default is 2.6, the same value used in
CDendro and CooRecorder.}
}
\value{
A rwl-style data.frame with the transformed ring width series.
}
\description{
This function implements Lars-Ake Larsson's (author of the popular CooRecorder and CDendro
programs) proportion of the last two years (limited) P2YrsL standardizing. This is calculated
as i = rw(c)/rw(c) + rw(c-1), with a cut off of 2.6x the SD for the whole series
(for details see: https://cdendro.se/wiki/index.php/Normalization).
The cut off is imposed to reduce (not eliminate!) the chances that you will get spurious
correlations because of the alignment of relatively tiny rings in a reference and the series
you are evaluating.

P2YrsL is similar in concept to computing AR residuals (a form of "prewhitening") - the default
in dplR. P2YrsL is the default in CDendro and CooRecorder. The essential element of both
approaches is that we remove the mid- to low-frequency variation and highlight the high-frequency
for cross-dating purposes only. The high-frequency variation facilitates checking both
statistical (correlations) and visual (squiggly line plots) correspondence of two series.
}
\details{
See https://cdendro.se/wiki/index.php/Normalization
}
\examples{
library(dplR)
data(ca533)
# before
ca533 |> spag.plot()
# after
ca533 |> p2yrsL() |> spag.plot()

}
\references{
Larsson & Larsson (2023) \emph{CDendro and CooRecorder programs of the CDendro package},
Cybis Elektronik & Data AB. https://www.cybis.se/forfun/dendro/index.htm
}
