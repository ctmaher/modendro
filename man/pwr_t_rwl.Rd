% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pwr_t_rwl.R
\name{pwr_t_rwl}
\alias{pwr_t_rwl}
\title{Find optimal power and apply Cook & Peters (1997) transformation on tree ring series}
\usage{
pwr_t_rwl(rwl = NULL, universal = FALSE, ID.group.substr = NULL)
}
\arguments{
\item{rwl}{A rwl object (read in by dplR's \code{\link[dplR]{read.rwl}}). Essentially a
data.frame with columns names as series IDs and years as rownames.}

\item{universal}{A logical vector indicating whether to compute a single "universal" optimal
power for all series within a group (if \code{TRUE} and arg \code{ID.group.substr} is unspecified, will
assume all series in the rwl). If \code{FALSE}, the function estimates a separate optimal power for
each individual series.}

\item{ID.group.substr}{A numeric vector of length 2 that determines a sub grouping for
calculating "universal" optimal power. Ulitmately passes to \code{\link[base]{substr}} as the
start (1st number) and stop (2nd number) args to split the series IDs in your rwl into groups
(performed in \code{\link{find_opt_pwr}}).}
}
\value{
A two-element list, 1 is the transformed series and 2 contains the messages about the
transformations
}
\description{
Function to apply power transformation to a tree ring series.
}
\details{
This function uses the estimated optimal power of transformation and a simple set of selection
criteria to transform a heteroscastic (variance changes with mean) tree ring series into a ±
homoscedastic one. The selection criteria are as follows: If the estimated power is ≤ 0.1, then
log10 transform. If greater than 1, don't transform the series. The set of estimated optimal
powers that is actually used (applied as: ring width ^ opt pwr) lies between 0.1 < opt pwr ≤ 1.

See \code{\link{find_opt_pwr}}, \code{\link{cp_detrend}},and Cook and Peters (1997) for more
details.

Some cleaning of the rwl data is performed before transformation: because log(0) is undefined,
0 ring width values are replaced with the minimum value possible given the resolution of the
data, which for tree ring data is typically 0.01 or 0.001 mm.
}
\examples{
library(dplR)
data("ca533")
pwr_t_rwl(rwl = ca533)
}
\references{
Cook, E. R., and Peters, K. (1997) Calculating unbiased tree-ring indices for the study of
climatic and environmental change.
\emph{The Holocene}, \strong{7}(3), 361-370.
}
\seealso{
\code{\link{find_opt_pwr}}, \code{\link{cp_detrend}}
}
