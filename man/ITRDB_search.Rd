% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ITRDB_search.R
\name{ITRDB_search}
\alias{ITRDB_search}
\title{Search the International Tree Ring Data Bank (ITRDB) for datasets by species and location}
\usage{
ITRDB_search(
  species = NULL,
  lon.range = NULL,
  lat.range = NULL,
  limit = 20,
  simplify = TRUE
)
}
\arguments{
\item{species}{A character vector specifying 1 or more species (4-letter codes; see examples
for a list of the NCEI species codes)}

\item{lon.range}{A numeric vector of length = 2 specifying the min and max (in that order)
bounds of longitude for the search, in decimal degrees.}

\item{lat.range}{A numeric vector of length = 2 specifying the min and max (in that order)
bounds of latitude for the search, in decimal degrees.}

\item{limit}{A numeric vector of length = 1 specifying a limit for how many datasets to return.
Default is an arbitrary 20. Set higher if your spatial bounds are reasonably small. You might
start with the default value and then increase gradually}

\item{simplify}{A logical vector that controls whether to return a simplified output
(recommended, & the default) or the raw output.}
}
\value{
A data.frame with the info for the datasets that meet the specified criteria.
}
\description{
This function harnesses the US National Oceanic & Atmospheric Administration's (NOAA) National
Center for Environmental Information (NCEI) API interface to find existing tree ring datasets
in the International Tree Ring Data Bank (ITRDB) based on species and location.

Finding a nearby tree ring dataset of the same (or similar) species as your study is an
important step in the cross-dating process for your samples.

This function does the search, but the downloads have to be done manually. A URL (& doi) are
included in the output data.frame.
}
\details{
The basic idea for this function derives from: https://www.ncei.noaa.gov/access/paleo-search/api
}
\examples{

# Get the species codes
spp.codes <- jsonlite::fromJSON(
"https://www.ncei.noaa.gov/access/paleo-search/study/params.json")

spp.codes <- c(spp.codes$species$NOAA[[1]], spp.codes$species$NOAA[[2]])
spp.codes <- spp.codes[order(spp.codes)]
spp.codes

# extract the abbreviated codes only
spp.codes.abv <- strsplit(spp.codes, split = ":") |>
  lapply(FUN = \(x) x[2]) |>
  unlist()

spp.codes.abv


# Find some spruce collections in northern Alaska
ITRDB_search(species = c("PCGL","PCMA"), # Picea glauca & Picea mariana
lon.range = c(-165, -140),
lat.range = c(64.5, 70),
limit = 10)

}
