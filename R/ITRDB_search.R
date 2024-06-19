#' Search the International Tree Ring Data Bank (ITRDB) for datasets by species and location
#'
#' @description
#' This function harnesses the US National Oceanic & Atmospheric Administration's (NOAA) National
#' Center for Environmental Information (NCEI) API interface to find existing tree ring datasets
#' in the International Tree Ring Data Bank (ITRDB) based on species and location.
#'
#' Finding a nearby tree ring dataset of the same (or similar) species as your study is an
#' important step in the cross-dating process for your samples.
#'
#' This function does the search, but the downloads have to be done manually. A URL (& doi) are
#' included in the output data.frame.
#'
#' @param species A character vector specifying 1 or more species (4-letter codes; see examples
#' for a list of the NCEI species codes)
#' @param lon.range A numeric vector of length = 2 specifying the min and max (in that order)
#' bounds of longitude for the search, in decimal degrees.
#' @param lat.range A numeric vector of length = 2 specifying the min and max (in that order)
#' bounds of latitude for the search, in decimal degrees.
#' @param limit A numeric vector of length = 1 specifying a limit for how many datasets to return.
#' Default is an arbitrary 20. Set higher if your spatial bounds are reasonably small. You might
#' start with the default value and then increase gradually
#' @param simplify A logical vector that controls whether to return a simplified output
#' (recommended, & the default) or the raw output.
#'
#'
#' @details
#' The basic idea for this function derives from: https://www.ncei.noaa.gov/access/paleo-search/api
#'
#' @return A data.frame with the info for the datasets that meet the specified criteria.
#'
#' @importFrom jsonlite fromJSON
#'
#' @export
#'
#' @examples
#'
#' # Get the species codes
#' spp.codes <- jsonlite::fromJSON(
#' "https://www.ncei.noaa.gov/access/paleo-search/study/params.json")
#'
#' spp.codes <- c(spp.codes$species$NOAA[[1]], spp.codes$species$NOAA[[2]])
#' spp.codes <- spp.codes[order(spp.codes)]
#' spp.codes
#'
#' # extract the abbreviated codes only
#' spp.codes.abv <- strsplit(spp.codes, split = ":") |>
#'   lapply(FUN = \(x) x[2]) |>
#'   unlist()
#'
#' spp.codes.abv
#'
#'
#' # Find some spruce collections in northern Alaska
#' ITRDB_search(species = c("PCGL","PCMA"), # Picea glauca & Picea mariana
#' lon.range = c(-165, -140),
#' lat.range = c(64.5, 70),
#' limit = 10)
#'

ITRDB_search <- function(species = NULL,
                         lon.range = NULL,
                         lat.range = NULL,
                         limit = 20,
                         simplify = TRUE) {

  stopifnot(
    "species arg must be a character vector of length 1 or more" =
      is.character(species) == TRUE &
      length(species) >= 1
  )

  stopifnot(
    "lon.range arg must be a numeric vector of length 2" =
      is.numeric(lon.range) == TRUE &
      length(lon.range) == 2
  )

  stopifnot(
    "lat.range arg must be a numeric vector of length 2" =
      is.numeric(lat.range) == TRUE &
      length(lon.range) == 2
  )

  stopifnot(
    "lon.range values must be between -180:180" =
      lon.range[1] >= -180 &
      lon.range[1] <= 180 &
      lon.range[2] >= -180 &
      lon.range[2] <= 180
  )

  stopifnot(
    "lat.range values must be between -80:80" =
      lat.range[1] >= -80 &
      lat.range[1] <= 80 &
      lat.range[2] >= -80 &
      lat.range[2] <= 80
  )

  stopifnot(
    "First element of lon.range arg must be the minimum value, second the maximum" =
      lon.range[1] < lon.range[2]
  )

  stopifnot(
    "First element of lat.range arg must be the minimum value, second the maximum" =
      lat.range[1] < lat.range[2]
  )


  stopifnot(
    "limit arg must be a numeric vector of length 1" =
      is.numeric(limit) == TRUE |
      is.null(limit)
  )

  stopifnot(
    "simplify must be a logical vector (TRUE or FALSE) of length 1" =
      is.logical(simplify) == TRUE &
      length(simplify) == 1
  )

  spp.codes <- jsonlite::fromJSON(
    "https://www.ncei.noaa.gov/access/paleo-search/study/params.json")

  spp.codes <- c(spp.codes$species$NOAA[[1]], spp.codes$species$NOAA[[2]])
  spp.codes <- spp.codes[order(spp.codes)]

  # extract the abbreviated codes only for comparison
  spp.codes.abv <- strsplit(spp.codes, split = ":") |>
    lapply(FUN = \(x) x[2]) |>
    unlist()

  stopifnot(
    "species code(s) don't match NCEI codes" =
      all(species %in% spp.codes.abv) == TRUE
  )

  ## specify API request, using values of arguments to paste together the parameters
  api.base = "https://www.ncei.noaa.gov/access/paleo-search/study/search.json?"

  req.params <- paste("dataPublisher=NOAA&dataTypeId=18",
                      paste0("minLat=", lat.range[1]),
                      paste0("maxLat=", lat.range[2]),
                      paste0("minLon=", lon.range[1]),
                      paste0("maxLon=", lon.range[2]),
                      ifelse(length(species) == 1,
                             paste0("species=", species),
                             paste0("species=",
                                    paste(strwrap(species, width = 4),
                                          collapse = "%7C"),
                                    "&speciesAndOr=or")
                      ),
                      ifelse(is.null(limit), "limit=1000",
                             paste0("limit=", limit)),
                      sep = "&"
  )

  req.str <- paste(api.base, req.params, sep="")

  json.data <- NULL
  try(
    json.data <- jsonlite::fromJSON(req.str)[[1]],
    silent = TRUE
  )

  stopifnot("No records found for specified search criteria" =
              !is.null(json.data)
  )

  json.data.abv <- json.data[,c("xmlId","NOAAStudyId","studyName",
                                "studyCode",
                                "earliestYearCE","mostRecentYearCE",
                                "doi","onlineResourceLink")]
  json.data.abv$lon <- lapply(json.data$site, FUN = \(x) {
    x$geo$geometry$coordinates[[1]][2]
  }) |> unlist() |> as.numeric()
  json.data.abv$lat <- lapply(json.data$site, FUN = \(x) {
    x$geo$geometry$coordinates[[1]][1]
  }) |> unlist() |> as.numeric()

  if (simplify == TRUE) {
    json.data.abv
  } else {
    json.data
  }
}
