#' Find the optimal power for transforming a tree ring series to approximate homoscedasticity
#'
#' @description
#' Use a spread versus level method to estimate the optimal power of transformation for a tree
#' ring series
#'
#'
#' @details
#' This function estimates the optimal power of transformation via the local spread versus level
#' relationship of each series, where local spread is defined as the absolute value of the first
#' differences, S, (|rwt - rwt-1|) and the local level is the arithmetic mean of each pair of
#' adjacent values, M, (rwt + rwt-1)/2.
#' The spread versus level relationship is then modeled in a simple linear regression as
#' log10(S) ~ log10(M). The optimal power of transformation is then estimated as p = |1-slope| of
#' this regression model.
#' See \code{\link{cp_detrend}} and Cook & Peters (1997) for more details.
#'
#' Note that this is a very simple function that only implements the spread vs. level estimation
#' of an optimal power for making tree ring series ± homoscedastic. Operationally in the `modendro`
#'  package, the final determination of how to transform series is performed in the
#'  \code{\link{pwr_t_rwl}} function, per Cook & Peters (1997).
#'
#' The `universal` and `ID.group.substr` arguments control how many series are included in the
#' spread vs. level estimation, allowing for the same transformation to be performed across all
#' tree ring series within each group (you determine the group). `universal` is a logical switch
#' while `ID.group.substr` is a 2-element numeric that passes to the `start` & `stop` args of
#' \code{\link[base]{substr}} which will act on the series IDs (i.e., the colnames in `rwl`). For
#' example, if the first 3 characters of your series IDs refers to the site where the series were
#' collected, you would specify `ID.group.substr = c(1, 3)`. Or, say the 4th element of the series
#' IDs refers to the species collected, you could specify `ID.group.substr = c(4, 4)` to estimate
#' optimal power by species. This assumes your series IDs contain this information in a systematic
#' way. If not, you could also just create separate rwls for each of the groups you are interested
#' in, set `universal = TRUE` and leave `ID.group.substr` unspecified.
#' The essential idea of a "universal power transformation" comes from a GitHub issues post in the
#' dplR repository by Stefan Klesse (2020). My approach here only scratches the surface of the idea
#' presented by Klesse, but is an incremental step in that direction.
#'
#'
#' @param rwl A rwl object (read in by dplR's  \code{\link[dplR]{read.rwl}}). Essentially a
#' data.frame with columns names as series IDs and years as rownames.
#' @param universal A logical vector indicating whether to compute a single "universal" optimal
#' power for all series within a group (if `TRUE` and arg `ID.group.substr` is unspecified,
#' will assume all series in the rwl). If `FALSE`, the function estimates a separate optimal power
#' for each individual series.
#' @param ID.group.substr A numeric vector of length 2 that determines a sub grouping for
#' calculating "universal" optimal power. Passes to \code{\link[base]{substr}} as the start
#' (1st number) and stop (2nd number) args to split the series IDs in your rwl into groups.
#'
#' @return A named numeric vector with the series IDs (colnames) and the estimated optimal power
#' of transformation for each series -or- a numeric vector of length 1.
#'
#' @references
#' Cook, E. R., and Peters, K. (1997) Calculating unbiased tree-ring indices for the study of
#' climatic and environmental change.
#' \emph{The Holocene}, \strong{7}(3), 361-370.
#'
#' Klesse, S. (2020) "universal" or "signal-free" power transformation
#' https://github.com/OpenDendro/dplR/issues/8
#'
#' @seealso \code{\link{pwr_t_rwl}}, \code{\link{cp_detrend}}
#'
#' @export
#'
#' @examples
#' library(dplR)
#' data("ca533")
#' find_opt_pwr(rwl = ca533)
#'
#' # Universal power transformation for all series in the rwl - a sensible choice because
#' # these are all Pinus longaeva series from the same site in the White Mountains, California
#' find_opt_pwr(rwl = ca533, universal = TRUE)


find_opt_pwr <- function(rwl, universal = FALSE, ID.group.substr = NULL) {
  # Error catching
  stopifnot(
    "rwl is not an object of class 'rwl', 'data.frame', or 'matrix'" =
      data.class(rwl) %in% "rwl" |
      data.class(rwl) %in% "data.frame" |
      data.class(rwl) %in% "matrix"
  )

  stopifnot(
    "rwl has no rownames (must be years only) or no colnames (must be series IDs only)" =
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
    these_are_NA <-
      colnames(rwl)[which(apply(rwl, MARGIN = 2, FUN = \(x) all(is.na(x))) == TRUE)]
    stop("The following series have no values (all NAs): " ,
         paste(these_are_NA, collapse = ", "))
  }

  orig.IDs <- colnames(rwl) # original series names in original order

  # absolute value of 1st differences.
  # set up lagged vector for 1st differences
  rwl.lag <-
    rbind(rep(NA, ncol(rwl)), rwl[1:(nrow(rwl) - 1), , drop = FALSE])

  rwl.lag <- rwl.lag[, orig.IDs]

  # absolute value of 1st differences
  # Cook & Peters equate this as the standard deviation (it is not the same, but 1st differences
  # are what they use).
  diffs <- abs(rwl - rwl.lag)
  lmean <- (rwl + rwl.lag) / 2

  diffs <- diffs[, orig.IDs]
  lmean <- lmean[, orig.IDs]

  # Any zero values in diffs or means should be avoided
  diffs.ind <- which(diffs < 0.001, arr.ind = TRUE)
  lmean.ind <- which(lmean < 0.001, arr.ind = TRUE)

  diffs[diffs.ind] <- NA
  diffs[lmean.ind] <- NA
  lmean[diffs.ind] <- NA
  lmean[lmean.ind] <- NA

  # Make into long-format dfs using the modendro function rwl_longer

  diffs.long <- rwl_longer(diffs, series.name = "series",
                           dat.name = "diff", trim = TRUE, na.warn = FALSE)

  lmean.long <- rwl_longer(lmean, series.name = "series",
                           dat.name = "lmean", trim = TRUE, na.warn = FALSE)

  diff.lmean <- merge(diffs.long, lmean.long, by = c("year","series"))
  diff.lmean$series <- factor(diff.lmean$series, levels = orig.IDs)
  diff.lmean <- diff.lmean[order(diff.lmean$year, diff.lmean$series),]


  if (universal == FALSE) {

    # Split the df into a list - i.e., dfs for each series
    diff.lmean.list <- split(diff.lmean, f = diff.lmean$series)

    # Get the individual slopes for each series
    slopes <- lapply(diff.lmean.list, FUN = \(x) {
      lm(log10(x[,"diff"]) ~ log10(x[,"lmean"]), na.action = "na.exclude")$coefficients[2][[1]]
    })

    # Make sure the order matches the original rwl
    slopes <- slopes[orig.IDs] |> unlist()

    # Determine optimal power of transformation for each series by subtracting slope from 1
    abs(1 - slopes)

  } else {

    if (is.null(ID.group.substr)) {

      univ.slope <- lm(log10(diff.lmean[,"diff"]) ~ log10(diff.lmean[,"lmean"]),
                       na.action = "na.exclude")$coefficients[2][[1]]

      # Universal estimate with random slope & intercept for each series, and random intercept
      # for each year across all series
      # univ.slopes <- lmer(log10(diff) ~ log10(lmean) + (log10(lmean)|series) + (1|year),
      # data = diff.lmean)
      # The median of all random slopes
      # univ.slope <- coef(univ.slopes)$series[,"log10(lmean)"] |> median()

      # Alternate ideas
      # slopes4 <- lmer(log10(diff) ~ log10(lmean) + (1|series) + (1|year), data = diff.lmean)
      # slopes4
      #
      # slopes5 <- lmer(log10(diff) ~ log10(lmean) + (1|year), data = diff.lmean)
      # slopes5

      # Determine optimal power of transformation by subtracting universal slope from 1,
      # repeat for all series and return a named numeric
      univ.slopes <- rep(abs(1 - univ.slope), length(orig.IDs))
      names(univ.slopes) <- orig.IDs
      univ.slopes[orig.IDs]

    } else {
      diff.lmean$group <- substr(diff.lmean$series,
                                 start = ID.group.substr[[1]],
                                 stop = ID.group.substr[[2]])

      # Split the df into a list - i.e., dfs for each group
      diff.lmean.list <- split(diff.lmean, f = diff.lmean$group)

      # Get the slopes for each group of series
      slopes <- lapply(diff.lmean.list, FUN = \(x) {
        lm(log10(x[,"diff"]) ~ log10(x[,"lmean"]), na.action = "na.exclude")$coefficients[2][[1]]
      }) |>
        unlist() |>
        as.data.frame()

      colnames(slopes) <- "slope"
      slopes$group <- row.names(slopes)


      # Match up the group slopes to the appropriate series IDs
      diff.lmean.slope <- merge(diff.lmean, slopes, by = "group")

      # Extract just the slopes as a named numeric
      slopes.num <- diff.lmean.slope[diff.lmean.slope$series %in% orig.IDs, "slope"]
      names(slopes.num) <- diff.lmean.slope[diff.lmean.slope$series %in% orig.IDs, "series"]

      # Determine optimal power of transformation for each series by subtracting slope from 1
      abs(1 - slopes.num)[orig.IDs]# Make double sure the order matches the original rwl
    }
  }
} # End of function
