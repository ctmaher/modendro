#' Read in CooRecorder .pos files
#'
#' @description
#' \code{\link{read_pos}} reads in .pos files from CooRecorder (Cybis Elektronik & Data AB, Larsson
#' & Larsson). \code{\link{read_pos}} can handle a single file or a directory containing multiple
#' .pos files. \code{\link{read_pos}} uses package \code{\link[parallel]} to process chunks of pos
#' files simultaneously, and is pretty fast even with 1000s of files. Performance will vary
#' depending on your machine.
#'
#' The motivation of this function is to replicate some of the basic operations of CDendro in R.
#' This allows the user to avoid saving tree-ring collections using some of the arcane tree-ring
#' file formats which don't allow users to take full advantage of the great features of CooRecorder,
#' such as point labels, comments, etc. It also calculates the ring widths with greater precision
#' (they are not rounded by default) than what CDendro exports. Additionally, seasonwood widths are
#' saved in the same step as whole ring widths, eliminating the extra steps you would have to do in
#' CDendro.
#'
#' \code{\link{read_pos}} is also a great tool for basic error checking in .pos files. The raw
#' coordinates are part of the output, so you can use faceted plots in ggplot2 to visually check
#' point placement (see Examples below). This is particularly useful for checking pith location
#' estimates.
#'
#' The outputs of \code{\link{read_pos}} are long-format data.tables/data.frames, which are much
#' more versatile than the old rwl-format. If you need to export this data as a .rwl, however,
#' you can use \code{\link{longer_rwl}} and dplR's \code{\link[dplR]{write.rwl}}. Note that for
#' \code{\link[dplR]{write.rwl}} series name lengths are limited.
#'
#' @param path file path to a single .pos file or a directory that may contain several .pos files.
#'
#' @param default.OD optional - a numeric vector of length = 1 to give an outer date (OD) year for
#' files that did not have a date assigned in CooRecorder.
#'
#' @details If `path` is a directory the function will search all sub-directories for .pos files, thus
#' accommodating a range of directory structures. The main outputs are 1) a "Ring widths" data.frame
#' containing the ring widths (whole ring, late wood, and early wood, as applicable) and any
#' year-specific point labels 2) an "Attributes" data.frame containing the distance to pith, the
#' outer data, the inner-most date, the radius (sum of all ring widths plus the distance to pith),
#' and the comment for whole series, 3) a "Raw coordinates" data.frame containing the original
#' coordinates from the .pos file after they have been converted to numeric values, and 4) a
#' "Not read" data.frame which lists the files that were not read in and gives an error message
#' describing a likely reason.
#'
#' The `default.OD` argument is for circumventing a common "Not read" situation in which the user
#' did not assign an OD year to the series in CooRecorder - but be careful!! Sometimes CooRecorder
#' will still give a number for the OD, and this will not be caught by `default.OD`. Always check
#' the outputs by looking at plots or the data themselves. If you see weird years
#' assigned to the ring widths, then go back to CooRecorder and make sure those series have the
#' correct OD year assigned.
#'
#' \code{\link{read_pos}} contains several error catching heuristics that try to minimize direct
#' interventions from the user or try to clearly guide the user to a fix they can make in
#' CooRecorder.
#'
#' Because the output from read_pos includes the `"Raw coordinates"`, you can check for these
#' errors in R. This is much more efficient than going back into to CooRecorder to check multiple
#' .pos files. See the examples below for a workflow using faceted plots in ggplot2.
#'
#' A possible error in CooRecorder is that points can be saved out of order. Usually CooRecorder
#' will give you an "Erroneous point order" message, but it will save the file anyhow. This can
#' wreak havoc on determining ring widths. \code{\link{read_pos}} has a simple way of determining
#' if there are out-of-order points and it will warn you if it finds something. Since false
#' warnings are likely with this approach, \code{\link{read_pos}} will read in the files anyway
#' with warnings the console and in the `message` column of the `"Attributes"` data.frame in the
#' output. I recommend plotting the coordinates, as in the Examples below. You can then decide if
#' these warnings are warranted or not.
#'
#' Another possible error in CooRecorder is that the pith location gets jumbled around and is no
#' longer valid. This will make any calculations using distance to pith wrong, like age estimates
#' or basal area/basal area increment! Unfortunately this doesn't cause any error to arise in
#' CooRecorder itself. In attempt to catch these errors, \code{\link{read_pos}} uses a similar logic
#' as for the erroneous points - if the pith veers off in a different direction from the last set
#' of points \emph{or} if the distance to pith from the last point is greater than 50% of the range
#' of the long axis of coordinates, you'll get a warning in the console and in the `message` column
#' of the `"Attributes"` data.frame in the output. There is a constraint in that this check will
#' ignore pith distances within the range of normal ring width distances. This may produce some
#' false warnings, but will keep the actual errors to a manageable level. Assuming that the small
#' pith distances are close to accurate (i.e., the true pith distance is \emph{not} really large),
#' then these small errors should have a small affect on the accuracy of estimates using distance
#' to pith.
#'
#' \code{\link{read_pos}} handles crack marking points (i.e., gap points) differently than CDendro
#' does. CDendro measures the distance between two ring boundary points and then subtracts the
#' distance between two gap points or between one gap point and the earlier (in time) ring boundary
#' point. This method assumes that all points fall on the same straight line. This is not always
#' possible or realistic. To remedy this and allow more flexible use of crack marking,
#' \code{\link{read_pos}} instead measures the distance from the later ring boundary to the first
#' gap point, and then from the second gap point to the earlier ring boundary. It then ignores the
#' distance between the two gap points, or if there is only one gap point, it ignores the distance
#' between the gap point and the earlier ring boundary. This gives the same result as CDendro if all
#' points are in line, but gives a sensible result when they are not in line.
#'
#' In CooRecorder, it is possible to have a seasonwood ("W") point or a gap point at the very end of
#' a measurement series. Since these are nonsensical for determining distances between points, these
#' points are removed by \code{\link{read_pos}} automatically.
#'
#' \code{\link{read_pos}} warns you if it finds replicate .pos file names. File names are series
#' names, so you will need to fix this in your file system. The warning lists the replicated file
#' names and these show up in the error.message column of the Attributes data.frame in the output
#' list.
#'
#' For the forseeable future, \code{\link{read_pos}} will only work with .pos files from CooRecorder
#' 7.8 or greater - which is the earliest version (I think!) to include the actual pith coordinates
#' in the .pos files.
#'
#' @return A list containing 4 data.tables for ring widths, attributes, the original coordinates,
#' and a 2-column data.table of the files not read and the associated error message.
#'
#' @references
#' Larsson & Larsson (2023) \emph{CDendro and CooRecorder programs of the CDendro package},
#'  Cybis Elektronik & Data AB. https://www.cybis.se/forfun/dendro/index.htm
#'
#' @seealso \code{\link{longer_rwl}}, \code{\link{rwl_longer}}
#'
#' @import stats
#' @import data.table
#' @import parallel
#'
#' @export
#'
#' @examples
#'
#' library(ggplot2)
#' # Read in some example .pos files that show normal files and behavior on files with errors.
#' ex.pos <- read_pos(system.file("extdata", package = "modendro"))
#' # We get two erroneous point order warnings - one is real the other is a false positive.
#' # Below we can see the difference between the two.
#'
# Check the contents of the output list
#' names(ex.pos)
#'
#' # Take a look at the ring widths
#' ex.pos[["Ring widths"]] |> head()
#'
#' # Take a look at the attributes
#' ex.pos[["Attributes"]]
#'
#' # "Not read" gives you a data.frame of files that were not read in and potentially why
#' ex.pos[["Not read"]]
#'
#' # Check the coordinates - this is an efficient way to check point order or pith location errors
#' ggplot(ex.pos[["Raw coordinates"]], aes(x, y)) +
#'   geom_path() +
#'   geom_point(aes(color = type)) +
#'   facet_wrap(~series, ncol = 1, scales = "free")
#' # Note that one file truly had erroneous point order - signified by the jagged black line from
#' # geom_path (which plots points in the order it receives them).
#' # This file you would want to fix in CooRecorder
#'
#' # take a look at the ring widths - what you came here for
#' ggplot(ex.pos[["Ring widths"]], aes(year, rw.mm)) +
#'   geom_line() +
#'   facet_wrap(~series, ncol = 1, scales = "free")
#' # The true erroneous order file has invalid ring widths.

read_pos <- function(path = NULL,
                     default.OD = NULL) {
  current.opciones <- options()
  on.exit(options(current.opciones), add = TRUE)
  options(digits = 8)

  ## Error checking
  # Is path a character?
  stopifnot(
    "path argument must be a character vector (a file path or paths)" =
      is.character(path) |
      length(path) >= 1
  )

  if (!is.null(default.OD)) {
    stopifnot(
      "default.OD argument must be a numeric vector (a year)" =
        is.numeric(default.OD) &
        length(default.OD) == 1
    )
  }

  # Determine if path is a directory or not using list.files, then slim down to a list of .pos
  # files
  pos.files <- list.files(
    path,
    pattern = "\\.pos{1}",
    recursive = TRUE,
    include.dirs = FALSE,
    full.names = TRUE
  )

  if (length(pos.files) == 0 &&
      !length((grep("\\.pos{1}", path))) == 0) { # For single files
    pos.files <- path

    stopifnot("path does not lead to a .pos file" =
                file.exists(pos.files))
  } else { # For multiple files

    # There must be at least one .pos file in the path - this doesn't check the single files
    stopifnot("path must lead to at least one .pos file" =
                length(pos.files[grep("\\.pos{1}", pos.files)]) >= 1)

    # Check for replicate series
    filenames <- sapply(pos.files, FUN = \(f) {
      path.pieces <- strsplit(f, split = "/", fixed = TRUE)[[1]]
      path.pieces[[length(path.pieces)]]
    }) # This isolates the .pos files only
    fdt <- data.table::data.table(filenames = filenames, ind = 1:length(filenames))
    #fdt.agg <- aggregate(ind ~ filenames, data = fdt, FUN = \(f) length(f))
    fdt.agg <- fdt[, .(ind = .N), by = filenames]
    these.are.dups <- fdt.agg[ind > 1]

    if (nrow(these.are.dups) != 0) {
      # Give warning
      warning(paste0(
        "path contains replicate .pos files (also see error.message column in 'Attributes'):\n",
        paste(these.are.dups$filenames, collapse = "\n")),
        call. = FALSE, immediate. = TRUE)
      # We can flag the duplicated ones in the output too
    }
  }

  ## Run the pos parsing through all files.
  # Determine and set up for parallel process
  n.cores <- parallel::detectCores()
  if (is.na(n.cores)) {n.cores <- 1L}
  n.workers <- max(1L, min(n.cores - 1L, length(pos.files)))

  if (n.workers > 1L) {
    cl <- parallel::makeCluster(n.workers) # PSOCK by default -> Windows-safe
    on.exit(parallel::stopCluster(cl), add = TRUE) # add=TRUE: keeps your options() on.exit
    #system.time(
    out.list <- parallel::parLapplyLB(cl, pos.files, fun = \(f) {

      # For internal testing
      #f <- pos.files
      #f <- pos.files[[2]]
      .w <- character(0)
      res <- withCallingHandlers(
        {
          # wrap in tryCatch so we can handle the error
          tryCatch({
            # Read in the "raw" data
            raw.input <- scan(
              file = f,
              what = "character",
              sep = ";",
              quote = "",
              quiet = TRUE
            )

            # Control for this being empty - have had issues with this.

            # Setup check for CooRecorder version
            CR.ver.string <- raw.input[grep("CooRecorder=", raw.input)] |>
              strsplit(split = "=", fixed = TRUE)

            CR.ver <- strsplit(CR.ver.string[[1]][2], split = " ", fixed = TRUE)[[1]][1]
            CR.ver <- gsub(
              pattern = ".",
              replacement = "",
              CR.ver,
              fixed = TRUE
            )
            CR.ver <- ifelse(nchar(CR.ver) < 3, paste0(CR.ver, "0"), CR.ver) |> as.numeric()

            ## Get the outer date, if it exists. If not replace with default.OD
            if (length(grep("DATED", raw.input)) >= 1) {
              date.line <- raw.input[grep("DATED", raw.input)] |>
                strsplit(split = " ")

              OD.which <- which(!is.na(suppressWarnings(sapply(
                date.line[[1]], as.numeric
              ))))
              OD <- date.line[[1]][OD.which] |> as.numeric()
            } else {
              if (is.null(default.OD)) {
                OD <- NULL
              } else {
                OD <- default.OD
              }
            }


            # File has to be a DENDRO flavored .pos file - top line has this info
            # has to be CooRecorder version 7.8 or higher,
            # and must have a DATED line
            if (length(grep("DENDRO", raw.input[1])) == 1 &&
                CR.ver >= 780 &&
                !is.null(OD)) {
              # The crucial data are the seriesID, the outer date, the PithCoordinates,
              # comments if they exist,
              # & the ring boundary coordinates

              ## Get the series ID
              # The top line has the full file path, which includes the file name.
              # top.line <- raw.input[1] |> strsplit(split = "\\", fixed = TRUE)
              # The path may be complex, but the file name always ends in .pos
              # seriesID <- strsplit(top.line[[1]][length(top.line[[1]])], split = ".pos")[[1]][1]
              # Get the series ID from the file name
              f.split <- strsplit(f, split = "/", fixed = TRUE)[[1]]
              seriesID <- strsplit(f.split[length(f.split)], split = ".pos")[[1]]


              ## Get the DPI - this is used by CooRecorder to get the numeric x and y scales in mm
              DPI <- raw.input[grep("DPI", raw.input)] |>
                gsub(pattern = "#DPI ", replacement = "") |>
                as.numeric()

              ## Get the pith coordinates
              # If there are pith coords...
              if (length(grep("PithCoordinates", raw.input)) != 0) {
                pith.coord.line <- raw.input[grep("PithCoordinates", raw.input)] |>
                  strsplit(split = "=")
                # The coords are the after the "=", and are themselves separated by ","
                pith.coords <- pith.coord.line[[1]][2] |> strsplit(split = ",")
                # 1st value is x, second is y
                pith.coords.dt <- data.table::data.table(
                  series = seriesID,
                  x = as.numeric(pith.coords[[1]][1]),
                  y = as.numeric(pith.coords[[1]][2]),
                  label = NA,
                  type = "pith"
                )

                pith.coords.dt$type <- factor(pith.coords.dt$type,
                                              levels = c("reg", "multi1", "multi2", "W", "gap", "pith"))
              } else {
                # If no pith coords
                pith.coords.dt <- NULL
              }


              ## Get the comment if it exists.
              ## DATED / Written / PithCoordinates / CalcRadius / CooRecorder / licensedTo, etc. are all
              ## '#C ' lines here too, so the comment is the one '#C ' line that ISN'T metadata.
              # I might get an exception here if I'm missing a possible metadata type
              cC.text  <- sub("^#C ", "", grep("^#C ", raw.input, value = TRUE))
              meta.pat <- "^(DATED |Written=|PithCoordinates=|CalcRadius=|CooRecorder=|licensedTo=|YearsToPith=|Radius=|DistanceToPith=)"
              comment.text <- cC.text[!grepl(meta.pat, cC.text)]
              comment <- if (length(comment.text) >= 1) paste(comment.text, collapse = "; ") else NA_character_


              ## Get the point coordinates
              # The license info seems to be the 2nd-to-last line before the rest of the coordinates
              # There is a blank line in between
              # This is usually true but not always! Maybe I need more of a positive ID of the
              # coordinates.
              # Find the first line that contains just coordinates. That is, two stings coercible to
              # numeric that are separated by a comma. The key is that there is nothing else in the
              # string.
              # coord.start <- grep("licensedTo", raw.input) + 2

              # This method will not count gaps or season wood - but I'm pretty sure it is not
              # possible to have these as the first point.
              suppressWarnings(
                coord.start <- which(!is.na(sapply(strsplit(raw.input, split = ","),
                                                   FUN = \(x) as.numeric(x[[1]]))))[1]
              )

              # get all the coords separated from the header
              ring.bound.raw <- raw.input[coord.start:length(raw.input)]

              # Any ZERO points have to handled first
              suppressWarnings(
                ring.bound.dt <- data.table::data.table(
                  ring.bound.raw = ring.bound.raw,
                  zeros.rep = 1 + as.numeric(sub(".*?#ZERO", "", ring.bound.raw))
                )
              )
              ring.bound.dt$zeros.rep <- ifelse(is.na(ring.bound.dt$zeros.rep),
                                                1,
                                                ring.bound.dt$zeros.rep)
              ring.bound <- rep(ring.bound.raw, times = ring.bound.dt$zeros.rep)



              ## ring.bound is a character vector, the following turns it into a data.table of coordinates

              ## Vectorized point parse (replaces the per-point mapply + rbindlist)
              x.raw <- unname(ring.bound)
              n <- length(x.raw)

              # D-prefixed points are seasonwood ("W") / gaps; strip the leading D only on those
              is.D           <- substr(x.raw, 1L, 1L) == "D"
              base.str       <- x.raw
              base.str[is.D] <- sub("D", "", x.raw[is.D])          # == noD, but only where needed

              # Coordinate portion = everything before the first " #", trimmed so a trailing
              # space on single points doesn't read as a second pair
              coords.string <- trimws(sub(" #.*", "", base.str))

              # A multipoint now = a non-D coord string with whitespace *between* two pairs
              is.multi <- !is.D & grepl("\\s", coords.string)

              # --- labels & types -------------------------------------------------------
              # Non-D label: text after the first '#', else NA (cleaned to NA later at the "" step)
              labelA <- ifelse(grepl("#", x.raw), sub(".*?#", "", x.raw), NA_character_)

              # D label/type: text after first '#'; a second '#' splits it into type#label
              labelsD   <- sub(".*?#", "", base.str)               # base.str == noD on D rows
              hasSecond <- grepl("#", labelsD)
              typeD     <- ifelse(hasSecond, sub("#.*",   "", labelsD), labelsD)
              labelD    <- ifelse(hasSecond, sub(".*?#",  "", labelsD), NA_character_)
              typeD     <- gsub("%", "", typeD)                    # drop the '%' from gaps
              labelD    <- gsub("^#+|#+$", "", labelD)             # strip stray leading/trailing '#'

              label.el <- ifelse(is.D, labelD, labelA)             # one label per input element

              # --- coordinates ----------------------------------------------------------
              first.pair  <- sub("\\s.*", "", coords.string)                        # "x1,y1"
              second.pair <- ifelse(is.multi, sub("^\\S+\\s+", "", coords.string),  # "x2,y2"
                                    NA_character_)

              x1 <- as.numeric(sub(",.*$",    "", first.pair))
              y1 <- as.numeric(sub("^[^,]*,", "", first.pair))
              x2 <- as.numeric(sub(",.*$",    "", second.pair))    # NA where not multi
              y2 <- as.numeric(sub("^[^,]*,", "", second.pair))

              # --- expand (multipoints -> two rows) and assemble ------------------------
              n.rows <- ifelse(is.multi, 2L, 1L)
              oi     <- rep.int(seq_len(n), n.rows)                # output row -> input element
              slot   <- sequence(n.rows)                           # 1 for singles; 1,2 within a multi

              ring.bound.dt <- data.table::data.table(
                series = seriesID,
                x      = data.table::fifelse(slot == 2L, x2[oi], x1[oi]),
                y      = data.table::fifelse(slot == 2L, y2[oi], y1[oi]),
                label  = label.el[oi],
                type   = data.table::fifelse(
                  is.D[oi], typeD[oi],
                  data.table::fifelse(is.multi[oi],
                                      data.table::fifelse(slot == 1L, "multi1", "multi2"),
                                      "reg"))
              )

              # Files that end in a gap or a W should be trimmed - I think this is rare but is possible
              # Really what we need is that the file ends in a reg point or a multi2 point
              if(!(ring.bound.dt$type[nrow(ring.bound.dt)] %in% c("reg","multi2"))) {
                ring.bound.dt <- ring.bound.dt[1:(nrow(ring.bound.dt) - 1),]
              }

              # Check it twice!
              if(!(ring.bound.dt$type[nrow(ring.bound.dt)] %in% c("reg","multi2"))) {
                ring.bound.dt <- ring.bound.dt[1:(nrow(ring.bound.dt) - 1),]
              }

              # Add the pith coordinates.
              all.coords <- rbind(ring.bound.dt[, c("series","x","y","label","type")], pith.coords.dt)

              # Sometimes the files have "Erroneous order" messages, (points are out of out of order).
              # We can order the points by their x or y position. Leave out pith for now if it's there -
              # it always should go at the end.
              # The following several lines fix the order if it is broken

              # Each set of coords needs a unique identifier
              #all.coords$orig.index <- 1:nrow(all.coords)

              # ggplot(all.coords) +
              #   geom_path(aes(x, y), inherit.aes = F) +
              #   geom_point(aes(x, y, color = type), inherit.aes = F) +
              #   coord_fixed()

              # A fairly simple way to determine if any points are out of order is to see if there are any
              # directional changes in the differences of BOTH the axis. 1 at a time is okay and
              # plausible.
              check.diffs <- all.coords
              check.diffs$x.diff <- c(NA, diff(check.diffs$x))
              check.diffs$y.diff <- c(NA, diff(check.diffs$y))

              # The top row has NA diffs and must be removed
              check.diffs <- check.diffs[-1,]

              # Define the direction of movement
              check.diffs$x.dir <- ifelse(check.diffs$x.diff < 0,
                                          "neg", "pos")
              check.diffs$y.dir <- ifelse(check.diffs$y.diff < 0,
                                          "neg", "pos")

              # 0 diffs should take on the neg or pos label of the preceding point for direction
              # Must also account for 0 diffs at the top of - take the next one instead of the previous
              # one
              for (i in seq_along(check.diffs$x.dir)) {
                if (check.diffs$x.diff[i] == 0) {
                  if (length(check.diffs$x.dir[i - 1]) == 0) {
                    check.diffs$x.dir[i] <- check.diffs$x.dir[i + 1]
                  } else {
                    check.diffs$x.dir[i] <- check.diffs$x.dir[i - 1]
                  }
                }
              }

              for (i in seq_along(check.diffs$y.dir)) {
                if (check.diffs$y.diff[i] == 0) {
                  if (length(check.diffs$x.dir[i - 1]) == 0) {
                    check.diffs$y.dir[i] <- check.diffs$y.dir[i + 1]
                  } else {
                    check.diffs$y.dir[i] <- check.diffs$y.dir[i - 1]
                  }
                }
              }

              # Have to exclude the multi2 diffs, as those will represent the diff between multi points.
              # This assumes that the multi points are in order - they will be with respect to each other,
              # because I parsed them.
              check.diffs <- check.diffs[!(check.diffs$type %in% "multi2"),]

              # Find the prevailing direction of each axis
              # This is imperfect in that there can be multiple "prevailing" directions!
              # Tree-ring series may twist in such a way that the short axis can have two prevailing
              # directions. There are a lot of observations in these cases though.
              # I err on the side of producing false negatives rather than missing true negatives.
              # Exclude the pith at this point
              x.prevailing <- aggregate(series ~ x.dir,
                                        data = check.diffs[!(check.diffs$type %in% "pith"),],
                                        length)
              y.prevailing <- aggregate(series ~ y.dir,
                                        data = check.diffs[!(check.diffs$type %in% "pith"),],
                                        length)

              check.diffs$x.head <- ifelse(check.diffs$x.dir %in%
                                             x.prevailing$x.dir[which.max(x.prevailing$series)],
                                           "norm", "div")
              check.diffs$y.head <- ifelse(check.diffs$y.dir %in%
                                             y.prevailing$y.dir[which.max(y.prevailing$series)],
                                           "norm", "div")

              # If we have any observations where BOTH heading components are divergent, then we may have
              # erroneous point order. Read these in anyway, but give a warning and a tag.



              # Set up the null data.table
              error.dt <- data.table::data.table(file = f, message = NA)

              # All points but pith
              if (any(check.diffs$x.head[!(check.diffs$type %in% "pith")] %in% "div" &
                      check.diffs$y.head[!(check.diffs$type %in% "pith")] %in% "div")) {

                warning(paste0("Check coordinates for ", unique(check.diffs$series), ".pos - ",
                               "possible erroneous point order"),
                        call. = FALSE, immediate. = TRUE)

                error.dt <- data.table::data.table(file = f,
                                                   message = paste0("Check coordinates - ",
                                                                    "possible erroneous point order"))

              }

              # If pith.coords.dt exists,
              # check if pith direction is the same as last set of points,
              # if the distance is really large,
              # or if the pith lies within the range of the long axis.
              if (!is.null(pith.coords.dt)) {
                # Give a separate warning if the pith direction is divergent & based on different
                # criteria: only one heading has to be divergent or if the long axis dist is more
                # than 50% of the range of the other points
                long.axis <- ifelse(which.max(c(
                  abs(diff(range(check.diffs$x[!(check.diffs$type %in% "pith")]))),
                  abs(diff(range(check.diffs$y[!(check.diffs$type %in% "pith")]))))
                ) == 1,
                "x", "y")

                short.axis <- ifelse(which.min(c(
                  abs(diff(range(check.diffs$x[!(check.diffs$type %in% "pith")]))),
                  abs(diff(range(check.diffs$y[!(check.diffs$type %in% "pith")]))))
                ) == 1,
                "x", "y")

                # long.axis.range <- range(check.diffs[!(check.diffs$type %in% "pith"),
                #                                      long.axis])
                long.axis.range <- range(check.diffs[!(check.diffs$type %in% "pith"), ..long.axis])
                long.axis.range.diff <- abs(diff(long.axis.range))

                max.diff <- abs(c(check.diffs$x.diff[!(check.diffs$type %in% "pith")],
                                  check.diffs$y.diff[!(check.diffs$type %in% "pith")])) |> max()

                short.axis.dir <- paste0(short.axis, ".dir")

                if (
                  # pith diffs are larger than max diff for ring widths - ie, ignore the small pith diffs
                  (max(abs(c(check.diffs$x.diff[check.diffs$type %in% "pith"],
                             check.diffs$y.diff[check.diffs$type %in% "pith"]))) >
                   max.diff)
                  &&
                  # direction to pith is different than the last two points
                  !(check.diffs[check.diffs$type %in% "pith", ..short.axis.dir] %in%
                    check.diffs[(nrow(check.diffs) - 1), ..short.axis.dir])
                  ||
                  # max pith dist is larger than half of the long.axis.range.diff
                  max(abs(c(check.diffs$x.diff[check.diffs$type %in% "pith"],
                            check.diffs$y.diff[check.diffs$type %in% "pith"]))) >=
                  0.5*long.axis.range.diff
                ) {

                  warning(paste0("Check pith location for ", unique(check.diffs$series), ".pos - ",
                                 "possible pith location error"),
                          call. = FALSE, immediate. = TRUE)

                  if (is.na(error.dt$message)) {
                    error.dt$message <- paste0("Check pith location - ",
                                               "possible pith location error")
                  } else {
                    error.dt$message <- paste(error.dt$message,
                                              paste0("Check pith location - ",
                                                     "possible pith location error"),
                                              sep = "; ")
                  }
                }
              }

              # The coordinates can look inverted vertically relative to what I see in CooRecorder.
              # This doesn't matter for the distances

              # Clean up the label
              all.coords$label <- ifelse(all.coords$label %in% "", NA, all.coords$label)

              # type as a factor
              all.coords$type <- factor(all.coords$type,
                                        levels = c("reg", "multi1", "multi2", "W", "gap", "pith"))

              # Calculate the distances in mm using Pythagorean Theorem
              all.coords$x.dist <- c(NA, diff(all.coords$x, lag = 1, differences = 1))
              all.coords$y.dist <- c(NA, diff(all.coords$y, lag = 1, differences = 1))
              all.coords$dist.mm <- sqrt(all.coords$x.dist^2 + all.coords$y.dist^2)

              ## Assigning years based on OD and the order of points
              # If there is season wood, then we need to add the two distances together for each year.
              # If there is a gap, I need to subtract (or skip) the distance indicated by the gap.
              # Gaps can be indicated by 2 points (we skip the distance between them) or by a single point
              # (we skip the distance between a regular point and the gap point - also would apply to
              # "W" points).
              # The order of points is the only thing that indicates the year each point is
              # associated with.
              #all.coords$year1 <- NA
              all.coords$year <- NA
              all.coords$year[all.coords$type %in%
                                c("reg", "multi1")] <- seq(from = OD,
                                                           to = (OD -
                                                                   (nrow(all.coords[all.coords$type
                                                                                    %in%
                                                                                      c("reg", "multi1"),
                                                                   ]) - 1)),
                                                           by = -1)


              for (i in seq_along(all.coords$year)) {
                if (is.na(all.coords$year[i])) {
                  all.coords$year[i] <- all.coords$year[i - 1]
                }
              }

              # Year needs to be shifted by 1 to match dist.mm, which are ahead by 1
              all.coords$year <- c(NA, all.coords$year[1:(nrow(all.coords) - 1)])


              all.coords$year[all.coords$type %in% "pith"] <- NA

              # Gaps need some special handling - two gaps in a row need new labels - gap1 & gap2.
              # We will ignore gap2 distances, but keep gap1. Single gaps will stay as "gap", and we will
              # ignore those. This means that any dist AFTER a single gap point will be deleted.
              # I'll need to label these points too.

              # Find where gaps are and if they are sequential - run length encoding
              rle.vec <- rle(as.character(all.coords$type))

              # Initialize an empty result vector
              result <- as.character(all.coords$type)

              # Keep track of the position in the original vector
              pos <- 1

              for (i in seq_along(rle.vec$lengths)) {
                # Check for consecutive "gap" points
                if (rle.vec$values[i] == "gap" && rle.vec$lengths[i] > 1) {
                  # Create the relabeling sequence
                  labels <- rep(c("gap1", "gap2"), length.out = rle.vec$lengths[i])
                  # Assign the new labels
                  result[pos:(pos + rle.vec$lengths[i] - 1)] <- labels
                }
                # Also alter the labels of the points following single gaps
                if (rle.vec$values[i] == "gap" && rle.vec$lengths[i] == 1) {
                  result[pos + rle.vec$lengths[i]] <- paste0("gap.", result[pos + rle.vec$lengths[i]])
                }

                # Move to the next group
                pos <- pos + rle.vec$lengths[i]
              }

              all.coords$new.type <- result


              # Whole ring width is now the sum of dist.mm within each year, excluding multi2 and gaps
              whole.ring.widths <- all.coords[
                !(new.type %in% c("multi2", "gap2", "pith")) &
                  !(substr(new.type, 1L, 4L) %in% "gap.") &
                  !is.na(year) & !is.na(dist.mm),          # replicates na.action = na.omit
                .(rw.mm = sum(dist.mm)),                    # == cbind(rw.mm = dist.mm) ... FUN = sum
                keyby = year                               # ascending, matching aggregate's sort
              ]

              # Add the seriesID
              whole.ring.widths$series <- seriesID

              # ggplot(whole.ring.widths, aes(year, rw.mm)) +
              #   geom_line()

              ## Get seasonwood widths if they exist
              # These already exist actually, but need to be labeled properly

              if (any(all.coords$type %in% "W")) {
                #
                seas.wood.widths1 <- all.coords[!(all.coords$new.type %in%
                                                    c("multi2", "gap2", "pith")), ]
                # Do some error catching here - check for multiple W points per year. If this happens,
                # skip the EW/LW points and give the user a message about potential errors
                # per.year.check <- aggregate(type ~ year,
                #                             data = seas.wood.widths1,
                #                             FUN = length,
                #                             drop = FALSE)

                per.year.check <- seas.wood.widths1[!is.na(year), .(type = .N), keyby = year]

                if (any(per.year.check$type > 2)) {
                  paste0("Possibly multiple seasonwood boundaries detected in ",
                         per.year.check[per.year.check$type > 2, year],". ",
                         "Seasonwood data not read.")
                  whole.ring.widths$lw.mm <- NA
                  whole.ring.widths$ew.mm <- NA
                } else {

                  # Split by year
                  seas.wood.widths <- lapply(split(seas.wood.widths1, f = seas.wood.widths1$year),
                                             FUN = \(this.year) {
                                               if (any(this.year$type %in% "W")) {
                                                 # have to deal with the gaps if they exist
                                                 if (any(this.year$type %in% "gap")) {
                                                   for (i in 1:nrow(this.year)) {
                                                     if (this.year$type[i] %in% "gap" &&
                                                         i < nrow(this.year)) {
                                                       # Add "gap" dist to the next row
                                                       this.year$dist.mm[i + 1] <-
                                                         this.year$dist.mm[i + 1] +
                                                         this.year$dist.mm[i]
                                                     }
                                                   }
                                                 }
                                                 this.year <- this.year[!(this.year$type %in% "gap"), ]
                                                 this.year$wood.portion <- c("LW", "EW")
                                               } else {
                                                 this.year$wood.portion <- "WR"
                                               }
                                               this.year

                                             }) |> data.table::rbindlist()


                  # Add the seasonwood points to the whole ring widths
                  lw <- seas.wood.widths[seas.wood.widths$wood.portion %in% "LW", c("year", "dist.mm")]
                  colnames(lw)[colnames(lw) %in% "dist.mm"] <- "lw.mm"

                  ew <- seas.wood.widths[seas.wood.widths$wood.portion %in% "EW", c("year", "dist.mm")]
                  colnames(ew)[colnames(ew) %in% "dist.mm"] <- "ew.mm"

                  # Merge together and account for years that don't have seas wood with all = TRUE
                  whole.ring.widths <- Reduce(f = \(x, y) base::merge(x, y, by = "year", all = TRUE),
                                              list(whole.ring.widths, lw, ew))

                }
              } else {
                whole.ring.widths$lw.mm <- NA
                whole.ring.widths$ew.mm <- NA

              }


              year.labels <- all.coords[
                !is.na(year),
                .(label = {
                  lab <- label[!is.na(label)]
                  if (length(lab)) paste(as.character(lab), collapse = "; ") else NA_character_
                }),
                keyby = year
              ]

              ## Make the main output data.tables
              whole.ring.widths1 <- whole.ring.widths[year.labels, on = "year", nomatch = NULL]

              # ggplot(whole.ring.widths1, aes(year, rw.mm)) +
              #   geom_line()

              attributes <- data.table::data.table(
                series = seriesID,
                img.DPI = NA,
                d2pith.mm = ifelse(is.null(pith.coords.dt),
                                   NA,
                                   all.coords$dist.mm[all.coords$type %in% "pith"]),
                out.date = OD,
                in.date = min(whole.ring.widths1$year, na.rm = TRUE),
                total.rw.mm = sum(all.coords$dist.mm[!(all.coords$type %in% "multi2")], na.rm = TRUE)
              )


              error.message <- error.dt$message

              if (!is.na(DPI)) {
                attributes$img.DPI <- DPI
                DPI.message <- NA
                if (DPI < 600) {
                  DPI.message <- "DPI is suspiciously low!"
                  warning(paste0("DPI is low for ",
                                 attributes$series,
                                 ".pos. Ring width measurements will ",
                                 "be incorrect if DPI is wrong."),
                          call. = FALSE, immediate. = TRUE)
                }
                error.message <- paste(error.dt$message, DPI.message, sep = "; ")
              }

              attributes$radius.mm <- ifelse(is.null(pith.coords.dt),
                                             NA,
                                             attributes$total.rw.mm + attributes$d2pith.mm)
              attributes$comment <- comment

              # Clean up the error message and attach it to attributes
              attributes$error.message <- ifelse(error.message %in% "NA; NA", NA, error.message)

              # Order the attributes by series
              attributes <- attributes[order(attributes$series, decreasing = FALSE),]

              # Clean up and order the ring widths
              rw <- whole.ring.widths1[, c("series", "year", "rw.mm",
                                           "ew.mm", "lw.mm", "label")]

              rw <- rw[order(rw$series, rw$year, decreasing = FALSE),]

              # Assemble the output list
              this.series.list <- list(
                "Ring widths" = rw,
                "Attributes" = attributes,
                "Raw coordinates" = all.coords[, c("series", "year", "x", "y",
                                                   "dist.mm", "type", "label")]
              )

              this.series.list

            } else {
              # file doesn't meet criteria above for DENDRO files, CooRecorder version, or is not
              # dated.
              # Get specific here for the individual file
              if (length(grep("DENDRO", raw.input[1])) != 1) {
                tbdr.dt <- data.table::data.table(file = f,
                                                  message = "This file could not be identified as a DENDRO file (another data type?)")
              }

              if (CR.ver < 780) {
                tbdr.dt <- data.table::data.table(file = f,
                                                  message = enc2utf8(".pos files must be from CooRecorder \u2265 7.8 (update CooRecorder & resave file)"))
              } else {

                if (is.null(OD)) {
                  tbdr.dt <- data.table::data.table(file = f,
                                                    message = "This file was not dated (no outer year assigned in CooRecorder)")
                }}

              tbdr.dt
            }
          }, error = function(e) {
            # message("Error: Unknown problem with .pos file")
            # Return tbdr.dt in case of an error
            tbdr.dt <- data.table::data.table(file = f,
                                              message = "Unknown problem with .pos file. Check in CooRecorder.")
            tbdr.dt
          }) # end of tryCatch
        },
        warning = function(w) {
          .w[[length(.w) + 1L]] <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        }
      )
      attr(res, "read_pos_warnings") <- .w
      res
    })
    #)
    for (x in out.list) for (msg in attr(x, "read_pos_warnings")) warning(msg, call. = FALSE, immediate. = TRUE)
  } else {
    #system.time(
    out.list <- lapply(pos.files, FUN = \(f) {

      # For internal testing
      #f <- pos.files
      #f <- pos.files[[2]]

      # wrap in tryCatch so we can handle the error
      tryCatch({
        # Read in the "raw" data
        raw.input <- scan(
          file = f,
          what = "character",
          sep = ";",
          quote = "",
          quiet = TRUE
        )

        # Control for this being empty - have had issues with this.

        # Setup check for CooRecorder version
        CR.ver.string <- raw.input[grep("CooRecorder=", raw.input)] |>
          strsplit(split = "=", fixed = TRUE)

        CR.ver <- strsplit(CR.ver.string[[1]][2], split = " ", fixed = TRUE)[[1]][1]
        CR.ver <- gsub(
          pattern = ".",
          replacement = "",
          CR.ver,
          fixed = TRUE
        )
        CR.ver <- ifelse(nchar(CR.ver) < 3, paste0(CR.ver, "0"), CR.ver) |> as.numeric()

        ## Get the outer date, if it exists. If not replace with default.OD
        if (length(grep("DATED", raw.input)) >= 1) {
          date.line <- raw.input[grep("DATED", raw.input)] |>
            strsplit(split = " ")

          OD.which <- which(!is.na(suppressWarnings(sapply(
            date.line[[1]], as.numeric
          ))))
          OD <- date.line[[1]][OD.which] |> as.numeric()
        } else {
          if (is.null(default.OD)) {
            OD <- NULL
          } else {
            OD <- default.OD
          }
        }


        # File has to be a DENDRO flavored .pos file - top line has this info
        # has to be CooRecorder version 7.8 or higher,
        # and must have a DATED line
        if (length(grep("DENDRO", raw.input[1])) == 1 &&
            CR.ver >= 780 &&
            !is.null(OD)) {
          # The crucial data are the seriesID, the outer date, the PithCoordinates,
          # comments if they exist,
          # & the ring boundary coordinates

          ## Get the series ID
          # The top line has the full file path, which includes the file name.
          # top.line <- raw.input[1] |> strsplit(split = "\\", fixed = TRUE)
          # The path may be complex, but the file name always ends in .pos
          # seriesID <- strsplit(top.line[[1]][length(top.line[[1]])], split = ".pos")[[1]][1]
          # Get the series ID from the file name
          f.split <- strsplit(f, split = "/", fixed = TRUE)[[1]]
          seriesID <- strsplit(f.split[length(f.split)], split = ".pos")[[1]]


          ## Get the DPI - this is used by CooRecorder to get the numeric x and y scales in mm
          DPI <- raw.input[grep("DPI", raw.input)] |>
            gsub(pattern = "#DPI ", replacement = "") |>
            as.numeric()

          ## Get the pith coordinates
          # If there are pith coords...
          if (length(grep("PithCoordinates", raw.input)) != 0) {
            pith.coord.line <- raw.input[grep("PithCoordinates", raw.input)] |>
              strsplit(split = "=")
            # The coords are the after the "=", and are themselves separated by ","
            pith.coords <- pith.coord.line[[1]][2] |> strsplit(split = ",")
            # 1st value is x, second is y
            pith.coords.dt <- data.table::data.table(
              series = seriesID,
              x = as.numeric(pith.coords[[1]][1]),
              y = as.numeric(pith.coords[[1]][2]),
              label = NA,
              type = "pith"
            )

            pith.coords.dt$type <- factor(pith.coords.dt$type,
                                          levels = c("reg", "multi1", "multi2", "W", "gap", "pith"))
          } else {
            # If no pith coords
            pith.coords.dt <- NULL
          }


          ## Get the comment if it exists.
          ## DATED / Written / PithCoordinates / CalcRadius / CooRecorder / licensedTo, etc. are all
          ## '#C ' lines here too, so the comment is the one '#C ' line that ISN'T metadata.
          # I might get an exception here if I'm missing a possible metadata type
          cC.text  <- sub("^#C ", "", grep("^#C ", raw.input, value = TRUE))
          meta.pat <- "^(DATED |Written=|PithCoordinates=|CalcRadius=|CooRecorder=|licensedTo=|YearsToPith=|Radius=|DistanceToPith=)"
          comment.text <- cC.text[!grepl(meta.pat, cC.text)]
          comment <- if (length(comment.text) >= 1) paste(comment.text, collapse = "; ") else NA_character_


          ## Get the point coordinates
          # The license info seems to be the 2nd-to-last line before the rest of the coordinates
          # There is a blank line in between
          # This is usually true but not always! Maybe I need more of a positive ID of the
          # coordinates.
          # Find the first line that contains just coordinates. That is, two stings coercible to
          # numeric that are separated by a comma. The key is that there is nothing else in the
          # string.
          # coord.start <- grep("licensedTo", raw.input) + 2

          # This method will not count gaps or season wood - but I'm pretty sure it is not
          # possible to have these as the first point.
          suppressWarnings(
            coord.start <- which(!is.na(sapply(strsplit(raw.input, split = ","),
                                               FUN = \(x) as.numeric(x[[1]]))))[1]
          )

          # get all the coords separated from the header
          ring.bound.raw <- raw.input[coord.start:length(raw.input)]

          # Any ZERO points have to handled first
          suppressWarnings(
            ring.bound.dt <- data.table::data.table(
              ring.bound.raw = ring.bound.raw,
              zeros.rep = 1 + as.numeric(sub(".*?#ZERO", "", ring.bound.raw))
            )
          )
          ring.bound.dt$zeros.rep <- ifelse(is.na(ring.bound.dt$zeros.rep),
                                            1,
                                            ring.bound.dt$zeros.rep)
          ring.bound <- rep(ring.bound.raw, times = ring.bound.dt$zeros.rep)



          ## ring.bound is a character vector, the following turns it into a data.table of coordinates

          ## Vectorized point parse (replaces the per-point mapply + rbindlist)
          x.raw <- unname(ring.bound)
          n <- length(x.raw)

          # D-prefixed points are seasonwood ("W") / gaps; strip the leading D only on those
          is.D           <- substr(x.raw, 1L, 1L) == "D"
          base.str       <- x.raw
          base.str[is.D] <- sub("D", "", x.raw[is.D])          # == noD, but only where needed

          # Coordinate portion = everything before the first " #", trimmed so a trailing
          # space on single points doesn't read as a second pair
          coords.string <- trimws(sub(" #.*", "", base.str))

          # A multipoint now = a non-D coord string with whitespace *between* two pairs
          is.multi <- !is.D & grepl("\\s", coords.string)

          # --- labels & types -------------------------------------------------------
          # Non-D label: text after the first '#', else NA (cleaned to NA later at the "" step)
          labelA <- ifelse(grepl("#", x.raw), sub(".*?#", "", x.raw), NA_character_)

          # D label/type: text after first '#'; a second '#' splits it into type#label
          labelsD   <- sub(".*?#", "", base.str)               # base.str == noD on D rows
          hasSecond <- grepl("#", labelsD)
          typeD     <- ifelse(hasSecond, sub("#.*",   "", labelsD), labelsD)
          labelD    <- ifelse(hasSecond, sub(".*?#",  "", labelsD), NA_character_)
          typeD     <- gsub("%", "", typeD)                    # drop the '%' from gaps
          labelD    <- gsub("^#+|#+$", "", labelD)             # strip stray leading/trailing '#'

          label.el <- ifelse(is.D, labelD, labelA)             # one label per input element

          # --- coordinates ----------------------------------------------------------
          first.pair  <- sub("\\s.*", "", coords.string)                        # "x1,y1"
          second.pair <- ifelse(is.multi, sub("^\\S+\\s+", "", coords.string),  # "x2,y2"
                                NA_character_)

          x1 <- as.numeric(sub(",.*$",    "", first.pair))
          y1 <- as.numeric(sub("^[^,]*,", "", first.pair))
          x2 <- as.numeric(sub(",.*$",    "", second.pair))    # NA where not multi
          y2 <- as.numeric(sub("^[^,]*,", "", second.pair))

          # --- expand (multipoints -> two rows) and assemble ------------------------
          n.rows <- ifelse(is.multi, 2L, 1L)
          oi     <- rep.int(seq_len(n), n.rows)                # output row -> input element
          slot   <- sequence(n.rows)                           # 1 for singles; 1,2 within a multi

          ring.bound.dt <- data.table::data.table(
            series = seriesID,
            x      = data.table::fifelse(slot == 2L, x2[oi], x1[oi]),
            y      = data.table::fifelse(slot == 2L, y2[oi], y1[oi]),
            label  = label.el[oi],
            type   = data.table::fifelse(
              is.D[oi], typeD[oi],
              data.table::fifelse(is.multi[oi],
                                  data.table::fifelse(slot == 1L, "multi1", "multi2"),
                                  "reg"))
          )

          # Files that end in a gap or a W should be trimmed - I think this is rare but is possible
          # Really what we need is that the file ends in a reg point or a multi2 point
          if(!(ring.bound.dt$type[nrow(ring.bound.dt)] %in% c("reg","multi2"))) {
            ring.bound.dt <- ring.bound.dt[1:(nrow(ring.bound.dt) - 1),]
          }

          # Check it twice!
          if(!(ring.bound.dt$type[nrow(ring.bound.dt)] %in% c("reg","multi2"))) {
            ring.bound.dt <- ring.bound.dt[1:(nrow(ring.bound.dt) - 1),]
          }

          # Add the pith coordinates.
          all.coords <- rbind(ring.bound.dt[, c("series","x","y","label","type")], pith.coords.dt)

          # Sometimes the files have "Erroneous order" messages, (points are out of out of order).
          # We can order the points by their x or y position. Leave out pith for now if it's there -
          # it always should go at the end.
          # The following several lines fix the order if it is broken

          # Each set of coords needs a unique identifier
          #all.coords$orig.index <- 1:nrow(all.coords)

          # ggplot(all.coords) +
          #   geom_path(aes(x, y), inherit.aes = F) +
          #   geom_point(aes(x, y, color = type), inherit.aes = F) +
          #   coord_fixed()

          # A fairly simple way to determine if any points are out of order is to see if there are any
          # directional changes in the differences of BOTH the axis. 1 at a time is okay and
          # plausible.
          check.diffs <- all.coords
          check.diffs$x.diff <- c(NA, diff(check.diffs$x))
          check.diffs$y.diff <- c(NA, diff(check.diffs$y))

          # The top row has NA diffs and must be removed
          check.diffs <- check.diffs[-1,]

          # Define the direction of movement
          check.diffs$x.dir <- ifelse(check.diffs$x.diff < 0,
                                      "neg", "pos")
          check.diffs$y.dir <- ifelse(check.diffs$y.diff < 0,
                                      "neg", "pos")

          # 0 diffs should take on the neg or pos label of the preceding point for direction
          # Must also account for 0 diffs at the top of - take the next one instead of the previous
          # one
          for (i in seq_along(check.diffs$x.dir)) {
            if (check.diffs$x.diff[i] == 0) {
              if (length(check.diffs$x.dir[i - 1]) == 0) {
                check.diffs$x.dir[i] <- check.diffs$x.dir[i + 1]
              } else {
                check.diffs$x.dir[i] <- check.diffs$x.dir[i - 1]
              }
            }
          }

          for (i in seq_along(check.diffs$y.dir)) {
            if (check.diffs$y.diff[i] == 0) {
              if (length(check.diffs$x.dir[i - 1]) == 0) {
                check.diffs$y.dir[i] <- check.diffs$y.dir[i + 1]
              } else {
                check.diffs$y.dir[i] <- check.diffs$y.dir[i - 1]
              }
            }
          }

          # Have to exclude the multi2 diffs, as those will represent the diff between multi points.
          # This assumes that the multi points are in order - they will be with respect to each other,
          # because I parsed them.
          check.diffs <- check.diffs[!(check.diffs$type %in% "multi2"),]

          # Find the prevailing direction of each axis
          # This is imperfect in that there can be multiple "prevailing" directions!
          # Tree-ring series may twist in such a way that the short axis can have two prevailing
          # directions. There are a lot of observations in these cases though.
          # I err on the side of producing false negatives rather than missing true negatives.
          # Exclude the pith at this point
          x.prevailing <- aggregate(series ~ x.dir,
                                    data = check.diffs[!(check.diffs$type %in% "pith"),],
                                    length)
          y.prevailing <- aggregate(series ~ y.dir,
                                    data = check.diffs[!(check.diffs$type %in% "pith"),],
                                    length)

          check.diffs$x.head <- ifelse(check.diffs$x.dir %in%
                                         x.prevailing$x.dir[which.max(x.prevailing$series)],
                                       "norm", "div")
          check.diffs$y.head <- ifelse(check.diffs$y.dir %in%
                                         y.prevailing$y.dir[which.max(y.prevailing$series)],
                                       "norm", "div")

          # If we have any observations where BOTH heading components are divergent, then we may have
          # erroneous point order. Read these in anyway, but give a warning and a tag.



          # Set up the null data.table
          error.dt <- data.table::data.table(file = f, message = NA)

          # All points but pith
          if (any(check.diffs$x.head[!(check.diffs$type %in% "pith")] %in% "div" &
                  check.diffs$y.head[!(check.diffs$type %in% "pith")] %in% "div")) {

            warning(paste0("Check coordinates for ", unique(check.diffs$series), ".pos - ",
                           "possible erroneous point order"),
                    call. = FALSE, immediate. = TRUE)

            error.dt <- data.table::data.table(file = f,
                                               message = paste0("Check coordinates - ",
                                                                "possible erroneous point order"))

          }

          # If pith.coords.dt exists,
          # check if pith direction is the same as last set of points,
          # if the distance is really large,
          # or if the pith lies within the range of the long axis.
          if (!is.null(pith.coords.dt)) {
            # Give a separate warning if the pith direction is divergent & based on different
            # criteria: only one heading has to be divergent or if the long axis dist is more
            # than 50% of the range of the other points
            long.axis <- ifelse(which.max(c(
              abs(diff(range(check.diffs$x[!(check.diffs$type %in% "pith")]))),
              abs(diff(range(check.diffs$y[!(check.diffs$type %in% "pith")]))))
            ) == 1,
            "x", "y")

            short.axis <- ifelse(which.min(c(
              abs(diff(range(check.diffs$x[!(check.diffs$type %in% "pith")]))),
              abs(diff(range(check.diffs$y[!(check.diffs$type %in% "pith")]))))
            ) == 1,
            "x", "y")

            # long.axis.range <- range(check.diffs[!(check.diffs$type %in% "pith"),
            #                                      long.axis])
            long.axis.range <- range(check.diffs[!(check.diffs$type %in% "pith"), ..long.axis])
            long.axis.range.diff <- abs(diff(long.axis.range))

            max.diff <- abs(c(check.diffs$x.diff[!(check.diffs$type %in% "pith")],
                              check.diffs$y.diff[!(check.diffs$type %in% "pith")])) |> max()

            short.axis.dir <- paste0(short.axis, ".dir")

            if (
              # pith diffs are larger than max diff for ring widths - ie, ignore the small pith diffs
              (max(abs(c(check.diffs$x.diff[check.diffs$type %in% "pith"],
                         check.diffs$y.diff[check.diffs$type %in% "pith"]))) >
               max.diff)
              &&
              # direction to pith is different than the last two points
              !(check.diffs[check.diffs$type %in% "pith", ..short.axis.dir] %in%
                check.diffs[(nrow(check.diffs) - 1), ..short.axis.dir])
              ||
              # max pith dist is larger than half of the long.axis.range.diff
              max(abs(c(check.diffs$x.diff[check.diffs$type %in% "pith"],
                        check.diffs$y.diff[check.diffs$type %in% "pith"]))) >=
              0.5*long.axis.range.diff
            ) {

              warning(paste0("Check pith location for ", unique(check.diffs$series), ".pos - ",
                             "possible pith location error"),
                      call. = FALSE, immediate. = TRUE)

              if (is.na(error.dt$message)) {
                error.dt$message <- paste0("Check pith location - ",
                                           "possible pith location error")
              } else {
                error.dt$message <- paste(error.dt$message,
                                          paste0("Check pith location - ",
                                                 "possible pith location error"),
                                          sep = "; ")
              }
            }
          }

          # The coordinates can look inverted vertically relative to what I see in CooRecorder.
          # This doesn't matter for the distances

          # Clean up the label
          all.coords$label <- ifelse(all.coords$label %in% "", NA, all.coords$label)

          # type as a factor
          all.coords$type <- factor(all.coords$type,
                                    levels = c("reg", "multi1", "multi2", "W", "gap", "pith"))

          # Calculate the distances in mm using Pythagorean Theorem
          all.coords$x.dist <- c(NA, diff(all.coords$x, lag = 1, differences = 1))
          all.coords$y.dist <- c(NA, diff(all.coords$y, lag = 1, differences = 1))
          all.coords$dist.mm <- sqrt(all.coords$x.dist^2 + all.coords$y.dist^2)

          ## Assigning years based on OD and the order of points
          # If there is season wood, then we need to add the two distances together for each year.
          # If there is a gap, I need to subtract (or skip) the distance indicated by the gap.
          # Gaps can be indicated by 2 points (we skip the distance between them) or by a single point
          # (we skip the distance between a regular point and the gap point - also would apply to
          # "W" points).
          # The order of points is the only thing that indicates the year each point is
          # associated with.
          #all.coords$year1 <- NA
          all.coords$year <- NA
          all.coords$year[all.coords$type %in%
                            c("reg", "multi1")] <- seq(from = OD,
                                                       to = (OD -
                                                               (nrow(all.coords[all.coords$type
                                                                                %in%
                                                                                  c("reg", "multi1"),
                                                               ]) - 1)),
                                                       by = -1)


          for (i in seq_along(all.coords$year)) {
            if (is.na(all.coords$year[i])) {
              all.coords$year[i] <- all.coords$year[i - 1]
            }
          }

          # Year needs to be shifted by 1 to match dist.mm, which are ahead by 1
          all.coords$year <- c(NA, all.coords$year[1:(nrow(all.coords) - 1)])


          all.coords$year[all.coords$type %in% "pith"] <- NA

          # Gaps need some special handling - two gaps in a row need new labels - gap1 & gap2.
          # We will ignore gap2 distances, but keep gap1. Single gaps will stay as "gap", and we will
          # ignore those. This means that any dist AFTER a single gap point will be deleted.
          # I'll need to label these points too.

          # Find where gaps are and if they are sequential - run length encoding
          rle.vec <- rle(as.character(all.coords$type))

          # Initialize an empty result vector
          result <- as.character(all.coords$type)

          # Keep track of the position in the original vector
          pos <- 1

          for (i in seq_along(rle.vec$lengths)) {
            # Check for consecutive "gap" points
            if (rle.vec$values[i] == "gap" && rle.vec$lengths[i] > 1) {
              # Create the relabeling sequence
              labels <- rep(c("gap1", "gap2"), length.out = rle.vec$lengths[i])
              # Assign the new labels
              result[pos:(pos + rle.vec$lengths[i] - 1)] <- labels
            }
            # Also alter the labels of the points following single gaps
            if (rle.vec$values[i] == "gap" && rle.vec$lengths[i] == 1) {
              result[pos + rle.vec$lengths[i]] <- paste0("gap.", result[pos + rle.vec$lengths[i]])
            }

            # Move to the next group
            pos <- pos + rle.vec$lengths[i]
          }

          all.coords$new.type <- result


          # Whole ring width is now the sum of dist.mm within each year, excluding multi2 and gaps
          whole.ring.widths <- all.coords[
            !(new.type %in% c("multi2", "gap2", "pith")) &
              !(substr(new.type, 1L, 4L) %in% "gap.") &
              !is.na(year) & !is.na(dist.mm),          # replicates na.action = na.omit
            .(rw.mm = sum(dist.mm)),                    # == cbind(rw.mm = dist.mm) ... FUN = sum
            keyby = year                               # ascending, matching aggregate's sort
          ]

          # Add the seriesID
          whole.ring.widths$series <- seriesID

          # ggplot(whole.ring.widths, aes(year, rw.mm)) +
          #   geom_line()

          ## Get seasonwood widths if they exist
          # These already exist actually, but need to be labeled properly

          if (any(all.coords$type %in% "W")) {
            #
            seas.wood.widths1 <- all.coords[!(all.coords$new.type %in%
                                                c("multi2", "gap2", "pith")), ]
            # Do some error catching here - check for multiple W points per year. If this happens,
            # skip the EW/LW points and give the user a message about potential errors
            # per.year.check <- aggregate(type ~ year,
            #                             data = seas.wood.widths1,
            #                             FUN = length,
            #                             drop = FALSE)

            per.year.check <- seas.wood.widths1[!is.na(year), .(type = .N), keyby = year]

            if (any(per.year.check$type > 2)) {
              paste0("Possibly multiple seasonwood boundaries detected in ",
                     per.year.check[per.year.check$type > 2, year],". ",
                     "Seasonwood data not read.")
              whole.ring.widths$lw.mm <- NA
              whole.ring.widths$ew.mm <- NA
            } else {

              # Split by year
              seas.wood.widths <- lapply(split(seas.wood.widths1, f = seas.wood.widths1$year),
                                         FUN = \(this.year) {
                                           if (any(this.year$type %in% "W")) {
                                             # have to deal with the gaps if they exist
                                             if (any(this.year$type %in% "gap")) {
                                               for (i in 1:nrow(this.year)) {
                                                 if (this.year$type[i] %in% "gap" &&
                                                     i < nrow(this.year)) {
                                                   # Add "gap" dist to the next row
                                                   this.year$dist.mm[i + 1] <-
                                                     this.year$dist.mm[i + 1] +
                                                     this.year$dist.mm[i]
                                                 }
                                               }
                                             }
                                             this.year <- this.year[!(this.year$type %in% "gap"), ]
                                             this.year$wood.portion <- c("LW", "EW")
                                           } else {
                                             this.year$wood.portion <- "WR"
                                           }
                                           this.year

                                         }) |> data.table::rbindlist()


              # Add the seasonwood points to the whole ring widths
              lw <- seas.wood.widths[seas.wood.widths$wood.portion %in% "LW", c("year", "dist.mm")]
              colnames(lw)[colnames(lw) %in% "dist.mm"] <- "lw.mm"

              ew <- seas.wood.widths[seas.wood.widths$wood.portion %in% "EW", c("year", "dist.mm")]
              colnames(ew)[colnames(ew) %in% "dist.mm"] <- "ew.mm"

              # Merge together and account for years that don't have seas wood with all = TRUE
              whole.ring.widths <- Reduce(f = \(x, y) base::merge(x, y, by = "year", all = TRUE),
                                          list(whole.ring.widths, lw, ew))

            }
          } else {
            whole.ring.widths$lw.mm <- NA
            whole.ring.widths$ew.mm <- NA

          }


          year.labels <- all.coords[
            !is.na(year),
            .(label = {
              lab <- label[!is.na(label)]
              if (length(lab)) paste(as.character(lab), collapse = "; ") else NA_character_
            }),
            keyby = year
          ]

          ## Make the main output data.tables
          whole.ring.widths1 <- whole.ring.widths[year.labels, on = "year", nomatch = NULL]

          # ggplot(whole.ring.widths1, aes(year, rw.mm)) +
          #   geom_line()

          attributes <- data.table::data.table(
            series = seriesID,
            img.DPI = NA,
            d2pith.mm = ifelse(is.null(pith.coords.dt),
                               NA,
                               all.coords$dist.mm[all.coords$type %in% "pith"]),
            out.date = OD,
            in.date = min(whole.ring.widths1$year, na.rm = TRUE),
            total.rw.mm = sum(all.coords$dist.mm[!(all.coords$type %in% "multi2")], na.rm = TRUE)
          )


          error.message <- error.dt$message

          if (!is.na(DPI)) {
            attributes$img.DPI <- DPI
            DPI.message <- NA
            if (DPI < 600) {
              DPI.message <- "DPI is suspiciously low!"
              warning(paste0("DPI is low for ",
                             attributes$series,
                             ".pos. Ring width measurements will ",
                             "be incorrect if DPI is wrong."),
                      call. = FALSE, immediate. = TRUE)
            }
            error.message <- paste(error.dt$message, DPI.message, sep = "; ")
          }

          attributes$radius.mm <- ifelse(is.null(pith.coords.dt),
                                         NA,
                                         attributes$total.rw.mm + attributes$d2pith.mm)
          attributes$comment <- comment

          # Clean up the error message and attach it to attributes
          attributes$error.message <- ifelse(error.message %in% "NA; NA", NA, error.message)

          # Order the attributes by series
          attributes <- attributes[order(attributes$series, decreasing = FALSE),]

          # Clean up and order the ring widths
          rw <- whole.ring.widths1[, c("series", "year", "rw.mm",
                                       "ew.mm", "lw.mm", "label")]

          rw <- rw[order(rw$series, rw$year, decreasing = FALSE),]

          # Assemble the output list
          this.series.list <- list(
            "Ring widths" = rw,
            "Attributes" = attributes,
            "Raw coordinates" = all.coords[, c("series", "year", "x", "y",
                                               "dist.mm", "type", "label")]
          )

          this.series.list

        } else {
          # file doesn't meet criteria above for DENDRO files, CooRecorder version, or is not
          # dated.
          # Get specific here for the individual file
          if (length(grep("DENDRO", raw.input[1])) != 1) {
            tbdr.dt <- data.table::data.table(file = f,
                                              message = "This file could not be identified as a DENDRO file (another data type?)")
          }

          if (CR.ver < 780) {
            tbdr.dt <- data.table::data.table(file = f,
                                              message = enc2utf8(".pos files must be from CooRecorder \u2265 7.8 (update CooRecorder & resave file)"))
          } else {

            if (is.null(OD)) {
              tbdr.dt <- data.table::data.table(file = f,
                                                message = "This file was not dated (no outer year assigned in CooRecorder)")
            }}

          tbdr.dt
        }
      }, error = function(e) {
        # message("Error: Unknown problem with .pos file")
        # Return tbdr.dt in case of an error
        tbdr.dt <- data.table::data.table(file = f,
                                          message = "Unknown problem with .pos file. Check in CooRecorder.")
        tbdr.dt
      }) # end of tryCatch
    } # End pos parsing
    )
    #)
  }

  if (length(out.list) == 1) {
    # If a single series, just return the out.list
    # Make an exception if it falls under "not read"
    if (any(colnames(out.list[[1]]) %in% "message")) {
      warning(paste0("File not read: ", out.list[[1]]$message),
              call. = FALSE, immediate. = TRUE)
    } else {
      out.list[[1]]
    }

  } else { # If many series

    not.read <- vapply(out.list, is.data.frame, logical(1))

    rw    <- lapply(out.list[!not.read], `[[`, "Ring widths")     |> data.table::rbindlist()
    att   <- lapply(out.list[!not.read], `[[`, "Attributes")      |> data.table::rbindlist()
    coord <- lapply(out.list[!not.read], `[[`, "Raw coordinates") |> data.table::rbindlist()
    tbdr  <- out.list[ not.read] |> data.table::rbindlist()

    rw  <- rw[order(series, year)]
    att <- att[order(series)]

    if (nrow(these.are.dups) != 0) {                     # flag AFTER att exists
      dup.series <- gsub(".pos", "", these.are.dups$filenames)
      att[series %in% dup.series,
          error.message := ifelse(is.na(error.message),
                                  "REPLICATED FILE",
                                  paste0("REPLICATED FILE; ", error.message))]
    }

    if (nrow(tbdr) >= 1L) {
      warning("Some files not read. See 'Not read' list for details.",
              call. = FALSE, immediate. = TRUE)
    }

    list("Ring widths" = rw, "Attributes" = att,
         "Raw coordinates" = coord, "Not read" = tbdr)
  }
} ## End of function

