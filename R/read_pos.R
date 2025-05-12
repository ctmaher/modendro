#' Read in CooRecorder .pos files
#'
#' @description
#' \code{\link{read_pos}} reads in .pos files from CooRecorder (Cybis Elektronik & Data AB, Larsson
#' & Larsson). \code{\link{read_pos}} can handle a single file or a directory containing multiple
#' .pos files. Note that reading in multiple .pos files can be slow if there are 1000s of files.
#'
#' The motivation of this function is to replicate some of the basic operations of CDendro in R.
#' This allows the user to avoid saving tree-ring collections using some of the arcane tree-ring
#' file formats which don't allow users to take full advantage of the great features of CooRecorder,
#' such as point labels, comments, etc. It also calculates the ring widths with greater precision
#' (they are not rounded by default) than what CDendro exports. The outputs of
#' \code{\link{read_pos}} are long-format data.frames, which are much more versatile than the old
#' rwl-format. If you need to export this data as a .rwl, for example, you can use
#' \code{\link{longer_rwl}} and dplR's \code{\link[dplR]{write.rwl}}. Note that for
#' \code{\link[dplR]{write.rwl}} series name lengths are limited.
#'
#' @param path file path to a single .pos file or a directory that may contain several .pos files.
#'
#' @param default.OD optional - a numeric vector of length = 1 to give an outer date (OD) year for
#' files that did not have a date assigned in CooRecorder.
#'
#' @details If path is a directory the function will search all sub-directories for .pos files, thus
#' accommodating a range of directory structures. The main outputs are 1) a "Ring widths" data.frame
#' containing the ring widths (whole ring, late wood, and early wood, as applicable) and any
#' year-specific point labels 2) an "Attributes" data.frame containing the distance to pith, the
#' outer data, the inner-most date, the radius (sum of all ring widths plus the distance to pith),
#' and the comment for whole series, 3) a "Raw coordinates" data.frame containing the orignal
#' coordinates from the .pos file after they have been converted to numeric values, and 4) a
#' "Not read" data.frame which lists the files that were not read in and gives an error message
#' describing a likely reason.
#'
#' The `default.OD` argument is for circumventing a common "Not read" situation in which the user
#' did not assign an OD year to the series in CooRecorder - but be careful!! Sometimes CooRecorder
#' will still give a number for the OD, and this will not be caught by `default.OD`. Always check
#' the outputs by looking at plots (e.g., spag.plot). If you see weird years assigned to the ring
#' widths, then go back to CooRecorder and make sure those series have the correct OD year assigned.
#'
#' \code{\link{read_pos}} contains several error catching heuristics that try to minimize direct
#' interventions from the user or try to clearly guide the user to a fix they can make in
#' CooRecorder.
#'
#' Note: because the output from read_pos includes the `"Raw coordinates"`, you can check for these
#' errors in R. This is much more efficient than going back into to CooRecorder to check multiple
#' .pos files. See the examples below for a workflow using faceted plots in ggplot2.
#'
#' A possible error in CooRecorder is that points can be saved out of order. Usually CooRecorder
#' will give you an "Erroneous point order" message, but it will save the file anyhow. This can
#' wreak havoc on determining ring widths. \code{\link{read_pos}} has a simple way of determining
#' if there are out of order points and it will warn you if it finds something. Since false
#' positives are likely with this approach, \code{\link{read_pos}} will read in the files anyway
#' with warnings the console and in the `message` column of the `"Attributes"` data.frame in the
#' output. This will produce false negatives.
#'
#' Another possible error in CooRecorder is that the pith location gets jumbled around and is no
#' longer valid. This will make any calculations using distance to pith (like age estimates or basal
#' area/basal area increment) wrong! Unfortunately this doesn't cause any error to arise in
#' CooRecorder itself. In attempt to catch these errors, \code{\link{read_pos}} uses a similar logic
#' as for the erroneous points - if the pith veers off in a different direction from the last set
#' of points \emph{or} if the distance to pith from the last point is greater than 50% of the range
#' of the long axis of coordinates, you'll get a warning in the console and in the `message` column
#' of the `"Attributes"` data.frame in the output. This will produce some false negatives.
#'
#' For the forseeable future, \code{\link{read_pos}} will only work with .pos files from CooRecorder
#' 7.8 or greater - which is the earliest version (I think!) to include the actual pith coordinates
#' in the .pos files.
#'
#' In CooRecorder, it is possible to have a seasonwood ("W") point or a gap point at the very end of
#' a measurement series. Since these are nonsensical for determining distances between points, these
#' points are removed by read_pos automatically.
#'
#' \code{\link{read_pos}} doesn't tolerate replicate .pos file names. File names are series names,
#' after all. If you have replicate file names in your path directory, you will get an error message
#' that lists the replicated files. Go fix these in your file system and run \code{\link{read_pos}}
#' again.
#'
#' @return A list containing 4 data.frames for ring widths, attributes, the original coordinates,
#' and a 2-column data.frame of the files not read and the associated error message.
#'
#' @references
#' Larsson & Larsson (2023) \emph{CDendro and CooRecorder programs of the CDendro package},
#'  Cybis Elektronik & Data AB. https://www.cybis.se/forfun/dendro/index.htm
#'
#' @seealso \code{\link{longer_rwl}}, \code{\link{rwl_longer}}
#'
#' @import stats
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
    fdf <- data.frame(filenames = filenames, ind = 1:length(filenames))
    fdf.agg <- aggregate(ind ~ filenames, data = fdf, FUN = \(f) length(f))
    these <- fdf.agg[which(fdf.agg$ind > 1),]

    if (nrow(these) != 0)
      stop(paste("path contains multiple .pos files with the same name:\n",
                 paste(these$filenames, collapse = "\n")))
  }

  out.list <- lapply(pos.files, FUN = \(f) {

    # For internal testing
    #f <- pos.files

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
          pith.coords.df <- data.frame(
            series = seriesID,
            x = as.numeric(pith.coords[[1]][1]),
            y = as.numeric(pith.coords[[1]][2]),
            label = NA,
            type = "pith"
          )

          pith.coords.df$type <- factor(pith.coords.df$type,
                                        levels = c("reg", "multi1", "multi2", "W", "gap", "pith"))
        } else {
          # If no pith coords
          pith.coords.df <- NULL
        }

        ## Get comments if they exist - should be between OD and pith coordinates
        # This needs to be the line after DATED and before Pith Coords or "Written"
        # Written will always exist, but pith coords may not
        # Sometimes these elements are in other places, so we have to control for that too.
        if (length(grep("DATED", raw.input)) >= 1) {
          if (is.null(pith.coords.df)) {
            if ((grep("Written", raw.input) - grep("DATED", raw.input)) > 1) {
              comment.lines <- raw.input[(grep("DATED", raw.input) + 1):(grep("Written",
                                                                              raw.input) - 1)]
              comment <- gsub(pattern = "#C ",
                              replacement = "",
                              comment.lines)
            } else {
              comment <- NA
            }
          } else {
            if ((grep("PithCoordinates", raw.input) - grep("DATED", raw.input)) > 1) {
              comment.lines <- raw.input[(grep("DATED", raw.input) + 1):(grep("PithCoordinates",
                                                                              raw.input) - 1)]
              comment <- gsub(pattern = "#C ",
                              replacement = "",
                              comment.lines)
            } else {
              comment <- NA
            }
          }
        } else { # Rely on the SCALE 1 line instead
          if (is.null(pith.coords.df)) {
            if ((grep("Written", raw.input) - grep("SCALE", raw.input)) > 1) {
              comment.lines <- raw.input[(grep("SCALE", raw.input) + 1):(grep("Written",
                                                                              raw.input) - 1)]
              comment <- gsub(pattern = "#C ",
                              replacement = "",
                              comment.lines)
            } else {
              comment <- NA
            }
          } else {
            if ((grep("PithCoordinates", raw.input) - grep("SCALE", raw.input)) > 1) {
              comment.lines <- raw.input[(grep("SCALE", raw.input) + 1):(grep("PithCoordinates",
                                                                              raw.input) - 1)]
              comment <- gsub(pattern = "#C ",
                              replacement = "",
                              comment.lines)
            } else {
              comment <- NA
            }
          }
        }
        comment <- ifelse(length(comment) > 1, NA, comment)

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
          ring.bound.df <- data.frame(
            ring.bound.raw = ring.bound.raw,
            zeros.rep = 1 + as.numeric(sub(".*?#ZERO", "", ring.bound.raw))
          )
        )
        ring.bound.df$zeros.rep <- ifelse(is.na(ring.bound.df$zeros.rep),
                                          1,
                                          ring.bound.df$zeros.rep)
        ring.bound <- rep(ring.bound.raw, times = ring.bound.df$zeros.rep)


        names(ring.bound) <- 1:length(ring.bound)

        ring.bound.df <- mapply(FUN = \(x, x.names) {
          if (!(substr(x, 1, 1) %in% "D")) {
            # for regular and multi coords
            # check for label, get it if exists
            label_match <- sub(".*?#", "", x)
            label <- ifelse(grepl("#", x), label_match, NA)

            coords.string <- sub(" #.*", "", x)  # Remove everything after ' #'
            # in the cases of multipoints, split by spaces. If single points this leaves it
            # unchanged.
            coords.string1 <- strsplit(coords.string, "\\s+")[[1]]

            if (length(coords.string1) == 1) { # "regular points"
              coords <- strsplit(coords.string1, split = ",")
              coords.df <- data.frame(
                series = seriesID,
                x = as.numeric(coords[[1]][1]),
                y = as.numeric(coords[[1]][2]),
                label = label,
                type = "reg"
              )
            } else { # multi points
              coords1 <- coords.string1[1] |> strsplit(split = ",")
              coords2 <- coords.string1[2] |> strsplit(split = ",")
              coords.df <- data.frame(
                series = seriesID,
                x = c(as.numeric(coords1[[1]][1]), as.numeric(coords2[[1]][1])),
                y = c(as.numeric(coords1[[1]][2]), as.numeric(coords2[[1]][2])),
                label = label,
                type = c("multi1", "multi2")
                #type = c(paste0("multi1_", x.names), paste0("multi2_", x.names))
              )
            }

          } else {
            # For season wood and gaps - these are always single points
            # And remove the "D"
            noD <- sub(pattern = "D", replacement = "", x)

            # check for label, get it if exists.
            # $ indicates go back fomr the end.
            # ^ is go forward from the start
            # Also have to deal with the "#W#" and "#%gap"
            labels <- sub(".*?#", "", noD)
            if (grepl("#", labels)) {
              # Split into primary and secondary labels
              type <- sub("#.*", "", labels)  # Before second #
              label <- sub(".*?#", "", labels)  # After second #
            } else {
              # Single label case - which will be a type (nothing after second #)
              type <- labels
              label <- NA
            }

            # Do some clean up
            # if a gap, remove the "%"
            type <- gsub(pattern = "%", replacement = "", type)
            label <- gsub("^#+|#+$", "", label)

            # Remove everything after ' #'
            coords.string <- sub(pattern = " #.*",
                                 replacement = "",
                                 x = noD)

            coords <- strsplit(coords.string, split = ",")
            coords.df <- data.frame(
              series = seriesID,
              x = as.numeric(coords[[1]][1]),
              y = as.numeric(coords[[1]][2]),
              label = label,
              type = type
            )

          }
          coords.df
        }, x = ring.bound, x.names = names(ring.bound), SIMPLIFY = FALSE) |>
          do.call(what = "rbind")

        # Files that end in a gap or a W should be trimmed - I think this is rare but is possible
        # Really what we need is that the file ends in a reg point or a multi2 point
        if(!(ring.bound.df$type[nrow(ring.bound.df)] %in% c("reg","multi2"))) {
          ring.bound.df <- ring.bound.df[1:(nrow(ring.bound.df) - 1),]
        }

        # Check it twice!
        if(!(ring.bound.df$type[nrow(ring.bound.df)] %in% c("reg","multi2"))) {
          ring.bound.df <- ring.bound.df[1:(nrow(ring.bound.df) - 1),]
        }

        # Add the pith coordinates.
        all.coords <- rbind(ring.bound.df[, c("series","x","y","label","type")], pith.coords.df)

        # Sometimes the files have "Erroneous order" messages, (points are out of out of order).
        # We can order the points by their x or y position. Leave out pith for now if it's there -
        # it always should go at the end.
        # The following several lines fix the order if it is broken

        # Each set of coords needs a unique identifier
        all.coords$orig.index <- 1:nrow(all.coords)

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

        # Give a separate warning if the pith direction is divergent & based on different
        # criteria: only one heading has to be divergent or if the long axis dist is more
        # than 50% of the range of the other points
        long.axis <- ifelse(which.max(c(
          abs(diff(range(check.diffs$x[!(check.diffs$type %in% "pith")]))),
          abs(diff(range(check.diffs$y[!(check.diffs$type %in% "pith")]))))
        ) == 1,
        "x", "y")
        long.axis.range <- abs(diff(range(check.diffs[!(check.diffs$type %in% "pith"),
                                                      long.axis])))

        # Set up the null data.frame
        error.df <- data.frame(file = f, message = NA)

        # All points but pith
        if (any(check.diffs$x.head[!(check.diffs$type %in% "pith")] %in% "div" &
                check.diffs$y.head[!(check.diffs$type %in% "pith")] %in% "div")) {

          warning(paste0("Check coordinates for ", unique(check.diffs$series), ".pos - ",
                         "possible erroneous point order"),
                  call. = FALSE, immediate. = TRUE)

          error.df <- data.frame(file = f,
                                 message = paste0("Check coordinates - ",
                                                  "possible erroneous point order"))

        }

        # Just the pith location - check if pith direction is the same as last set of points
        if (!(check.diffs$x.head[check.diffs$type %in% "pith"] %in%
              check.diffs$x.head[(nrow(check.diffs) - 1)]) ||
            !(check.diffs$y.head[check.diffs$type %in% "pith"] %in%
              check.diffs$y.head[(nrow(check.diffs) - 1)]) ||
            abs(check.diffs[check.diffs$type %in% "pith", paste0(long.axis, ".diff")]) >=
            0.5*long.axis.range) {

          warning(paste0("Check pith location for ", unique(check.diffs$series), ".pos - ",
                         "possible pith location error"),
                  call. = FALSE, immediate. = TRUE)

          if (is.na(error.df$message)) {
            error.df$message <- paste0("Check pith location - ",
                                       "possible pith location error")
          } else {
            error.df$message <- paste(error.df$message,
                                      paste0("Check pith location - ",
                                             "possible pith location error"),
                                      sep = "; ")
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
        # ignore those. This will ignore the distance between the gap point and the point before it
        # (as long as dist.mm starts with NA) - this is incorrect, although this is inconsistent
        # behavior. This means that any dist AFTER a single gap point needs to be deleted.
        # I'll need to label these points too.

        # Find where gaps are and if they are sequential
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
        whole.ring.widths <- aggregate(dist.mm ~ year,
                                       data = all.coords[!(all.coords$new.type %in%
                                                             c("multi2", "gap2", "pith")) &
                                                           !substr(all.coords$new.type, 1, 4) %in%
                                                           "gap.", ],
                                       FUN = sum,
                                       drop = FALSE)

        colnames(whole.ring.widths)[colnames(whole.ring.widths) %in% "dist.mm"] <- "rw.mm"

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
          per.year.check <- aggregate(type ~ year,
                                      data = seas.wood.widths1,
                                      FUN = length,
                                      drop = FALSE)
          if (any(per.year.check[,"type"] > 2)) {
            paste0("Possibly multiple seasonwood boundaries detected in ",
                   per.year.check[per.year.check$type > 2, "year"],".",
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

                                       }) |> do.call(what = "rbind")


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

        # Also need to get the labels back for the appropriate years
        # Prep the labels
        year.labels <- aggregate(
          label ~ year,
          all.coords,
          na.action = na.pass,
          drop = FALSE,
          FUN = \(x) {
            all.labels <- c(x)
            if (any(!is.na(all.labels))) {
              all.labels <- paste(c(as.character(stats::na.omit(
                all.labels
              ))), collapse = "; ")
            } else {
              all.labels <- NA
            }
            all.labels
          }
        )

        ## Make the main output data.frames

        whole.ring.widths1 <- merge(
          whole.ring.widths,
          year.labels,
          by = "year",
          # by.y = "year1",
          all.y = FALSE
        )

        # ggplot(whole.ring.widths1, aes(year, rw.mm)) +
        #   geom_line()

        attributes <- data.frame(
          series = seriesID,
          img.DPI = NA,
          d2pith.mm = ifelse(is.null(pith.coords.df),
                             NA,
                             all.coords$dist.mm[all.coords$type %in% "pith"]),
          out.date = OD,
          in.date = min(whole.ring.widths1$year, na.rm = TRUE),
          total.rw.mm = sum(all.coords$dist.mm[!(all.coords$type %in% "multi2")], na.rm = TRUE)
        )


        error.message <- error.df$message

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
          error.message <- paste(error.df$message, DPI.message, sep = "; ")
        }

        attributes$radius.mm <- ifelse(is.null(pith.coords.df),
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
          tbdr.df <- data.frame(file = f,
                                message = "This file could not be identified as a DENDRO file (another data type?)")
        }

        if (CR.ver < 780) {
          tbdr.df <- data.frame(file = f,
                                message = enc2utf8(".pos files must be from CooRecorder \u2265 7.8 (update CooRecorder & resave file)"))
        } else {

          if (is.null(OD)) {
            tbdr.df <- data.frame(file = f,
                                  message = "This file was not dated (no outer year assigned in CooRecorder)")
          }}

        tbdr.df
      }
    }, error = function(e) {
      # message("Error: Unknown problem with .pos file")
      # Return tbdr.df in case of an error
      tbdr.df <- data.frame(file = f,
                            message = "Unknown problem with .pos file. Check in CooRecorder.")
      tbdr.df
    }) # end of tryCatch
  })

  if (length(out.list) == 1) {
    # If a single series, just return the out.list
    # Make an exception if it falls under "not read"
    if (any(colnames(out.list[[1]]) %in% "message")) {
      warning(paste0("File not read: ", out.list[[1]]$message),
              call. = FALSE, immediate. = TRUE)
    } else {
      out.list[[1]]
    }

  } else {
    # If many series, bind the ring widths and
    # attributes into single data.frames

    rw <- lapply(out.list[which(sapply(out.list, class) %in% "list")], FUN = \(x) {
      x[["Ring widths"]]
    }) |>
      do.call(what = "rbind")

    rw <- rw[order(rw$series, rw$year, decreasing = FALSE),]

    att <- lapply(out.list[which(sapply(out.list, class) %in% "list")], FUN = \(x) {
      x[["Attributes"]]
    }) |>
      do.call(what = "rbind")

    att <- att[order(att$series, decreasing = FALSE),]

    coord <- lapply(out.list[which(sapply(out.list, class) %in% "list")], FUN = \(x) {
      x[["Raw coordinates"]]
    }) |>
      do.call(what = "rbind")

    tbdr <- out.list[which(sapply(out.list, class) %in% "data.frame")] |>
      do.call(what = "rbind")

    new.out.list <- list(rw, att, coord, tbdr)

    names(new.out.list) <- c("Ring widths", "Attributes", "Raw coordinates", "Not read")

    if (length(tbdr) >= 1) {
      warning("Some files not read. See 'Not read' list for details.",
              call. = FALSE, immediate. = TRUE)
    }
    new.out.list
  }
} ## End of function
