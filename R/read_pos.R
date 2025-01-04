#' Read in CooRecorder .pos files
#'
#' @description
#' This function reads in .pos files from CooRecorder (Cybis Elektronik & Data AB, Larsson &
#' Larsson). `read_pos` can handle a single file or a directory containing multiple .pos files.
#'
#' Note that in the current implementation, the multiple .pos files mode can be pretty slow.
#'
#' The motivation of this function is to replicate some of the basic operations of CDendro in R.
#' This allows the user to avoid saving tree-ring collections using some of the arcane tree-ring
#' file formats which don't allow users to take full advantage of the great features of CooRecorder,
#' such as point labels, comments, etc. It also calculates the ring widths with greater precision
#' (they are not rounded by default) than what CDendro exports.
#'
#' For the forseeable future, `read_pos` will only work with .pos files from CooRecorder 7.8 or
#' greater - which is the earliest version (I think!) to include the actual pith coordinates in the
#' .pos files. There is also the stipulation that the file *must* contain the outer year of the ring
#' width measurements. I'm not sure which version of CooRecorder started including this.
#'
#' @param path file path to a single .pos file or a directory that may contain several .pos files.
#'
#' @details If path is a directory the function will search all sub-directories for .pos files, thus
#' accomodating a range of directory structures. The output is slightly different if `path` leads to
#' a single file - the "raw" coordinates are also returned in this case.
#' The main outputs are a "Ring widths" data.frame containing the ring widths (whole ring, late wood
#' , and early wood, as applicable) and any year-specific point labels and an "Attributes"
#' data.frame containing the distance to pith, the outer data, the inner-most date, the radius (sum
#' of all ring widths plus the distance to pith), and the comment for whole series.
#'
#' @return A list containing 2-3 data.frames for ring widths, attributes, and the original
#' coordinates (for single .pos files only). The ring widths are in "long format", which is the more
#' flexible format for most analyses. If you need to export this data as a .rwl, for example, you
#' can use \code{\link{longer_rwl}} and dplR's \code{\link[dplR]{write.rwl}}.
#'
#' @references
#' Larsson & Larsson (2023) \emph{CDendro and CooRecorder programs of the CDendro package},
#'  Cybis Elektronik & Data AB. https://www.cybis.se/forfun/dendro/index.htm
#'
#' @seealso \code{\link{longer_rwl}}, \code{\link{rwl_longer}}, \code{\link{p2yrsL}}
#'
#' @export
#'
#' @examples
#' # example code
#'

read_pos <- function(path = NULL) {
  ## Error checking
  # Is path a character?
  stopifnot(
    "path argument must be a character vector (a file path or paths)" =
      is.character(path) |
      length(path) >= 1
  )

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
      !length((grep("\\.pos{1}", path))) == 0) {
    pos.files <- path
  }

  # There must be at least one .pos file in the path
  stopifnot("path must lead to at least one .pos file" =
              length(pos.files[grep("\\.pos{1}", pos.files)]) >= 1)


  out.list <- lapply(pos.files, FUN = \(f) {

    raw.input <- scan(
      file = f,
      what = "character",
      sep = ";",
      quote = "",
      quiet = TRUE
    )

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

    # File has to be a DENDRO flavored .pos file - top line has this info
    # has to be CooRecorder version 7.8 or higher,
    # and must have a DATED line
    if (length(grep("DENDRO", raw.input[1])) == 1 &&
        CR.ver >= 780 &&
        length(grep("DATED", raw.input)) >= 1) {
      # The crucial data are the seriesID, the outer date, the PithCoordinates,
      # comments if they exist,
      # & the ring boundary coordinates

      ## Get the series ID
      # The top line has the full file path, which includes the file name.
      top.line <- raw.input[1] |> strsplit(split = "\\", fixed = TRUE)
      # The path may be complex, but the file name always ends in .pos
      seriesID <- strsplit(top.line[[1]][length(top.line[[1]])], split = ".pos")[[1]][1]

      ## Get the outer date
      date.line <- raw.input[grep("DATED", raw.input)] |>
        strsplit(split = " ")

      OD.which <- which(!is.na(suppressWarnings(sapply(
        date.line[[1]], as.numeric
      ))))
      OD <- date.line[[1]][OD.which] |> as.numeric()

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
      } else {
        # If no pith coords
        pith.coords.df <- NULL
      }

      ## Get comments if they exist - should be between OD and pith coordinates
      # This needs to be the line after DATED and before Pith Coords or "Written"
      # Written will always exist, but pith coords may not
      # Sometimes these elements are in other places, so we have to control for that too.
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

      comment <- ifelse(length(comment) > 1, NA, comment)

      ## Get the actual coordinates
      # The license info seems to be the 2nd-to-last line before the rest of the coordinates
      # There is a blank line in between
      coord.start <- grep("licensedTo", raw.input) + 2

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


      ring.bound.df <- lapply(ring.bound, FUN = \(x) {
        if (!(substr(x, 1, 1) %in% "D")) {
          # for regular and multi coords
          # check for label, get it if exists
          label_match <- sub(".*?#", "", x)
          label <- ifelse(grepl("#", x), label_match, NA)

          coords.string <- sub(" #.*", "", x)  # Remove everything after ' #'
          # in the cases of multipoints, split by spaces. If single points this leaves it unchanged.
          coords.string1 <- strsplit(coords.string, "\\s+")[[1]]

          if (length(coords.string1) == 1) {
            coords <- strsplit(coords.string1, split = ",")
            coords.df <- data.frame(
              series = seriesID,
              x = as.numeric(coords[[1]][1]),
              y = as.numeric(coords[[1]][2]),
              label = label,
              type = "reg"
            )
          } else {
            coords1 <- coords.string1[1] |> strsplit(split = ",")
            coords2 <- coords.string1[2] |> strsplit(split = ",")
            coords.df <- data.frame(
              series = seriesID,
              x = c(as.numeric(coords1[[1]][1]), as.numeric(coords2[[1]][1])),
              y = c(as.numeric(coords1[[1]][2]), as.numeric(coords2[[1]][2])),
              label = label,
              type = c("multi1", "multi2")
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
      }) |> do.call(what = "rbind")

      all.coords <- rbind(ring.bound.df, pith.coords.df)

      # Clean up the label
      all.coords$label <- ifelse(all.coords$label %in% "", NA, all.coords$label)

      all.coords$type <- factor(all.coords$type,
                                levels = c("reg", "multi1", "multi2", "W", "gap", "pith"))


      # The coordinates look inverted vertically relative to what I see in CooRecorder,
      # but this doesn't matter for the distances

      # Calculate the distances in mm using Pythagorean Theorem
      all.coords$x.dist <- c(NA, diff(all.coords$x))
      all.coords$y.dist <- c(NA, diff(all.coords$y))
      all.coords$dist.mm <- sqrt(all.coords$x.dist^2 + all.coords$y.dist^2)

      # If there is season wood, then we need to add the two distances together for each year.
      # If there is a gap, I need to subtract (or skip) the distance indicated by the gap.
      # Gaps can be indicated by 2 points (we skip the distance between them) or by a single point
      # (we skip the distance between a regular point and the gap point - also would apply to
      # "W" points).
      # The order of points is the only thing that indicates the year each point is associated with.
      all.coords$year1 <- NA
      all.coords$year1[all.coords$type %in%
                         c("reg", "multi1")] <- seq(from = OD,
                                                    to = (OD -
                                                            (nrow(all.coords[all.coords$type
                                                                             %in%
                                                                               c("reg", "multi1"),
                                                            ]) - 1)),
                                                    by = -1)
      # Assigning the correct year to the season wood, gap, and the multi2 points -
      # what I can think of so far is to get the year of the nearest reg point *before* each
      # I don't know how to do that though.
      all.coords$index <- as.numeric(rownames(all.coords))

      # If the series ends in a pith, a multi2, or a gap2, or a combination, we have problems
      year.vec <- mapply(
        FUN = \(x, y) {
          rep(x, each = y)
        },
        x = all.coords$year1[!is.na(all.coords$year1)],
        y = c(diff(all.coords$index[!is.na(all.coords$year1)]), NA),
        SIMPLIFY = FALSE
      ) |>
        do.call(what = "c")

      if (is.null(pith.coords.df) &&
          all.coords$type[nrow(all.coords)] %in% c("reg","multi1")){

        all.coords$year <- year.vec

      } else {
        if (all.coords$type[(nrow(all.coords) - 1)] %in% c("gap2","multi2")) {
          all.coords$year <- c(NA, year.vec, NA)
        } else {
          all.coords$year <- c(NA, year.vec)
        }
      }

      all.coords$year1 <- ifelse(is.na(all.coords$year1),
                                 all.coords$year,
                                 all.coords$year1)

      all.coords$year[all.coords$type %in% "pith"] <-
        all.coords$year1[all.coords$type %in% "pith"] <-
        NA

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
      # Will this work for single gap points? Currently the wrong side of the gap is
      # getting deleted.
      whole.ring.widths <- aggregate(dist.mm ~ year,
                                     data = all.coords[!(all.coords$new.type %in%
                                                           c("multi2", "gap2", "pith")) &
                                                         !substr(all.coords$new.type, 1, 4) %in%
                                                         "gap.", ], sum)

      colnames(whole.ring.widths)[colnames(whole.ring.widths) %in% "dist.mm"] <- "rw.mm"

      # Add the seriesID
      whole.ring.widths$series <- seriesID

      ## Get seasonwood widths if they exist
      # These already exist actually, but need to be labeled properly

      if (any(all.coords$type %in% "W")) {
        #
        seas.wood.widths1 <- all.coords[!(all.coords$new.type %in%
                                            c("multi2", "gap2", "pith")), ]
        #

        #
        seas.wood.widths <- lapply(split(seas.wood.widths1, f = seas.wood.widths1$year),
                                   FUN = \(this.year) {
                                     if (any(this.year$type %in% "W")) {
                                       # have to deal with the gaps if they exist
                                       if (any(this.year$type %in% "gap")) {
                                         for (i in 1:nrow(this.year)) {
                                           if (this.year$type[i] %in% "gap" &&
                                               i < nrow(this.year)) {
                                             # Add "gap" dist to the next row
                                             this.year$dist.mm[i + 1] <- this.year$dist.mm[i + 1] +
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

        whole.ring.widths <- Reduce(f = \(x, y) base::merge(x, y, by = "year"),
                                    list(whole.ring.widths, lw, ew))


      } else {
        whole.ring.widths$lw.mm <- NA
        whole.ring.widths$ew.mm <- NA

      }

      # Also need to get the labels back for the appropriate years
      # Prep the labels
      year.labels <- aggregate(
        label ~ year1,
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
        by.x = "year",
        by.y = "year1",
        all.y = FALSE
      )


      attributes <- data.frame(
        series = seriesID,
        d2pith.mm = ifelse(is.null(pith.coords.df),
                           NA,
                           all.coords$dist.mm[all.coords$type %in% "pith"]),
        out.date = OD,
        in.date = min(whole.ring.widths1$year, na.rm = TRUE),
        total.rw.mm = sum(all.coords$dist.mm[!(all.coords$type %in% "multi2")], na.rm = TRUE)
      )

      attributes$radius.mm <- ifelse(is.null(pith.coords.df),
                                     NA,
                                     attributes$total.rw.mm + attributes$d2pith.mm)
      attributes$comment <- comment


      # also include automatic rwl-format conversion as well?
      # This dosn't make sense because this is just one series.

      this.series.list <- list(
        "Ring widths" = whole.ring.widths1[, c("series", "year", "rw.mm",
                                               "ew.mm", "lw.mm", "label")],
        "Attributes" = attributes,
        "Raw coordinates" = all.coords
      )

      this.series.list
    } else {
      # file doesn't meet criteria above for DENDRO files & CooRecorder version

      warning(
        "Some .pos files in path were not DENDRO files, were not CooRecorder 7.8 or greater, or
  where not dated (no outer year assigned in CooRecorder) - these were skipped"
      )

      this.series.list <- NA
      this.series.list
    }
  })

  if (length(out.list) == 1) {
    # If a single series, just return the out.list
    out.list[[1]]
  } else {
    # If many series, bind the ring widths and
    # attributes into single data.frames

    rw <- lapply(out.list[!is.na(out.list)], FUN = \(x) {
      x[["Ring widths"]]
    }) |>
      do.call(what = "rbind")

    att <- lapply(out.list[!is.na(out.list)], FUN = \(x) {
      x[["Attributes"]]
    }) |>
      do.call(what = "rbind")

    new.out.list <- list(rw, att)

    names(new.out.list) <- c("Ring widths", "Attributes")

    new.out.list
  }
} ## End of function
