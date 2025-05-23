% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_pos.R
\name{read_pos}
\alias{read_pos}
\title{Read in CooRecorder .pos files}
\usage{
read_pos(path = NULL, default.OD = NULL)
}
\arguments{
\item{path}{file path to a single .pos file or a directory that may contain several .pos files.}

\item{default.OD}{optional - a numeric vector of length = 1 to give an outer date (OD) year for
files that did not have a date assigned in CooRecorder.}
}
\value{
A list containing 4 data.frames for ring widths, attributes, the original coordinates,
and a 2-column data.frame of the files not read and the associated error message.
}
\description{
\code{\link{read_pos}} reads in .pos files from CooRecorder (Cybis Elektronik & Data AB, Larsson
& Larsson). \code{\link{read_pos}} can handle a single file or a directory containing multiple
.pos files. Note that reading in multiple .pos files can be slow if there are 1000s of files.

The motivation of this function is to replicate some of the basic operations of CDendro in R.
This allows the user to avoid saving tree-ring collections using some of the arcane tree-ring
file formats which don't allow users to take full advantage of the great features of CooRecorder,
such as point labels, comments, etc. It also calculates the ring widths with greater precision
(they are not rounded by default) than what CDendro exports. Additionally, seasonwood widths are
saved in the same step as whole ring widths, eliminating the extra steps you would have to do in
CDendro.

\code{\link{read_pos}} is also a great tool for basic error checking in .pos files. The raw
coordinates are part of the output, so you can use faceted plots in ggplot2, for example, to
visually check point placement. This is particularly useful for checking pith location estimates.
See details and examples below.

The outputs of
\code{\link{read_pos}} are long-format data.frames, which are much more versatile than the old
rwl-format. If you need to export this data as a .rwl, for example, you can use
\code{\link{longer_rwl}} and dplR's \code{\link[dplR]{write.rwl}}. Note that for
\code{\link[dplR]{write.rwl}} series name lengths are limited.
}
\details{
If path is a directory the function will search all sub-directories for .pos files, thus
accommodating a range of directory structures. The main outputs are 1) a "Ring widths" data.frame
containing the ring widths (whole ring, late wood, and early wood, as applicable) and any
year-specific point labels 2) an "Attributes" data.frame containing the distance to pith, the
outer data, the inner-most date, the radius (sum of all ring widths plus the distance to pith),
and the comment for whole series, 3) a "Raw coordinates" data.frame containing the orignal
coordinates from the .pos file after they have been converted to numeric values, and 4) a
"Not read" data.frame which lists the files that were not read in and gives an error message
describing a likely reason.

The \code{default.OD} argument is for circumventing a common "Not read" situation in which the user
did not assign an OD year to the series in CooRecorder - but be careful!! Sometimes CooRecorder
will still give a number for the OD, and this will not be caught by \code{default.OD}. Always check
the outputs by looking at plots (e.g., spag.plot). If you see weird years assigned to the ring
widths, then go back to CooRecorder and make sure those series have the correct OD year assigned.

\code{\link{read_pos}} contains several error catching heuristics that try to minimize direct
interventions from the user or try to clearly guide the user to a fix they can make in
CooRecorder.

Note: because the output from read_pos includes the \code{"Raw coordinates"}, you can check for these
errors in R. This is much more efficient than going back into to CooRecorder to check multiple
.pos files. See the examples below for a workflow using faceted plots in ggplot2.

A possible error in CooRecorder is that points can be saved out of order. Usually CooRecorder
will give you an "Erroneous point order" message, but it will save the file anyhow. This can
wreak havoc on determining ring widths. \code{\link{read_pos}} has a simple way of determining
if there are out of order points and it will warn you if it finds something. Since false
positives are likely with this approach, \code{\link{read_pos}} will read in the files anyway
with warnings the console and in the \code{message} column of the \code{"Attributes"} data.frame in the
output. This will produce false negatives.

Another possible error in CooRecorder is that the pith location gets jumbled around and is no
longer valid. This will make any calculations using distance to pith (like age estimates or basal
area/basal area increment) wrong! Unfortunately this doesn't cause any error to arise in
CooRecorder itself. In attempt to catch these errors, \code{\link{read_pos}} uses a similar logic
as for the erroneous points - if the pith veers off in a different direction from the last set
of points \emph{or} if the distance to pith from the last point is greater than 50\% of the range
of the long axis of coordinates, you'll get a warning in the console and in the \code{message} column
of the \code{"Attributes"} data.frame in the output. There is a constraint in that this check will
ignore pith distances that within the range of normal ring width distances. This may produce some
false positives, but will keep the false negatives to a manageable level. Assuming that the small
pith distances are close to accurate (i.e., the true pith distance is \emph{not} really large),
then these small errors should have a small affect on the accuracy of estimates using distance
to pith.

For the forseeable future, \code{\link{read_pos}} will only work with .pos files from CooRecorder
7.8 or greater - which is the earliest version (I think!) to include the actual pith coordinates
in the .pos files.

In CooRecorder, it is possible to have a seasonwood ("W") point or a gap point at the very end of
a measurement series. Since these are nonsensical for determining distances between points, these
points are removed by \code{\link{read_pos}} automatically.

\code{\link{read_pos}} doesn't tolerate replicate .pos file names. File names are series names,
after all. If you have replicate file names in your path directory, you will get an error message
that lists the replicated files. Go fix these in your file system and run \code{\link{read_pos}}
again.
}
\examples{

library(ggplot2)
# Read in some example .pos files that show normal files and behavior on files with errors.
ex.pos <- read_pos(system.file("extdata", package = "modendro"))
# We get two erroneous point order warnings - one is real the other is a false positive.
# Below we can see the difference between the two.

names(ex.pos)

# Take a look at the ring widths
ex.pos[["Ring widths"]] |> head()

# Take a look at the attributes
ex.pos[["Attributes"]]

# "Not read" gives you a data.frame of files that were not read in and potentially why
ex.pos[["Not read"]]

# Check the coordinates - this is an efficient way to check point order or pith location errors
ggplot(ex.pos[["Raw coordinates"]], aes(x, y)) +
  geom_path() +
  geom_point(aes(color = type)) +
  facet_wrap(~series, ncol = 1, scales = "free")
# Note that one file truly had erroneous point order - signified by the jagged black line from
# geom_path (which plots points in the order it receives them).
# This file you would want to fix in CooRecorder

# take a look at the ring widths - what you came here for
ggplot(ex.pos[["Ring widths"]], aes(year, rw.mm)) +
  geom_line() +
  facet_wrap(~series, ncol = 1, scales = "free")
# The true erroneous order file has invalid ring widths.
}
\references{
Larsson & Larsson (2023) \emph{CDendro and CooRecorder programs of the CDendro package},
Cybis Elektronik & Data AB. https://www.cybis.se/forfun/dendro/index.htm
}
\seealso{
\code{\link{longer_rwl}}, \code{\link{rwl_longer}}
}
