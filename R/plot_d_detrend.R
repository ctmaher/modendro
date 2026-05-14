#' Plot the d_detrend processes & final results
#'
#' @description
#' Takes the nested output list from \code{\link{d_detrend}} and makes a list of plots that show
#' the disturbance detrending process for each series.
#'
#' @param x A list produced by the \code{\link{d_detrend}} function
#'
#' @details
#' These plots are designed to illustrate how the \code{\link{d_detrend}} process works & visualize
#' the final results for each tree ring series. They are based on the plots in Druckenbrod et al.
#' (2024).
#'
#'
#'
#' @return A nested list of 3 output plots per series illustrating the disturbance detection process
#' (`"Detection plots"`), the curve fitting and removal steps for each event (`"Dist. detrending"`)
#' & the final result of the disturbance-detrended series befre and after final age detrending
#' (`"Result plots"`).
#'
#' @seealso \code{\link{d_detrend}}
#'
#' @import ggplot2
#' @import stats
#' @importFrom cowplot plot_grid
#' @export
#'
#' @examples
#' examples to come

plot_d_detrend <- function(x = NULL) {
  # Control digits for plotting, but restore to what user had before
  current.opciones <- base::options()
  base::on.exit(base::options(current.opciones), add = TRUE)
  base::options(digits = 2)

  ### Error catching
  base::stopifnot("x is not a list output from d_detrend()" =
                    base::class(x) %in% c("list", "d_detrend"))

  ### Basic design is to split the information by series so that we get a set of plots for each series,
  # and separately
  rw.pgc.split <- base::split(x[["PGC"]], f = x[["PGC"]][, "series"])
  events.split <- base::split(x[["Events"]], f = x[["Events"]][, "series"])
  ddtrd.split <- base::split(x[["Dist. detrending"]], f = x[["Dist. detrending"]][, "series"])

  ### Assign variables etc. that don't change series to series
  win.len <- x[["PGC"]]$win.len[1]
  pgc.thresh <- x[["Events"]]$pgc.thresh[1]
  d.detrend.method <- x[["Dist. detrending"]]$d.detrend.method[1]
  detrend.method <- x[["PGC"]]$detrend.method[1]

  # ID the events on the actual ring widths plot.
  pgc.name <- base::paste0("% growth change\n(", win.len, "-year win.)")
  # Get the series names
  series.names <- base::names(rw.pgc.split)

  ##### Run through the series
  series.plots <- base::lapply(series.names, FUN = \(this.series) {
    # separate out the series data.frames once here
    pgc_s  <- rw.pgc.split[[this.series]]
    evt_s  <- events.split[[this.series]]
    ddt_s  <- ddtrd.split[[this.series]]


    ### Make the original series plots with event detection via PGC

    # Helper function to recompute color and scale ranges so that end points are
    # always full red or blue. Ie., this creates an asymmetrical gradient each side of 0.
    make_anchor_scales <- function(data_col, pgc.thresh, pgc.name,
                                   low = "red", mid = "grey70", high = "blue") {
      rng <- base::range(data_col, na.rm = TRUE)
      rng[1] <- base::min(rng[1], -pgc.thresh)
      rng[2] <- base::max(rng[2],  pgc.thresh)
      vals_rescaled <- (base::c(rng[1], 0, rng[2]) - rng[1]) / base::diff(rng)
      base::list(
        color = ggplot2::scale_color_gradientn(name = pgc.name, colors = base::c(low,mid,high),
                                               values = vals_rescaled, limits = rng, na.value = "black"),
        fill = ggplot2::scale_fill_gradientn(name = pgc.name, colors = base::c(low,mid,high),
                                             values = vals_rescaled, limits = rng, na.value = "black")
      )
    }

    ## Show events on the PGC plot

    # Make a data.frame for the threshold hlines
    thresh.df <- base::data.frame(pgc = base::c(-pgc.thresh, pgc.thresh))

    # Also make a new data.frame for the event marking points (triangles), so we can have them pointing
    # in opposite directions. Use an event.type factor for this.
    event.mark.df <- pgc_s[evt_s$max.ind, ]
    event.mark.df$event.type <- base::ifelse(event.mark.df$pgc > 0, "release", "suppression")
    # Also give the event marking points a bit of spacing from the pgc line
    pgc.range <- base::diff(base::range(pgc_s$pgc, na.rm = TRUE))
    event.mark.df$pgc.pos <- base::ifelse(
      event.mark.df$pgc > 0,
      event.mark.df$pgc + pgc.range / 7,
      event.mark.df$pgc - pgc.range / 7
    )

    pgc.plot <- ggplot2::ggplot(pgc_s, ggplot2::aes(year, pgc, color = pgc)) +
      ggplot2::geom_vline(data = event.mark.df,
                          ggplot2::aes(xintercept = year),
                          linetype = 3) +
      ggplot2::geom_line(na.rm = TRUE) +
      ggplot2::geom_point(
        data = event.mark.df,
        ggplot2::aes(year, pgc.pos, fill = pgc, shape = event.type),
        size = 3,
        color = "black"
      ) +
      ggplot2::scale_shape_manual(
        breaks = c("release", "suppression"),
        values = c(25, 24),
        guide = "none"
      ) +
      ggplot2::geom_hline(data = thresh.df, ggplot2::aes(yintercept = pgc, color = pgc)) +
      make_anchor_scales(pgc_s$pgc, pgc.thresh, pgc.name) +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        legend.position = "top",
        legend.position.inside = c(0.5, 0),
        legend.direction = "horizontal",
        axis.title.x = ggplot2::element_blank()
      ) +
      ggplot2::scale_x_continuous(n.breaks = 10) +
      ggplot2::ylab("Percent growth change") +
      ggplot2::ggtitle(this.series, subtitle = paste0(pgc.thresh, "% threshold"))

    # Set up the ends of the rw plot to be black - denoting the lack of coverage of PGC
    get.these.rows1 <- c(1:win.len, (nrow(pgc_s) - win.len):nrow(pgc_s))
    get.these.rows2 <- c(1:(win.len - 1), (nrow(pgc_s) - (win.len - 1)):nrow(pgc_s))
    pgc_s$rw.pgc <- pgc_s$rw
    pgc_s$rw.pgc[-get.these.rows1] <- NA

    rw.pgc.plot <- ggplot2::ggplot() +
      ggplot2::geom_vline(data = event.mark.df,
                          ggplot2::aes(xintercept = year),
                          linetype = 3) +
      ggplot2::geom_line(
        data = pgc_s,
        na.rm = TRUE,
        aes(year, rw.pgc),
        color = "black"
      ) +
      ggplot2::geom_line(data = pgc_s[-get.these.rows2, ],
                         ggplot2::aes(year, rw, color = pgc),
                         na.rm = TRUE) +
      make_anchor_scales(pgc_s$pgc, pgc.thresh, pgc.name) +
      ggplot2::scale_shape_manual(values = c(25, 24), guide = "none") +
      ggplot2::theme(
        panel.background = element_blank(),
        legend.position = "none",
        legend.position.inside = c(0.5, 0),
        legend.direction = "horizontal"
      ) +
      ggplot2::scale_x_continuous(n.breaks = 10) +
      ggplot2::ylab("Ring width (mm)") +
      ggplot2::xlab("Year")

    # Make a grid of these two plots
    detection.plots <- cowplot::plot_grid(
      pgc.plot,
      rw.pgc.plot,
      ncol = 1,
      align = "v",
      axis = "lr"
    )

    ### Make the transformed series plots with iterative event removal

    # there might not be any events
    if (evt_s$message[1] %in% "No disturbance events detected") {
      d.iter.plots <- paste0("No disturbance events detected in series ", this.series)
    } else {
      # Make eventID (the year of the event) a factor for easy faceting
      ddt_s$eventID <- factor(ddt_s$eventID, levels = unique(ddt_s$eventID))

      # Make a labeller to add "release" or "suppression" to the facet label
      ddtrd.split.unique <- base::unique(ddt_s[, c("eventID", "event.type")])
      # Look up "table"
      event_type_lookup <- stats::setNames(
        paste0(
          ddtrd.split.unique$eventID, "\n",
          ddtrd.split.unique$event.type
        ),
        ddtrd.split.unique$eventID
      )

      d.iter.plots <- ggplot2::ggplot(data = ddt_s) +
        ggplot2::geom_line(
          ggplot2::aes(year, pt.rw.ddtrd.i),
          color = "grey40",
          linetype = 1,
          na.rm = TRUE
        ) +
        ggplot2::geom_line(
          ggplot2::aes(year, pt.rw.i),
          color = "black",
          na.rm = TRUE
        ) +
        ggplot2::geom_line(
          ggplot2::aes(year, curve, color = event.type),
          na.rm = TRUE
        ) +
        ggplot2::scale_color_manual(
          name = "Event type",
          breaks = c("release", "suppression"),
          values = c("blue", "red")
        ) +
        ggplot2::facet_wrap(
          ~ eventID,
          ncol = 1,
          strip.position = "right",
          labeller = ggplot2::labeller(eventID = event_type_lookup)
        ) +
        ggplot2::theme(
          panel.background = ggplot2::element_blank(),
          legend.position = "none",
          legend.position.inside = c(0.5, 0),
          legend.direction = "horizontal"
        ) +
        ggplot2::scale_x_continuous(n.breaks = 10) +
        ggplot2::ylab("Transformed ring width") +
        ggplot2::xlab("Year") +
        ggplot2::ggtitle(label = this.series, subtitle = base::bquote({
          D
        }[t] * .(
          base::paste0(" event detrending (", ddt_s$d.detrend.method[1], ")")
        )))
    }

    ### Make the summary ring width plot with detrended series compared to the original, plus show
    # overall detrend curve

    # Make new data.frames in long-format for this one so that we can have proper legends
    if (detrend.method %in% c("none","None")) {

      new.df <- pgc_s[, c("year", "rw", "rw.ddtrd")]

      new.df.long <- base::data.frame(
        year = rep(new.df$year, 2),
        rw   = c(new.df$rw, new.df$rw.ddtrd),
        type = base::rep(c("Original","Dist. detrended"), each = base::nrow(new.df))
      )
      new.df.long$type <- factor(new.df.long$type,
                                 levels = c("Original","Dist. detrended","Age detrend curve"))
    } else {

      new.df <- pgc_s[, c("year", "rw", "rw.ddtrd", "rw.ddtrd.curve")]

      new.df.long <- base::data.frame(
        year = rep(new.df$year, 3),
        rw   = c(new.df$rw, new.df$rw.ddtrd, new.df$rw.ddtrd.curve),
        type = base::rep(c("Original","Dist. detrended","Age detrend curve"), each = base::nrow(new.df))
      )
      new.df.long$type <- factor(new.df.long$type,
                                 levels = c("Original","Dist. detrended","Age detrend curve"))
    }

    # Make dfs for the ribbons too
    ribbon.df <- pgc_s
    ribbon.df$rw <- base::round(ribbon.df$rw, digits = 2)
    ribbon.df$rw.ddtrd <- base::round(ribbon.df$rw.ddtrd, digits = 2)
    ribbon.df$type <- base::ifelse(
      ribbon.df$rw > ribbon.df$rw.ddtrd,
      "release",
      ifelse(ribbon.df$rw < ribbon.df$rw.ddtrd, "suppression", NA)
    )
    ribbon.df.rel <- ribbon.df
    ribbon.df.rel[ribbon.df.rel$type %in% "suppression", c("rw", "rw.ddtrd")] <- NA
    ribbon.df.sup <- ribbon.df
    ribbon.df.sup[ribbon.df.rel$type %in% "release", c("rw", "rw.ddtrd")] <- NA

    # Nicer-looking version of the detrend method
    pretty.d.method <- base::bquote({
      A
    }[t] * .(base::paste0(
      " trend (", pgc_s$detrend.method[1], ")"
    )))

    res.orig.rw.plot <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(
        data = ribbon.df.rel,
        ggplot2::aes(
          year,
          ymin = rw.ddtrd,
          ymax = rw,
          fill = "release",
          alpha = "release"
        ),
        na.rm = TRUE
      ) +
      ggplot2::geom_ribbon(
        data = ribbon.df.sup,
        ggplot2::aes(
          year,
          ymin = rw,
          ymax = rw.ddtrd,
          fill = "suppression",
          alpha = "suppression"
        ),
        na.rm = TRUE
      ) +
      ggplot2::scale_fill_manual(
        name = NULL,
        breaks = c("release", "suppression"),
        values = c("blue", "red"),
        na.value = "black"
      ) +
      ggplot2::scale_alpha_manual(
        name = NULL,
        breaks = c("release", "suppression"),
        values = c(0.6, 0.35, 0)
      ) +
      ggplot2::geom_line(data = new.df.long,
                         ggplot2::aes(year, rw, color = type),
                         na.rm = TRUE) +
      ggplot2::scale_color_manual(
        name = NULL,
        values = c("grey40", "black", "darkorange"),
        labels = c("Original", bquote({
          D
        }[t] * ' detrended'), pretty.d.method)
      ) +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        legend.position = "top",
        legend.position.inside = c(0.5, 0),
        legend.direction = "horizontal",
        axis.title.x = ggplot2::element_blank()
      ) +
      ggplot2::scale_x_continuous(n.breaks = 10) +
      ggplot2::ylab("Ring width (mm)") +
      ggplot2::ggtitle(label = this.series)

    ### Make the final plot showing the results of all the detrending.
    if (pgc_s$detrend.method[1] %in% c("none", "None")) {
      rwi.plot <- ggplot2::ggplot() +
        ggplot2::geom_line(data = pgc_s,
                           ggplot2::aes(year, rw.ddtrd),
                           na.rm = TRUE) +
        ggplot2::theme(
          panel.background = ggplot2::element_blank(),
          legend.position = "bottom",
          legend.position.inside = c(0.5, 0),
          legend.direction = "horizontal"
        ) +
        ggplot2::scale_x_continuous(n.breaks = 10) +
        ggplot2::ylab(bquote({
          D
        }[t] * " detrended ring width (mm)")) +
        ggplot2::xlab("Year")

    } else {

      rwi.plot <- ggplot2::ggplot() +
        ggplot2::geom_line(data = pgc_s,
                           ggplot2::aes(year, ddtrd.index),
                           na.rm = TRUE) +
        ggplot2::theme(
          panel.background = ggplot2::element_blank(),
          legend.position = "bottom",
          legend.position.inside = c(0.5, 0),
          legend.direction = "horizontal"
        ) +
        ggplot2::scale_x_continuous(n.breaks = 10) +
        ggplot2::ylab(bquote({
          A
        }[t] * " & " * {
          D
        }[t] * " detrended RWI")) +
        ggplot2::xlab("Year")
    }

    final.plots <- cowplot::plot_grid(
      res.orig.rw.plot,
      rwi.plot,
      ncol = 1,
      align = "v",
      axis = "lr"
    )

    out.plots <- list(detection.plots, d.iter.plots, final.plots)
    names(out.plots) <- c("Detection plots", "Dist. detrending", "Result plots")
    return(out.plots)
  })

  names(series.plots) <- series.names

  return(series.plots)

} ## End of plot_d_detrend
