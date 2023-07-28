#' Plot the `ci_detect()` processes
#'
#' @description
#' Takes the nested output list from `ci_detect()` and makes a list of plots that show the process for each series.
#'
#' @param ci_output A list produced by the `ci_detect()` function
#'
#' @details
#' These plots are designed to give the user a window into how the `ci_detect()` process works.
#'
#' @return A nested list of output plots illustrating the disturbance detection & removal iterations (1st list element)
#' for each tree ring series & the final result of the disturbance-free series after all iterations (2nd list element).
#'
#' @import ggplot2
#' @import cowplot
#' @import tidyr
#' @export

plot_ci_detect <- function(ci_output) {

  # Split up the ci_detect output into parts
  or_series <- apply(ci_output[["Original series"]], MARGIN = 2, FUN = \(x) {
    data.frame(year = as.numeric(names(x)), rw = x, type = "Original") |> na.omit()
  }, simplify = FALSE)

  ci_series <- apply(ci_output[["Disturbance-free series"]], MARGIN = 2, FUN = \(x) {
    data.frame(year = as.numeric(names(x)), rw = x, type = "Disturb.-free") |> na.omit()
  }, simplify = FALSE)


  # Map() applies the rbind action across the 2 lists. It is a wrapper for mapply().
  combined_list <- Map(rbind.data.frame, or_series, ci_series)

  # Now to extract the relevant information about the outliers
  # The structure of the outlier iterations list is 1) the iterations, 2a) the entire corrected rwi series for all IDs,
  # 2b) the dataframes of outlier metadata (for just the outlier period) for all IDs or a message that says "No outliers detected".
  # Take the outlier metadata dataframes for each series and add the iteration number, then bind all of them together
  # for each series ID.

  outlier_iter_len <- length(ci_output[["Outlier removal iterations"]])

  whole_series <-
    lapply(ci_output[["Outlier removal iterations"]][2:outlier_iter_len], FUN = \(x) {
      x[[1]]
    })

  # This pulls out just the outlier removal iterations with the specific outlier metadata
  out_curves_only <-
    lapply(ci_output[["Outlier removal iterations"]][2:outlier_iter_len], FUN = \(x) {
      x[[2]]
    })

  out_detection_only <-
    lapply(ci_output[["Outlier removal iterations"]][2:outlier_iter_len], FUN = \(x) {
      x[[3]]
    })
  # The outermost layer of these lists represents the iterations. Within each of those inner lists,
  # there are data.frames or a character vector for each of the series.
  # Add iteration & series variables to the data.frames, then rbind together. Then we use series
  # as a splitting variable later.

  # use a list of the names to pass the iterations to the data

  out_curves_df <- Map(f = \(x, y) {

    # This inner Map() does
    Map(f = \(i, e, y) {
      if (is.character(i)){
        # Skip the "No outliers detected" ones
        i <- data.frame(rwi = NA, year = NA, curve = NA, rwi.cor = NA, dir = NA)
        i$series <- e
        i$iter <- y
        i

      } else {
        i$series <- e
        i$iter <- y
        i
      }
    }, i = x, e = names(x), y = y) |>
      do.call(what = "rbind")

  }, x = out_curves_only, y = names(out_curves_only)) |>
    do.call(what = "rbind")


  # The detection & removal iterations must contain the specific outlier periods
  # (curves & dist-free series) and the entire series as it was at the start of this iteration.
  # This is included in the whole_series list.
  out_detection_df <- Map(
    f = \(x, y, z) {
      Map(
        f = \(i, e, y, z) {
          if (is.character(i)) {
            # Skip the "No outliers detected" ones
            years <- as.numeric(names(z))
            z <-
              data.frame(value = z,
                         type = "Transformed RW")
            z$year <- years
            i <- data.frame(value = NA,
                            type = NA,
                            year = NA)
            i <- rbind(i, z)
            i$series <- e
            i$iter <- y
            i
          } else {
            years <- as.numeric(names(z))
            z <-
              data.frame(value = z,
                         type = "Transformed RW")
            z$year <- years
            i$year <- rep(years, times = 5)
            i <- rbind(i, z)
            i$series <- e
            i$iter <- y
            i
          }
        },
        i = x,
        e = names(x),
        y = y,
        z = z
      ) |>
        do.call(what = "rbind")

    },
    x = out_detection_only,
    y = names(out_detection_only),
    z = whole_series
  ) |>
    do.call(what = "rbind")

  # Convert year var in the curves df to numeric for plotting (already numeric for detection df)
  out_curves_df$year <- out_curves_df$year |> as.numeric()

  # Now split the dfs into lists by series.
  out_curves_split <- split(out_curves_df, out_curves_df$series)
  out_detection_split <- split(out_detection_df, out_detection_df$series)

  ## The outlier detection & removal iteration plots
  out_det_rem_plots <- Map(f = \(det, rem) {
    # Add a single faceting variable to make the process labels
    det$process <- "Detection"
    rem$process <- "Removal"

    # Convert iterations into an ordered factor
    det$iter <- as.numeric(det$iter)
    rem$iter <- as.numeric(rem$iter)

    # Split into lists for the iterations
    det_split <- split(det, det$iter)
    rem_split <- split(rem, rem$iter)

    # Reduce the lists to contain just the iterations with detected disturbances
    # Series with no detected disturbances get a message instead of a plot
    iter_trim <- lapply(rem_split, nrow) > 1
    det_split <- det_split[iter_trim]
    rem_split <- rem_split[iter_trim]

    # If either of these end up with 0 length, then no disturbances were detected
    if (length(det_split) == 0) {
      "No disturbances detected"

    } else {
      Map(f = \(det_iter, rem_iter) {
        det_no_transRW <-
          det_iter[!c(det_iter$type %in% "Transformed RW"),]
        type_fac <- unique(det_no_transRW$type)
        stable_lev <- c("AR residuals", "TBRM", "Detection thresh."," ")
        det_no_transRW$type <- factor(det_no_transRW$type, levels = c(stable_lev, type_fac[!c(type_fac %in% stable_lev)]))

        # set up x-axis breaks for the detection & removal plots:
        dist_year <- min(rem_iter$year, na.rm = TRUE)
        x_breaks <- c(dist_year,
                      labeling::extended(range(det_no_transRW$year)[1],
                                         range(det_no_transRW$year)[2],
                                         m = 5))
        # Remove any base breaks that are too near the disturbance start year
        take_out <- which(abs(dist_year - x_breaks) < 10 &
                            abs(dist_year - x_breaks) > 0)
        if (length(take_out) == 0) {
          clean_breaks <- x_breaks
        } else {
          clean_breaks <- x_breaks[-take_out]
        }

        ## Detection plots
        out_det_plot <-
          ggplot(na.omit(det_no_transRW),
                 aes(year, value, color = type, linetype = type)) +
          scale_color_manual(
            name = NULL,
            values = c("grey80", "black", "black", "black", "orange")
          ) +
          scale_linetype_manual(name = NULL,
                                values = c(1, 1, 3, 3, 1)) +
          annotate(
            geom = "rect",
            xmin = min(rem_iter$year, na.rm = TRUE),
            xmax = max(rem_iter$year, na.rm = TRUE),
            ymin = min(det_iter$value, na.rm = TRUE),
            ymax = max(det_iter$value, na.rm = TRUE),
            color = "grey30"
          ) +
          geom_line() +
          ylab("AR residuals") +
          scale_x_continuous(breaks = clean_breaks) +
          theme(
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            legend.key = element_blank(),
            legend.position = "top"
          ) +
          #coord_fixed(ratio = 45) +
          facet_wrap(~ process, strip.position = "right") + # Use the single factor level to label the plot
          ggtitle(paste0(
            "Core ID: ",
            det_iter$series,
            ", Iteration: ",
            det_iter$iter
          ))


        ## The removal plots
        out_long <-
          pivot_longer(
            rem_iter[, !c(colnames(rem_iter) %in% "rwi.cor")],
            cols = c("curve", "rwi"),
            names_to = "type",
            values_to = "value"
          )
        # add 1 year each before and after to the rwi - this is for plotting aesthetics
        # without this, the outlier/disturbance section appears to float above/below the disturbance-free series.
        cor_series_iter <- det_iter[det_iter$type %in% "Transformed RW",]
        extra_years <- c(min(out_long$year, na.rm = TRUE) - 1,
                         max(out_long$year, na.rm = TRUE) + 1)
        extra_rwi <- cor_series_iter$value[cor_series_iter$year %in% extra_years]

        rem_series <-
          rbind(out_long[, c("value", "type", "year", "series", "iter", "process")],
                data.frame(value = extra_rwi, type = "rwi", year = extra_years,
                           series = out_long$series[1:2], iter = out_long$iter[1:2], process = "Removal"),
                det_iter[det_iter$type %in% "Transformed RW",])
        #rem_series

        # Code the curve fit label & color according to the direction of disturbance
        curve_lab <- ifelse(rem_iter[,"dir"][1] %in% "pos", "Fitted release curve", "Fitted suppression curve")
        curve_col <- ifelse(rem_iter[,"dir"][1] %in% "pos", "blue", "red")

        rem_series$type <-
          ifelse(
            rem_series$type %in% "Transformed RW",
            "Disturb.-free",
            ifelse(rem_series$type %in% "rwi", "Original",
                   curve_lab)
          )

        rem_series$type <-
          factor(rem_series$type,
                 levels = c("Disturb.-free", "Original", curve_lab))

        # process has to be reasserted here to make everything "Removal"
        rem_series$process <-
          "Removal"

        out_rem_plot <-
          ggplot(na.omit(rem_series),
                 aes(year, value, color = type, linewidth = type)) +
          scale_color_manual(name = element_blank(),
                             values = c("black", "grey20", curve_col)) +
          scale_linewidth_manual(name = element_blank(), values = c(0.5, 0.25, 0.5)) +
          scale_linetype_manual(name = element_blank(), values = c(1, 1, 1)) +
          geom_line() +
          ylab("Transformed RW") +
          scale_x_continuous(breaks = clean_breaks) +
          theme(
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.title.x = element_blank(),
            legend.key = element_blank(),
            legend.position = "top"
          ) +
          #coord_fixed(ratio = 45) +
          facet_wrap(~ process, strip.position = "right") # Use a single factor level to label the plot

        # Plot the two panels together
        plot_grid(
          out_det_plot,
          out_rem_plot,
          align = "v",
          axis = "lr",
          ncol = 1
        )
      },
      det_iter = det_split,
      rem_iter = rem_split)
    }
  }, det = out_detection_split, rem = out_curves_split)

  ###
  ## The final output plots (after all iterations)
  # output of this process is a list
  ## use the collapsed iterations in out_detection_split to derive the first year of each disturbance

  out_start_dir <- lapply(out_curves_split, FUN = \(x){
    # Control for the series with no detected disturbances
    # & only make lines for the iterations with detected disturbances
    xNA <- na.omit(x)
    if (nrow(xNA) == 0){
      data.frame(iter = numeric(0), dir = character(0), year = numeric(0))
    } else {
      aggregate(year ~ iter + dir, data = xNA, min)
    }
  })

  final_plots <- Map(f = \(dat, ID, dist) {

    if (nrow(dist) == 0) { # Just give a message if there are no disturbances.
      "No disturbances detected"
    } else {

      dat$type <-
        factor(
          dat$type,
          levels = c(
            "Disturb.-free",
            "Original"
          )
        )

      dist <- dist[order(dist$year),]
      dist$dir <- ifelse(dist$dir %in% "neg", "Suppresions", "Releases")

      # add rw values (for yend) for the geom_segment data
      dist_free <- dat[dat$type %in% "Disturb.-free",]
      dist <- merge(dist, dist_free[,c("year","rw")], by = "year")

      # Make the plot
      ggplot(dat, aes(
        year,
        rw,
        linewidth = type,
      )) +
        geom_segment(data = dist,
                     aes(x = year, xend = year, y = -Inf, yend = rw, color = dir),
                     linewidth = 0.75,
                     inherit.aes = FALSE) +
        scale_color_manual(name = "Disturb. type",
                           values = c("blue","red")) +
        geom_line() +
        scale_linewidth_manual(name = element_blank(), values = c(0.75, 0.25)) +
        ylab("Ring width (mm)") +
        scale_x_continuous(breaks = dist$year) +
        theme(
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.key = element_blank()
        ) +
        ggtitle(paste(ID, "intervention detection results"))
    }
  },
  # To add here: some kind of coord_fixed() with dimensions based on Tufte's rule of thumb of a mean (median?)
  # slope of 45° in the graph. coord_fixed() is not right, we want a relatively fixed ylim, but a flexible xlim

  dat = combined_list, ID = names(combined_list), dist = out_start_dir)

  plot_list <- list(out_det_rem_plots, final_plots)
  names(plot_list) <- c("Disturbance detection & removal plots", "Final disturbance-free series plots")
  plot_list
}
