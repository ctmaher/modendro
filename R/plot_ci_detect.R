#' Plot the ci_detect processes & final results
#'
#' @description
#' Takes the nested output list from \code{\link{ci_detect}}  and makes a list of plots that show
#' the process for each series.
#'
#' @param ci_output A list produced by the \code{\link{ci_detect}} function
#'
#' @details
#' These plots are designed to illustrate how the \code{\link{ci_detect}} process works & visualize
#' the final results for each tree ring series. Correspondingly, there are two kinds of plots (if
#' no disturbances were detected, the function returns a message instead of a plot). The 1st plot
#' type demonstrates the iterative disturbance detection and removal steps. The top panel of this
#' plot shows the detection step for the current disturbance (iteration number is displayed in the
#' plot title). The grey rectangle underlies the time period corresponding to the disturbance, with
#' the raw autoregressive residuals as the grey line, the moving window mean in orange, and the
#' Tukey Biweight Robust Mean and detection thresholds as the horizontal black line and dotted
#' lines, respectively. The bottom panel shows the disturbance removal steps of curve fitting and
#' subtraction on the detrended & transformed ring width series. The blue (releases) or red
#' (suppressions) line segment represents the fitted curve. The thin line segment represents the
#' original series the curve was fitted to. The thicker black line is the resulting
#' "disturbance-free" series after the fitted curve is subtracted. See \code{\link{ci_detect}}
#' for more details on the processes. The shared x-axis for both panels marks evenly placed years
#' and the estimated starting year of the disturbance.
#'
#' The second plot type shows the entire final "disturbance-free" series as a thick black line, the
#' original series as a thin grey line, and releases & disturbances as blue & red vertical lines,
#' respectively. Only the estimated disturbance start years are shown on the x-axis. These results
#' are shown in the original units (ring width), with the long term age/size trend reintroduced &
#' the power transformation reversed.
#'
#'
#' @return A nested list of output plots illustrating the disturbance detection & removal
#' iterations (1st list element) for each tree ring series & the final result of the
#' disturbance-free series after all iterations (2nd list element).
#'
#' @seealso \code{\link{ci_detect}}, \code{\link{dist_det_rem}}, \code{\link{cp_detrend}},
#' \code{\link{plot_cp_detrend}}
#'
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @import labeling
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @examples
#' library(dplR)
#' data("ca533")
#' # Note that this will be somewhat slow, depending on your machine
#'
#' ca533_ci <- ci_detect(rwl = ca533)
#' ca533_ci_plots <- plot_ci_detect(ca533_ci)
#' names(ca533_ci_plots) # See what each list element contains
#' # The first element contains the disturbance detection & removal processes & their iterations
#' # this will display all iterations - scroll backward through your plotting window to see
#' # them all.
#' ca533_ci_plots[[1]][['CAM011']]
#'
#' # The second element contains a plot of the final series compared to the original,
#' # with the disturbances indicated
#' ca533_ci_plots[[2]][['CAM011']]
#'
#' # If you wanted to write these plots to pdf to browse them more freely,
#' # you could do the following steps (will write to your working directory):
#'
#' library(ggplot2)
#'
#' # For the iteration plots:
#' \dontrun{
#' dir.create("ca533_ci_iter_plots") # create a folder first
#' # Then make plot sheets for each series
#' # choose wise values for width & height to make sure your plots aren't squished
#' mapply(FUN = \(p, n) {
#' ggsave(
#' filename = paste0("ca533_ci_iter_plots/", n, "_ci_iter_plots.pdf"),
#' plot = marrangeGrob(p, nrow = length(p), ncol = 1),
#' width = 10, height = length(p)*4
#' )
#' }, p = ca533_ci_plots[[1]][sapply(ca533_ci_plots[[1]], FUN = \(x) !is.character(x))],
#' n = names(ca533_ci_plots[[1]][sapply(ca533_ci_plots[[1]], FUN = \(x) !is.character(x))])
#' )
#'
#' # For the final plots:
#' # choose wise values for width & height to make sure your plots aren't squished
#' ggsave(
#' filename = "ca533_final_ci_plots.pdf",
#' plot = marrangeGrob(ca533_ci_plots[[2]][sapply(ca533_ci_plots[[2]],
#' FUN = \(x) !is.character(x))], nrow=1, ncol=1),
#' width = 10,
#' height = 4
#' )
#' }
#'
#' # You can also use the `modendro` function `plot_cp_detrend()` to make plots of just initial
#' # detrending & transformaton steps.
#' # The required input is the 4th element in the `ci_detect()` output.
#' ca533_cp_plots <- plot_cp_detrend(ca533_ci[[5]])
#' ca533_cp_plots[[1]]


plot_ci_detect <- function(ci_output) {
  # Split up the ci_detect output into parts
  or_series <-
    apply(ci_output[["Cook & Peters detrend"]][["Raw ring widths"]],
          MARGIN = 2,
          FUN = \(x) {
            data.frame(year = as.numeric(names(x)),
                       rw = x,
                       type = "Original") |> na.omit()
          },
          simplify = FALSE)

  orig.IDs <- names(or_series) # capture the original names in their original order

  ci_series <-
    apply(ci_output[["Disturbance-free series"]],
          MARGIN = 2,
          FUN = \(x) {
            data.frame(year = as.numeric(names(x)),
                       rw = x,
                       type = "Disturb.-free") |> na.omit()
          },
          simplify = FALSE)

  ci_series <- ci_series[orig.IDs] # make sure the order matches


  # Map() applies the rbind action across the 2 lists. It is a wrapper for mapply().
  combined_list <- Map(rbind.data.frame, or_series, ci_series)

  # Now to extract the relevant information about the disturbances
  # The structure of the disturbance iterations list is 1) the iterations, 2a) the entire corrected
  # rwi series for all IDs, 2b) the dataframes of disturbance metadata (for just the disturbance
  # period) for all IDs or a message that says "No disturbances detected".
  # Take the disturbance metadata dataframes for each series and add the iteration number,
  # then bind all of them together for each series ID.

  whole_series <-
    lapply(ci_output[["Disturbance removal iterations"]], FUN = \(x) {
      x[["Corrected RWI"]][orig.IDs] # make sure the order matches
    })


  # This pulls out just the disturbance removal iterations with the specific disturbance metadata
  dist_curves_only <-
    lapply(ci_output[["Disturbance removal iterations"]], FUN = \(x) {
      x[["Disturbance curves"]][orig.IDs] # make sure the order matches
    })

  dist_detection_only <-
    lapply(ci_output[["Disturbance removal iterations"]], FUN = \(x) {
      x[["Disturbance detection"]][orig.IDs] # make sure the order matches
    })
  # The outermost layer of these lists represents the iterations. Within each of those inner lists,
  # there are data.frames or a character vector for each of the series.
  # Add iteration & series variables to the data.frames, then rbind together. Then we use series
  # as a splitting variable later.

  # use a list of the names to pass the iterations to the data

  dist_curves_df <- Map(f = \(x, y) {
    # This inner Map() does
    Map(
      f = \(i, e, y) {
        if (is.character(i)) {
          # Skip the "No disturbances detected" ones
          i <-
            data.frame(
              rwi = NA,
              year = NA,
              dir = NA,
              dur = NA,
              x = NA,
              curve = NA,
              eq = NA,
              rwi.cor = NA
            )
          i$series <- e
          i$iter <- y
          i

        } else {
          i$series <- e
          i$iter <- y
          i
        }
      },
      i = x,
      e = names(x),
      y = y
    ) |>
      do.call(what = "rbind")

  },
  x = dist_curves_only,
  y = names(dist_curves_only)) |>
    do.call(what = "rbind")


  # The detection & removal iterations must contain the specific disturbance periods
  # (curves & dist-free series) and the entire series as it was at the start of this iteration.
  # This is included in the whole_series list.
  dist_detection_df <- Map(
    f = \(x, y, z) {
      Map(
        f = \(i, e, y, z) {
          if (is.character(i)) {
            # Skip the "No disturbances detected" ones
            years <- as.numeric(names(z))
            z <-
              data.frame(value = z,
                         type = "Detrended resids.")
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
                         type = "Detrended resids.")
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
    x = dist_detection_only,
    y = names(dist_detection_only),
    z = whole_series
  ) |>
    do.call(what = "rbind")

  # Convert year var in the curves df to numeric for plotting (already numeric for detection df)
  dist_curves_df$year <- dist_curves_df$year |> as.numeric()

  # Now split the dfs into lists by series.
  # make sure the order matches
  dist_curves_split <- split(dist_curves_df, dist_curves_df$series)[orig.IDs]
  # make sure the order matches
  dist_detection_split <- split(dist_detection_df, dist_detection_df$series)[orig.IDs]


  ## The disturbance detection & removal iteration plots
  dist_det_rem_plots <- Map(f = \(det, rem) {
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
          det_iter[!c(det_iter$type %in% "Detrended resids."), ]
        type_fac <- unique(det_no_transRW$type)
        stable_lev <-
          c("AR residuals", "TBRM", "Detection thresh.", " ")
        det_no_transRW$type <-
          factor(det_no_transRW$type,
                 levels = c(stable_lev, type_fac[!c(type_fac %in% stable_lev)]))

        # set up x-axis breaks for the detection & removal plots:
        dist_year <- min(rem_iter$year, na.rm = TRUE)
        x_breaks <- c(dist_year,
                      labeling::extended(
                        range(det_no_transRW$year)[1],
                        range(det_no_transRW$year)[2],
                        m = 5
                      ))
        # Remove any base breaks that are too near the disturbance start year
        take_out <- which(abs(dist_year - x_breaks) < 10 &
                            abs(dist_year - x_breaks) > 0)
        if (length(take_out) == 0) {
          clean_breaks <- x_breaks
        } else {
          clean_breaks <- x_breaks[-take_out]
        }

        ## Detection plots
        # Have to reference variables in an odd way so that check() doesn't throw a note
        x_val <- "year"
        y_val <- "value"
        col_lt_val <- "type"
        dist_det_plot <-
          ggplot2::ggplot(na.omit(det_no_transRW),
                 aes(.data[[x_val]], .data[[y_val]],
                     color = .data[[col_lt_val]], linetype = .data[[col_lt_val]])) +
          ggplot2::scale_color_manual(
            name = NULL,
            values = c("grey80", "black", "black", "black", "orange")
          ) +
          ggplot2::scale_linetype_manual(name = NULL,
                                values = c(1, 1, 3, 3, 1)) +
          ggplot2::annotate(
            geom = "rect",
            xmin = min(rem_iter$year, na.rm = TRUE),
            xmax = (rem_iter$dur - 1) + min(rem_iter$year, na.rm = TRUE),
            ymin = min(det_iter$value, na.rm = TRUE),
            ymax = max(det_iter$value, na.rm = TRUE),
            color = "grey30"
          ) +
          ggplot2::geom_line() +
          ggplot2::ylab("AR resids.") +
          ggplot2::scale_x_continuous(breaks = clean_breaks) +
          ggplot2:: theme(
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            legend.key = element_blank(),
            legend.position = "top"
          ) +
          # Use the single factor level to label the plot
          ggplot2::facet_wrap( ~ process, strip.position = "right") +
          ggplot2::ggtitle(paste0(
            "Series ID: ",
            det_iter$series,
            ", Iteration: ",
            det_iter$iter
          ))


        ## The removal plots
        dist_long <-
          tidyr::pivot_longer(
            rem_iter[,!c(colnames(rem_iter) %in% "rwi.cor")],
            cols = c("curve", "rwi"),
            names_to = "type",
            values_to = "value"
          ) |> as.data.frame()

        # add 1 year each before and after to the rwi - this is for plotting aesthetics
        # without this, the disturbance section appears to float above/below the disturbance-free
        # series. However, don't do this if the disturbance period reaches the beginning or end of
        # the original series
        cor_series_iter <-
          det_iter[det_iter$type %in% "Detrended resids.", ]
        min_dist_year1 <- min(dist_long$year, na.rm = TRUE) - 1
        max_dist_year1 <- max(dist_long$year, na.rm = TRUE) + 1
        min_series_year <- min(cor_series_iter$year, na.rm = TRUE)
        max_series_year <- max(cor_series_iter$year, na.rm = TRUE)
        extra_years <-
          c(
            ifelse(min_dist_year1 < min_series_year, NA, min_dist_year1),
            ifelse(max_dist_year1 > max_series_year, NA, max_dist_year1)
          ) |>
          na.omit() |>
          as.numeric()

        if (length(extra_years) == 0) {
          rem_series <-
            rbind(dist_long[, c("value", "type", "year", "series", "iter", "process")],
                  det_iter[det_iter$type %in% "Detrended resids.", ])
        } else {
          extra_rwi <-
            cor_series_iter$value[cor_series_iter$year %in% extra_years]
          rem_series <-
            rbind(
              dist_long[, c("value", "type", "year", "series", "iter", "process")],
              data.frame(
                value = extra_rwi,
                type = "rwi",
                year = extra_years,
                series = dist_long$series[1],
                iter = dist_long$iter[1],
                process = "Removal"
              ),
              det_iter[det_iter$type %in% "Detrended resids.", ]
            )
        }

        # Code the curve fit label & color according to the direction of disturbance
        curve_lab <-
          ifelse(rem_iter[, "dir"][1] %in% "pos",
                 "Fitted release curve",
                 "Fitted suppression curve")
        curve_col <-
          ifelse(rem_iter[, "dir"][1] %in% "pos", "blue", "red")

        rem_series$type <-
          ifelse(
            rem_series$type %in% "Detrended resids.",
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

        # Have to reference variables in an odd way so that check() doesn't throw a note
        x_val <- "year"
        y_val <- "value"
        col_lt_val <- "type"
        # Plot
        dist_rem_plot <-
          ggplot2::ggplot(na.omit(rem_series),
                 aes(.data[[x_val]], .data[[y_val]],
                     color = .data[[col_lt_val]], linewidth = .data[[col_lt_val]])) +
          ggplot2::scale_color_manual(name = element_blank(),
                             values = c("black", "grey20", curve_col)) +
          ggplot2::scale_linewidth_manual(name = element_blank(), values = c(0.75, 0.25, 0.5)) +
          ggplot2::scale_linetype_manual(name = element_blank(), values = c(1, 1, 1)) +
          ggplot2::geom_hline(yintercept = 0, linetype = 3) +
          ggplot2::geom_line() +
          ggplot2::ylab("Detrended resids.") +
          ggplot2::scale_x_continuous(breaks = clean_breaks) +
          ggplot2::theme(
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.title.x = element_blank(),
            legend.key = element_blank(),
            legend.position = "top"
          ) +
          ggplot2::labs(subtitle = ifelse(
            substr(dist_long$eq[[1]], start = 1, stop = 1) %in% "l",
            dist_long$eq[[1]],
            parse(text = dist_long$eq[[1]])
          )) +
          #coord_fixed(ratio = 45) +
          # Use a single factor level to label the plot
          ggplot2::facet_wrap( ~ process, strip.position = "right")

        # Plot the two panels together
        cowplot::plot_grid(
          dist_det_plot,
          dist_rem_plot,
          align = "v",
          axis = "lr",
          ncol = 1
        )
      },
      det_iter = det_split,
      rem_iter = rem_split)
    }
  }, det = dist_detection_split, rem = dist_curves_split)

  ###
  ## The final output plots (after all iterations)
  # output of this process is a list
  # use the collapsed iterations in dist_detection_split to derive
  # the first year of each disturbance

  dist_start_dir <- lapply(dist_curves_split, FUN = \(x) {
    # Control for the series with no detected disturbances
    # & only make lines for the iterations with detected disturbances
    xNA <- na.omit(x)
    if (nrow(xNA) == 0) {
      data.frame(
        iter = numeric(0),
        dir = character(0),
        year = numeric(0),
        dur = numeric(0)
      )
    } else {
      df <- aggregate(year ~ iter + dir, data = xNA, min)
      df$dur <- nrow(xNA)
      df
    }
  })

  final_plots <- Map(
    f = \(dat, ID, dist) {
      if (nrow(dist) == 0) {
        # Just give a message if there are no disturbances.
        "No disturbances detected"
      } else {
        dat$type <-
          factor(dat$type,
                 levels = c("Disturb.-free",
                            "Original"))

        dist <- dist[order(dist$year), ]
        dist$dir <-
          ifelse(dist$dir %in% "neg", "Suppresions", "Releases")

        # add rw values (for yend) for the geom_segment data
        dist_free <- dat[dat$type %in% "Disturb.-free", ]
        dist <- merge(dist, dist_free[, c("year", "rw")], by = "year")
        # Also need to control the colors for all 3 possible cases
        # (1 release, 1 suppression, or both)
        if (length(unique(dist$dir)) == 2) {
          col_val <- c("blue", "red")
        } else {
          col_val <- ifelse(unique(dist$dir) %in% "Releases", "blue", "red")
        }

        # Make the plot
        # Have to reference variables in an odd way so that check() doesn't throw a note
        x_val <- "year"
        y_val <- "rw"
        col_lt_val <- "type"
        col_var2 <- "dir"
        # Plot
        ggplot2::ggplot(dat, aes(.data[[x_val]],
                        .data[[y_val]],
                        linewidth = .data[[col_lt_val]],)) +
          ggplot2::geom_segment(
            data = dist,
            aes(
              x = .data[[x_val]],
              xend = .data[[x_val]],
              y = -Inf,
              yend = .data[[y_val]],
              color = .data[[col_var2]]
            ),
            linewidth = 0.75,
            inherit.aes = FALSE
          ) +
          ggplot2::scale_color_manual(name = "Disturb. type",
                             values = col_val) +
          ggplot2::geom_line() +
          ggplot2::scale_linewidth_manual(name = "Series", values = c(0.75, 0.25)) +
          ggplot2::ylab("Ring width (mm)") +
          ggplot2::scale_x_continuous(breaks = dist$year) +
          ggplot2::theme(
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "top",
            legend.direction = "horizontal",
            legend.key = element_blank(),
            axis.text.x = element_text(angle = 45)
          ) +
          ggplot2::guides(
            linewidth = guide_legend(title.position = "top", order = 1),
            color = guide_legend(title.position = "top", order = 2)
          ) +
          ggplot2::ggtitle(paste("Intervention detection results for series ID:", ID))
      }
    },

    dat = combined_list,
    ID = names(combined_list),
    dist = dist_start_dir
  )

  plot_list <- list(dist_det_rem_plots[orig.IDs], final_plots[orig.IDs])
  names(plot_list) <-
    c("Disturbance detection & removal plots",
      "Final disturbance-free series plots")
  plot_list
} ## End of function
