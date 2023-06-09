#' Plot the `cp_detrend()` process
#'
#' @description
#' Takes the output from `cp_detrend()` and makes a list of plots that show the process for each series.
#'
#'
#' @param cp_out A list produced by the `cp_detrend()` function
#'
#' @return A list of output plots showing the power transformation and detrending processes for each series.
#'
#' @import ggplot2
#' @export

plot_cp_detrend <- function(cp_out) {
  # If we did full detrending
  if (length(cp_out) > 3) {
    # convert all rwl format data.frames to long format, and then rbind them together.
    rw <- as.data.frame(cp_out[["Raw ring widths"]])
    rw[, "year"] <-
      rownames(cp_out[["Raw ring widths"]]) |> as.numeric()
    long.rw <- tidyr::pivot_longer(rw,
                                   cols = -year,
                                   names_to = "series" ,
                                   values_to = "value")
    long.rw$type <- "rw"

    trans <- as.data.frame(cp_out[["Transformed ring widths"]])
    trans[, "year"] <-
      rownames(cp_out[["Transformed ring widths"]]) |> as.numeric()
    long.trans <- tidyr::pivot_longer(
      trans,
      cols = -year,
      names_to = "series" ,
      values_to = "value"
    )
    long.trans$type <- "pwr.t_cu"

    curv <- as.data.frame(cp_out[["Detrend curves"]])
    curv[, "year"] <-
      rownames(cp_out[["Detrend curves"]]) |> as.numeric()
    long.curv <- tidyr::pivot_longer(
      curv,
      cols = -year,
      names_to = "series" ,
      values_to = "value"
    )
    long.curv$type <- "pwr.t_cu"

    detr <- as.data.frame(cp_out[["Resid. detrended series"]])
    detr[, "year"] <-
      rownames(cp_out[["Resid. detrended series"]]) |> as.numeric()
    long.detr <- tidyr::pivot_longer(
      detr,
      cols = -year,
      names_to = "series" ,
      values_to = "value"
    )
    long.detr$type <- "de"

    all.df <- rbind(long.rw, long.trans, long.detr) |>
      na.omit() |> as.data.frame()

    all.df$year <- as.numeric(all.df$year)
    all.df$label <- ifelse(
      all.df$type %in% "rw",
      paste("Raw ring widths - series ID:", all.df$series),
      ifelse(
        all.df$type %in% "pwr.t_cu",
        paste0("Transformed series with fitted trend (", detrend.method, ")"),
        "Residual detrended series"
      )
    )

    # add factor and year to the long.curv data.frame
    long.curv$type <- factor(long.curv$type, levels = c("pwr.t_cu"))
    long.curv$year <- as.numeric(long.curv$year)

    # merge the message data with the long.curv data.
    long.curv <-
      merge(long.curv, cp_out[["Metadata about transformations"]], by = "series")

    all.list <- split(all.df, f = as.factor(all.df$series))
    curv.list <- split(long.curv, f = as.factor(long.curv$series))


    mapply(
      FUN = \(x, y) {
        # Set up for the plot
        max.trans <-
          x[x[, "type"] %in% "pwr.t_cu", "value"] |> max(na.rm = TRUE)
        this.order <-
          c(unique(x[, "label"])[startsWith(unique(x[, "label"]),
                                            prefix = "Ra")],
            unique(x[, "label"])[startsWith(unique(x[, "label"]),
                                            prefix = "Tr")],
            unique(x[, "label"])[startsWith(unique(x[, "label"]),
                                            prefix = "Re")])
        x[, "type"] <- factor(x[, "label"],
                              levels = this.order)
        y[, "type"] <-
          factor(unique(x[, "label"])[startsWith(unique(x[, "label"]),
                                                 prefix = "Tr")],
                 levels = this.order)

        z <- y[1, ]
        z[, "value"] <- max.trans
        z[, "year"] <- median(y[, "year"])
        z$trans.message <- ifelse(
          z$action %in% "Power transformed",
          paste0(z$action, "; power = ", round(z$optimal.pwr, digits = 3)),
          ifelse(z$action %in% "log10 transformed", "No transformation")
        )

        x_axis_params <- seq(min(na.omit(x[, "year"])),
                             max(na.omit(x[, "year"])),
                             length.out = 6) |>
          round(digits = -1)

        ggplot2::ggplot(x, aes(year, value, color = type)) +
          scale_color_manual(values = c("black", "black", "blue"),
                             guide = "none") +
          geom_line(linewidth = 0.25, na.rm = TRUE) +
          facet_wrap( ~ type,
                      ncol = 1,
                      scales = "free_y") +
          geom_line(data = y,
                    color = "blue",
                    na.rm = TRUE) +
          geom_text(data = z,
                    size = 3,
                    aes(label = trans.message)) +
          xlab("Year") +
          scale_x_continuous(breaks = x_axis_params,
                             limits = range(x_axis_params)) +
          theme(
            strip.background = element_blank(),
            strip.text = element_text(hjust = 0),
            strip.clip = "off",
            axis.title.y = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank()
          )
      },
      x = all.list,
      y = curv.list,
      SIMPLIFY = FALSE
    )


  } else {
    # If we did only transformations

    # convert all rwl format data.frames to long format, and then rbind them together.
    rw <- as.data.frame(cp_out[["Raw ring widths"]])
    rw[, "year"] <-
      rownames(cp_out[["Raw ring widths"]]) |> as.numeric()
    long.rw <- tidyr::pivot_longer(rw,
                                   cols = -year,
                                   names_to = "series" ,
                                   values_to = "value")
    long.rw$type <- "rw"

    trans <- as.data.frame(cp_out[["Transformed ring widths"]])
    trans[, "year"] <-
      rownames(cp_out[["Transformed ring widths"]]) |> as.numeric()
    long.trans <- tidyr::pivot_longer(
      trans,
      cols = -year,
      names_to = "series" ,
      values_to = "value"
    )
    long.trans$type <- "pwr.t"


    all.df <- rbind(long.rw, long.trans) |>
      na.omit() |> as.data.frame()

    all.df$year <- as.numeric(all.df$year)
    all.df$label <- ifelse(
      all.df$type %in% "rw",
      paste("Raw ring widths - series ID:", all.df$series),
      "Transformed series"
    )

    # merge the message data with the rest of the data.
    all.df <-
      merge(all.df, cp_out[["Metadata about transformations"]], by = c("series"))

    all.list <- split(all.df, f = as.factor(all.df$series))


    lapply(X = all.list, FUN = \(x) {
      max.trans <-
        x[x[, "type"] %in% "pwr.t", "value"] |> max(na.rm = TRUE)
      this.order <-
        c(unique(x[, "label"])[startsWith(unique(x[, "label"]),
                                          prefix = "Ra")],
          unique(x[, "label"])[startsWith(unique(x[, "label"]),
                                          prefix = "Tr")])
      x[, "type"] <- factor(x[, "label"],
                            levels = this.order)

      z <- x[startsWith(x[, "label"], prefix = "Tr"), ][1, ]
      z[, "value"] <- max.trans
      z[, "year"] <- median(x[, "year"])
      z$trans.message <- ifelse(
        z$action %in% "Power transformed",
        paste0(z$action, "; power = ", round(z$optimal.pwr, digits = 3)),
        ifelse(z$action %in% "log10 transformed", "No transformation")
      )

      x_axis_params <- seq(min(na.omit(x[, "year"])),
                           max(na.omit(x[, "year"])),
                           length.out = 6) |>
        round(digits = -1)

      ggplot2::ggplot(x, aes(year, value)) +
        geom_line(linewidth = 0.25, na.rm = TRUE) +
        facet_wrap( ~ type,
                    ncol = 1,
                    scales = "free_y") +
        geom_text(data = z,
                  size = 2.5,
                  aes(label = trans.message)) +
        xlab("Year") +
        scale_x_continuous(breaks = x_axis_params,
                           limits = range(x_axis_params)) +
        theme(
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0),
          strip.clip = "off",
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank()
        )
    })
  }
} ## End of function
