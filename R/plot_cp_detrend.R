#' Plot the cp_detrend process
#'
#' @description
#' Takes the output from \code{\link{cp_detrend}} and makes a list of plots that show the Cook &
#' Peters (1997) process for each series.
#'
#'
#' @param cp_out A list produced by the \code{\link{cp_detrend}} function
#'
#' @return A list of output plots showing the power transformation and detrending processes for
#' each series.
#'
#' @references
#' Cook, E. R., and Peters, K. (1997) Calculating unbiased tree-ring indices for the study of
#' climatic and environmental change.
#' \emph{The Holocene}, \strong{7}(3), 361-370.
#'
#' @seealso \code{\link{cp_detrend}}, \code{\link{pwr_t_rwl}}, \code{\link{find_opt_pwr}}
#'
#' @import ggplot2
#' @export
#'
#' @examples
#'
#' library(dplR)
#' data("ca533")
#' ca533_cp <- cp_detrend(ca533, detrend.method = "AgeDepSpline")
#' ca533_cp_plots <- plot_cp_detrend(ca533_cp)
#' # View the first plot
#' ca533_cp_plots[[1]]
#'

plot_cp_detrend <- function(cp_out) {

    # convert all rwl format data.frames to long format, and then rbind them together.
    rw <- as.data.frame(cp_out[["Raw ring widths"]])
    orig.IDs <- colnames(rw)
    rw[, "year"] <-
      rownames(cp_out[["Raw ring widths"]]) |> as.numeric()

    long.rw <- rwl_longer(rwl = rw[,!(colnames(rw) %in% "year")],
                          series.name = "series",
                          dat.name = "value",
                          trim = TRUE,
                          new.val.internal.na = NULL)

    long.rw$type <- "rw"



    trans <- as.data.frame(cp_out[["Transformed ring widths"]][,orig.IDs])
    trans[, "year"] <-
      rownames(cp_out[["Transformed ring widths"]]) |> as.numeric()
    long.trans <- rwl_longer(rwl = trans[,!(colnames(trans) %in% "year")],
                          series.name = "series",
                          dat.name = "value",
                          trim = TRUE,
                          new.val.internal.na = NULL)

    long.trans$type <- "pwr.t_cu"

    curv <- as.data.frame(cp_out[["Detrend curves"]][,orig.IDs])
    curv[, "year"] <-
      rownames(cp_out[["Detrend curves"]]) |> as.numeric()
    long.curv <- rwl_longer(rwl = curv[,!(colnames(curv) %in% "year")],
                             series.name = "series",
                             dat.name = "value",
                             trim = TRUE,
                             new.val.internal.na = NULL)

    long.curv$type <- "pwr.t_cu"

    detr <- as.data.frame(cp_out[["Resid. detrended series"]][,orig.IDs])
    detr[, "year"] <-
      rownames(cp_out[["Resid. detrended series"]]) |> as.numeric()
    long.detr <- rwl_longer(rwl = detr[,!(colnames(detr) %in% "year")],
                            series.name = "series",
                            dat.name = "value",
                            trim = TRUE,
                            new.val.internal.na = NULL)

    long.detr$type <- "de"

    # rbind all the long-format data.frames together
    all.df <- rbind(long.rw, long.trans, long.detr) |>
      na.omit() |> as.data.frame()

    all.df$year <- as.numeric(all.df$year)
    all.df$label <-
      ifelse(
        all.df$type %in% "rw",
        "Ring width (mm)",
        ifelse(
          all.df$type %in% "pwr.t_cu",
          "Transformed RW",
          "Detrended resids."
        )
      )

    # add factor and year to the long.curv data.frame
    long.curv$type <- factor(long.curv$type, levels = c("pwr.t_cu"))
    long.curv$year <- as.numeric(long.curv$year)

    # merge the message data with the long.curv data.
    long.curv <-
      merge(long.curv,
            do.call("rbind", cp_out[["Transformation metadata"]]),
            by = "series")

    long.curv$optimal.pwr <- as.numeric(long.curv$optimal.pwr)

    all.list <- split(all.df, f = as.factor(all.df$series))
    curv.list <- split(long.curv, f = as.factor(long.curv$series))


    Map(
      f = \(x, y, d) {
        # Set up for the plot
        max.trans <-
          x[x[, "type"] %in% "pwr.t_cu", "value"] |> max(na.rm = TRUE)
        this.order <-
          c(unique(x[, "label"])[startsWith(unique(x[, "label"]),
                                            prefix = "Ri")],
            unique(x[, "label"])[startsWith(unique(x[, "label"]),
                                            prefix = "Tr")],
            unique(x[, "label"])[startsWith(unique(x[, "label"]),
                                            prefix = "De")])
        x[, "type"] <- factor(x[, "label"],
                              levels = this.order)
        y[, "type"] <-
          factor(unique(x[, "label"])[startsWith(unique(x[, "label"]),
                                                 prefix = "Tr")],
                 levels = this.order)

        # Set up the messages about transformation & detrending
        z <- y[1, ]
        if (d[1,"method"][[1]] %in% "AgeDepSpline" | d[1,"method"][[1]] %in% "Spline") {
          z$trans.message <- paste0(ifelse(
            z$action %in% "Power transformed",
            paste0(z$action, " (power = ", round(z$optimal.pwr, digits = 3),")"),
            ifelse(z$action %in% "log10 transformed", "log10 transformed", "No transformation")
          ), "\nDetrend method = ", d[1,"method"][[1]]," (nyrs = ", d[1,"nyrs"][[1]],")")
        } else {
          z$trans.message <- paste0(ifelse(
            z$action %in% "Power transformed",
            paste0(z$action, " (power = ", round(z$optimal.pwr, digits = 3),")"),
            ifelse(z$action %in% "log10 transformed", "log10 transformed", "No transformation")
          ), "\nDetrend method = ", d[1,"method"][[1]])
        }

        x_axis_params <- seq(min(na.omit(x[, "year"])),
                             max(na.omit(x[, "year"])),
                             length.out = 6) |>
          round(digits = -1)


        # Have to reference variables in an odd way so that check() doesn't throw a note
        # Strange work-around to "no visible binding for global variable" NOTE
        x_val <- "year"
        y_val <- "value"
        col_val <- "type"
        # Plot
        ggplot2::ggplot(x, aes(.data[[x_val]], .data[[y_val]], color = .data[[col_val]])) +
          ggplot2::scale_color_manual(values = c("black", "black", "blue"),
                             guide = "none") +
          ggplot2::geom_line(linewidth = 0.4, na.rm = TRUE) +
          ggplot2::facet_wrap( ~ type,
                      ncol = 1,
                      scales = "free_y",
                      strip.position = "left") +
          ggplot2::geom_line(data = y,
                    color = "blue",
                    na.rm = TRUE) +
          ggplot2::scale_x_continuous(breaks = x_axis_params,
                             limits = range(x_axis_params)) +
          ggplot2::theme(
            strip.placement = "outside",
            strip.background = element_blank(),
            strip.text = element_text(size = 8),
            strip.clip = "off",
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank()
          ) +
          ggplot2::ggtitle(paste0("C&P transform & detrend for series ID: ", x[1,"series"]),
                  subtitle = z$trans.message)
      },
      x = all.list,
      y = curv.list,
      d = cp_out[["Detrending metadata"]]
    )
} ## End of function

