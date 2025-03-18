#' Summary plotting function for n_mon_corr output
#'
#' @description
#' Exploratory data analysis (EDA) function to compute correlations between tree ring data and a
#' monthly climate variable aggregated for every combination (lengths 1:12) of consecutive months
#' going back a specified number of years.
#'
#' The results plots show 1) the percentage of tree-ring series that had statistically significant
#' correlations and 2) the mean correlation coefficient (ALL correlations, not just significant
#' ones) with all possible moving window combinations. Moving windows are represented by their start
#' month (the x-axis) and window length (represented by the lines of different colors). Each plot
#' is split by the direction of the relationships: positive (coef > 0) and negative (coef < 0).
#' The x-axis is extended to the left according to the number of lag years the user specifies (with
#' `max.lag`). Lag years are indicated by labeled rectangles just above the x-axis ticks.
#' "0" represents the current year.
#'
#' @param x the n_mon_corr output list object
#'
#'
#' @return 2 plots
#'
#' @import ggplot2
#' @importFrom grDevices hcl.colors
#'
#' @export



plot_n_mon_corr <- function(x = NULL) {
  # Subset out just the significant correlations
  # But there may not always be significant correlations

  # Theoretical - just one place holder for the possible combinations
  sig.only.th <- expand.grid(month = unique(x[["Correlation results"]]$month),
                             win.len = unique(x[["Correlation results"]]$win.len),
                             lag = unique(x[["Correlation results"]]$lag),
                             dir = unique(x[["Correlation results"]]$dir))

  # Empirical - what actually exists (could be none)
  sig.only.em <- x[["Correlation results"]][x[["Correlation results"]]$p <= 0.05, ]

  # Merge them
  sig.only <- merge(sig.only.th, sig.only.em,
                    by = c("month","win.len","lag","dir"),
                    all = TRUE)

  sig.only$lag <- as.character(sig.only$lag)

  res.agg <- aggregate(coef ~ month + win.len + lag + dir,
                       data = sig.only,
                       FUN = \(x) {
                         length(na.omit(x))
                       },
                       drop = FALSE, # make sure the blank combinations still appear
                       na.action = na.pass
  )

  lag.levels <- res.agg$lag |> unique()

  res.agg$lag <- factor(res.agg$lag, levels = lag.levels[order(as.numeric(lag.levels))])

  # Calculate the percentage of significant correlations
  res.agg$prop.sig <- (res.agg$coef / length(unique(x[["Correlation results"]][, "series"]))) *
    100

  # Make a new composite x-axis that combines month and lag
  # Have to control the order of lag and month
  res.agg <- res.agg[order(res.agg$lag, res.agg$month),]
  res.agg$comb.x <- paste(res.agg$lag, res.agg$month, sep = "_")
  res.agg$comb.x <- factor(res.agg$comb.x, levels = unique(res.agg$comb.x))
  res.agg$comb.x.num <- as.numeric(res.agg$comb.x)

  mean.coef.agg <- aggregate(coef ~ month + win.len + lag + dir,
                             data = x[["Correlation results"]],
                             FUN = \(x) mean(x),
                             drop = FALSE) # make sure the blank combinations still appear
  colnames(mean.coef.agg)[colnames(mean.coef.agg) %in% "coef"] <- "mean.coef"

  res.agg <- merge(res.agg, mean.coef.agg, by = c("month", "win.len", "lag", "dir"))
  #res.agg$lag <- factor(res.agg$lag, levels = lag.levels[order(as.numeric(lag.levels))])
  res.agg$dir <- factor(res.agg$dir, levels = c("Pos.", "Neg."))
  res.agg <- res.agg[order(res.agg$lag, res.agg$month),]

  # Make a plot.
  # These 4 lines are to deal with "no visible binding" NOTEs from check()
  x_var <- "comb.x.num"
  x_lab <- "month"
  y_var <- "prop.sig"
  y_var2 <- "mean.coef"
  col_var <- "win.len"
  x.intercept <- "xint"
  y.intercept <- "yint"

  # For geom_tile & geom_text
  lag.lab.df1 <- res.agg[, c("lag", "month", "comb.x", "comb.x.num")] |> unique()
  lag.lab.df1$dir <- factor("Neg.", levels = c("Pos.","Neg."))

  lag.lab.df <- aggregate(comb.x.num ~ lag + dir, data = lag.lab.df1, min)
  colnames(lag.lab.df)[colnames(lag.lab.df) %in% "comb.x.num"] <- "x.min"
  lag.lab.df$x.max <- aggregate(comb.x.num ~ lag + dir, data = lag.lab.df1, max)[,3] + 0.85
  lag.lab.df$y.min <- -10; lag.lab.df$y.max <- -1
  # Define some ylim values here
  ylim.val <- max(abs(na.omit(res.agg$mean.coef))) + 0.1
  lag.lab.df$y.min2 <- -ylim.val; lag.lab.df$y.max2 <- -(ylim.val + 0.15*ylim.val)
  lag.lab.df$lag <- factor(lag.lab.df$lag, levels = lag.levels[order(lag.levels)])

  x_min <- "x.min"
  x_max <- "x.max"
  y_min <- "y.min"
  y_max <- "y.max"
  y_min2 <- "y.min2"
  y_max2 <- "y.max2"
  lag_lab <- "lag"

  # Derive some of the arguments from the output
  clim.var <- unique(x[["Correlation results"]][, "clim.var"])
  agg.fun <- unique(x[["Climate data (raw)"]][, "agg.fun"])
  prewhiten <- !is.null(x[["Climate data (prewhitened)"]])
  max.win <- max(as.numeric(x[["Correlation results"]][, "win.len"]))
  win.align <- unique(x[["Climate data (raw)"]][, "win.align"])
  corr.method <- unique(x[["Correlation results"]][, "corr.method"])
  hemisphere <- unique(x[["Correlation results"]][, "hemisphere"])
  pretty.corr.method <- ifelse(corr.method %in% "pearson",
                               "Pearson ",
                               ifelse(corr.method %in% "spearman",
                                      "Spearman ",
                                      "Kendall "))

  per.plot <-
    ggplot2::ggplot(res.agg, ggplot2::aes(.data[[x_var]], .data[[y_var]],
                                          color = as.factor(.data[[col_var]]))) +
    ggplot2::ggtitle(paste0(pretty.corr.method,
                            "correlations between tree-ring and ",
                            clim.var,
                            " series"),
                     subtitle = paste0("Transformation: ",
                                       ifelse(prewhiten == TRUE, "ARIMA resids.", "None"),
                                       "; Aggregation of climate moving windows: ",
                                       paste0(agg.fun, "s"))) +
    ggplot2::scale_color_manual("Moving\nwindow\nlength\n(n months)",
                                values = grDevices::hcl.colors(max.win, palette = "Spectral")) +
    ggplot2::geom_rect(data = lag.lab.df,
                       ggplot2::aes(xmin = .data[[x_min]],
                                    xmax = .data[[x_max]],
                                    ymin = .data[[y_min]],
                                    ymax = .data[[y_max]]),
                       inherit.aes = FALSE,
                       fill = "grey15") +
    ggplot2::geom_text(data = lag.lab.df,
                       ggplot2::aes(x = (.data[[x_min]] + .data[[x_max]]) / 2,
                                    y = (.data[[y_min]] + .data[[y_max]]) / 2,
                                    label = .data[[lag_lab]]),
                       inherit.aes = FALSE,
                       color = "white",
                       size = 3) +
    ggplot2::geom_line(linewidth = 0.75,
                       ggplot2::aes(group = factor(.data[[col_var]], levels = rev(1:max.win)))) +
    ggplot2::facet_wrap( ~ dir, ncol = 1, strip.position = "right") +
    ggplot2::scale_x_continuous(breaks = unique(res.agg$comb.x.num),
                                labels = rep(1:12,
                                             times = length(unique(res.agg$comb.x.num))/12)) +
    ggplot2::ylab(
      paste(
        "Percentage of",
        length(unique(x[["Correlation results"]][, "series"])),
        "total series\nrecording significant correlations"
      )
    ) +
    ggplot2::xlab(ifelse(win.align %in% "right", "End month", "Start month")) +
    ggplot2::coord_cartesian(ylim = c(-5, 100),
                             xlim = c(min(res.agg$comb.x.num), max(res.agg$comb.x.num))) +
    ggplot2::theme_dark() +
    ggplot2::theme(
      panel.spacing.x = ggplot2::unit(-0.1, "lines"),
      panel.background = ggplot2::element_rect(fill = "black"),
      plot.background = ggplot2::element_rect(fill = "black"),
      legend.background = ggplot2::element_rect(fill = "black"),
      panel.grid = ggplot2::element_line(color = "grey40"),
      legend.text = ggplot2::element_text(color = "white"),
      legend.title = ggplot2::element_text(color = "white"),
      axis.title = ggplot2::element_text(color = "white"),
      axis.text = ggplot2::element_text(color = "white"),
      plot.title = ggplot2::element_text(color = "white"),
      plot.subtitle = ggplot2::element_text(color = "white"),
      legend.position = "right"
    )

  ### Mean coef plot

  mean.plot <-
    ggplot2::ggplot(res.agg, ggplot2::aes(.data[[x_var]], .data[[y_var2]],
                                          color = as.factor(.data[[col_var]]))) +
    ggplot2::ggtitle(paste0(pretty.corr.method,
                            "correlations between tree-ring and ",
                            clim.var,
                            " series"),
                     subtitle = paste0("Transformation: ",
                                       ifelse(prewhiten == TRUE, "ARIMA resids.", "None"),
                                       "; Aggregation of climate moving windows: ",
                                       paste0(agg.fun,"s"))) +
    ggplot2::scale_color_manual("Moving\nwindow\nlength\n(n months)",
                                values = grDevices::hcl.colors(max.win, palette = "Spectral")) +
    # Invisible hline to set the ylim
    ggplot2::geom_hline(
      data = data.frame(
        yint = c(ylim.val,0,0,-ylim.val), # This defines the ylim
        dir = factor(c("Pos.","Pos.","Neg.","Neg."), levels = c("Neg.","Pos."))
      ),
      ggplot2::aes(yintercept = .data[[y.intercept]]),
      color = NA
    ) +
    ggplot2::geom_rect(data = lag.lab.df,
                       ggplot2::aes(xmin = .data[[x_min]],
                                    xmax = .data[[x_max]],
                                    ymin = .data[[y_min2]],
                                    ymax = .data[[y_max2]]),
                       inherit.aes = FALSE,
                       fill = "grey15") +
    ggplot2::geom_text(data = lag.lab.df,
                       ggplot2::aes(x = (.data[[x_min]] + .data[[x_max]]) / 2,
                                    y = (.data[[y_min2]] + .data[[y_max2]]) / 2,
                                    label = .data[[lag_lab]]),
                       inherit.aes = FALSE,
                       color = "white",
                       size = 3) +
    ggplot2::geom_line(linewidth = 0.75,
                       ggplot2::aes(group = factor(.data[[col_var]], levels = rev(1:12)))) +
    ggplot2::facet_wrap( ~ dir, ncol = 1, strip.position = "right", scales = "free_y") +
    # ggplot2::geom_vline(
    #   data = data.frame(
    #     xint = ifelse(hemisphere == "S",
    #                   res.agg$comb.x.num[res.agg$comb.x %in% paste0("+1_", gro.period.end)],
    #                   res.agg$comb.x.num[res.agg$comb.x %in% paste0("0_",gro.period.end)])
    #   ),
    #   ggplot2::aes(xintercept = .data[[x.intercept]]),
    #   color = "white"
    # ) +
    ggplot2::scale_x_continuous(breaks = unique(res.agg$comb.x.num),
                                labels = rep(1:12,
                                             times = length(unique(res.agg$comb.x.num))/12)) +
    ggplot2::scale_y_continuous(breaks = seq(-1, 1, by = 0.2)) +
    ggplot2::ylab(
      paste(
        "Mean correlation coefficent of\n",
        length(unique(x[["Correlation results"]][, "series"])),
        "total series"
      )
    ) +
    ggplot2::xlab(ifelse(win.align %in% "right", "End month", "Start month")) +
    ggplot2::coord_cartesian(xlim = c(min(res.agg$comb.x.num), max(res.agg$comb.x.num))) +
    ggplot2::theme_dark() +
    ggplot2::theme(
      panel.spacing.x = ggplot2::unit(-0.1, "lines"),
      panel.background = ggplot2::element_rect(fill = "black"),
      plot.background = ggplot2::element_rect(fill = "black"),
      legend.background = ggplot2::element_rect(fill = "black"),
      panel.grid = ggplot2::element_line(color = "grey40"),
      legend.text = ggplot2::element_text(color = "white"),
      legend.title = ggplot2::element_text(color = "white"),
      axis.title = ggplot2::element_text(color = "white"),
      axis.text = ggplot2::element_text(color = "white"),
      plot.title = ggplot2::element_text(color = "white"),
      plot.subtitle = ggplot2::element_text(color = "white"),
      legend.position = "right"
    )

  return(c(per.plot, mean.plot))

} # End of function
