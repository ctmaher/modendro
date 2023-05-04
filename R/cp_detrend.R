#' Detrend raw ring widths using Cook & Peters' (1997) method.
#'
#' @description
#' Compute variance-stabilized & residual detrended ring width indices that minimize the series end effect problems
#' that can result from using ratios to derive detrended indices.
#'
#'
#' @details
#' For decades, the default method of removing long-term size/age trends from tree ring width series has been to fit a curve
#' (e.g., a negative exponential equation) to the data, then divide the raw data by the predicted (trend) values.
#' The resulting indices have ± stable variance for their duration and a mean of 1, and are ready to be aggregated
#' together into chronologies with other series standardized via the same method. However, this method can introduce artificial
#' trends in locations along the series where the trend line has a poor fit with the data, especially when the predicted values drop below 0.5.
#' This is common at the end of series. The risk is that the researcher might obtain biased estimates of growth trends - i.e., enhanced growth
#' increases or decreases.
#'
#' Cook and Peters (1997) proposed the method employed here as a fix this problem.
#' The basic steps are to estimate the optimal power of transformation via the local spread versus level relationship of each series,
#' where local spread is defined as the absolute value of the first
#' differences, S, (|rwt - rwt-1|) and the local level is the arithmetic mean of each pair of adjacent values, M, (rwt + rwt-1)/2.
#' The spread versus level relationship is then modeled in a simple linear regression as log10(S) ~ log10(M). The optimal power of transformation
#' is then estimated as p = |1-slope| of this regression model. If p ≤ 0.1 (i.e., slope is near 1), a log10() transformation is applied.
#' If p > 1 (i.e., slope is near 0), no transformation is applied. Otherwise, p is applied to all the raw ring widths as rw^p.
#'
#' With no special accommodation, this power transformation process can produce odd results for 0 values (missing rings).
#' Log transformation will obviously fail in these circumstances. To reduce these problems, 0 ring width values are replaced
#' with the minimum value possible given the resolution of the data, which for tree ring data is typically 0.01 or 0.001 mm.
#' If you have many 0 ring width values in your data, make sure to take a look at the output plots (you should do this anyway) to
#' check for any weirdness.
#'
#'
#' @param rwl A rwl object (read in by dplR's `read.rwl()`)
#' @param detrend.method Character string of the detrending method to use. Passes to dplR's `detrend()` function.
#' @param nyrs Numeric vector, used in dplR's `detrend()` function for the "Spline" and "AgeDepSpline" methods.
#'
#' @return A list with the following elements:
#' 1) the residual detrended power transformed series.
#' 2) the curves fitted to the transformed series.
#' 3) the power transformed series (before detrending).
#' 4) messages for each series about the type of transformation applied.
#' 5) a list of output plots showing the power transformation and detrending processes for each series.
#' If no detrending is performed, then only elements 3-5 are in the output list.
#'
#' @references
#' Cook, E. R., and Peters, K. (1997) Calculating unbiased tree-ring indices for the study of climatic and environmental change.
#' \emph{The Holocene}, \strong{7}(3), 361-370.
#'
#' @import dplR
#' @import ggplot2
#' @import tidyr
#' @export
#'
#' @examples
#' library(dplR)
#' library(ggplot2)
#' # Bristlecone pine tree ring collection from Campito mountain, White Mountains, CA
#' # Many 0 value rings in this collection
#' data("ca533")
#' ca533_cp <- cp_detrend(rwl = ca533,
#'                        detrend.method = "AgeDepSpline",
#'                        nyrs = 50)
#' # Use indexing to look at plots
#' ca533_cp[[5]][[1]] # The plot of the first series
#'
#' # If you want, you can save the plots (they are ggplot objects)
#' # to disk with something like this:
#'# Create a new directory for the plots (there could be a lot!)
#' cp_plot_dir <- paste0(getwd(),"/CP_plots")
#' dir.create(cp_plot_dir)
#' # Apply `ggsave()` to each plot:
#' mapply(FUN = \(x, y) {
#'         ggsave(filename = paste0(cp_plot_dir, "/", y, ".pdf"),
#'         plot = x,
#'         device = "pdf",
#'         width = 20,
#'         height = 14,
#'         units = "cm")
#' }, x = ca533_cp[[5]], y = names(ca533_cp[[5]]))

cp_detrend <-
  function(rwl,
           detrend.method = "none",
           nyrs = NULL) {
    # Power transform series prior to detrending using methods of Cook and Peters (1997)

    # Error catching
    stopifnot("rwl is not an object of class 'rwl', 'data.frame', or 'matrix'" =
                data.class(rwl) %in% "rwl" |
                data.class(rwl) %in% "data.frame" |
                data.class(rwl) %in% "matrix"
    )

    stopifnot("rwl has no rownames (must be years only) or no colnames (must be series IDs only)" =
                !is.null(rownames(rwl)) |
                !is.null(colnames(rwl))
    )

    if (apply(rwl, MARGIN = 2, FUN = \(x) all(is.na(x))) |> any() == TRUE) {
      these_are_NA <- colnames(rwl)[which(apply(rwl, MARGIN = 2, FUN = \(x) all(is.na(x))) == TRUE)]
      stop("The following series have no values (all NAs): " , paste(these_are_NA, collapse = ", "))
    }

    stopifnot("detrend.method is not a known option. See ?dplR::detrend.series for options" =
                detrend.method %in% c("none","None","Spline", "ModNegExp", "Mean",
                                      "Ar", "Friedman", "ModHugershoff","AgeDepSpline"))

    if (detrend.method %in% c("none", "None")) {
      message("Proceding with the default of no detrending - series will be transformed only.
            Make sure you want this, otherwise chose a detrending method wisely")
    }

    if (any(rwl < 0, na.rm = TRUE)) {
      rwl[which(rwl < 0),] <- NA
      warning("One or more negative values detected in the rwl file. These were replaced with NAs.
              Check your data to make sure this is correct (e.g., sometimes -9999 is used as a place-holder value for 0 rings.")
    }

    # Find the optimal power of transformation and transform the series
    trans.list <- pwr_t_rwl(rwl) # Returns a list
    trans <- trans.list[[1]]
    messages <- trans.list[[2]]

    # detrend the transformed series -or- don't
    if (!c(is.null(detrend.method) | detrend.method %in% "None" | detrend.method %in% "none")) {

      detr.result <- detrend(trans, method = detrend.method, make.plot = FALSE,
                             difference = TRUE, nyrs = nyrs, return.info = TRUE)

      curv <- detr.result$curves
      detr <- detr.result$series


      # Make the plots

      # convert all rwl format data.frames to long format, and then rbind them together.
      rw <- as.data.frame(rwl)
      rw[, "year"] <- rownames(rwl) |> as.numeric()
      long.rw <- tidyr::pivot_longer(rw, cols = -year,
                                     names_to = "series" , values_to = "value")
      long.rw$type <- "rw"

      trans <- as.data.frame(trans)
      trans[, "year"] <- rownames(curv) |> as.numeric()
      long.trans <- tidyr::pivot_longer(trans, cols = -year,
                                        names_to = "series" , values_to = "value")
      long.trans$type <- "pwr.t_cu"

      curv <- as.data.frame(curv)
      curv[, "year"] <- rownames(curv) |> as.numeric()
      long.curv <- tidyr::pivot_longer(curv, cols = -year,
                                       names_to = "series" , values_to = "value")
      long.curv$type <- "pwr.t_cu"

      detr <- as.data.frame(detr)
      detr[, "year"] <- rownames(detr) |> as.numeric()
      long.detr <- tidyr::pivot_longer(detr, cols = -year,
                                       names_to = "series" , values_to = "value")
      long.detr$type <- "de"

      all.df <- rbind(long.rwl, long.trans, long.detr) |>
        na.omit() |> as.data.frame()

      all.df$year <- as.numeric(all.df$year)
      all.df$label <- ifelse(all.df$type %in% "rw",
                             paste("Original ring widths - series ID:", all.df$series),
                             ifelse(all.df$type %in% "pwr.t_cu",
                                    paste0("Transformed series with fitted trend (",detrend.method,")"),
                                    "Residual detrended series"))

      # add factor and year to the long.curv data.frame
      long.curv$type <- factor(long.curv$type, levels = c("pwr.t_cu"))
      long.curv$year <- as.numeric(long.curv$year)

      # merge the message data with the long.curv data.
      mess.df <- data.frame(series = names(messages),
                            message = messages,
                            type = factor("pwr.t_cu", levels = "pwr.t_cu"))

      long.curv <- merge(long.curv, mess.df, by = c("series","type"))

      all.list <- split(all.df, f = as.factor(all.df$series))
      curv.list <- split(long.curv, f = as.factor(long.curv$series))


      # Make the plots
      plot.list <- mapply(FUN = \(x, y) {

        max.trans <- x[x[, "type"] %in% "pwr.t_cu", "value"] |> max(na.rm = TRUE)
        this.order <- c(unique(x[, "label"])[startsWith(unique(x[, "label"]),
                                                        prefix = "Or")],
                        unique(x[, "label"])[startsWith(unique(x[, "label"]),
                                                        prefix = "Tr")],
                        unique(x[, "label"])[startsWith(unique(x[, "label"]),
                                                        prefix = "Re")])
        x[, "type"] <- factor(x[, "label"],
                              levels = this.order)
        y[, "type"] <- factor(unique(x[, "label"])[startsWith(unique(x[, "label"]),
                                                              prefix = "Tr")],
                              levels = this.order)

        z <- y[1,]
        z[, "value"] <- max.trans
        z[, "year"] <- median(y[, "year"])

        x_axis_params <- seq(min(na.omit(x[, "year"])),
                             max(na.omit(x[, "year"])),
                             length.out = 5) |>
          round(digits = -1)

        ggplot2::ggplot(x, aes(year, value, color = type)) +
          scale_color_manual(values = c("black","black","blue"),
                             guide = "none") +
          geom_line(linewidth = 0.25) +
          facet_wrap(~type,
                     ncol = 1,
                     scales = "free_y") +
          geom_line(data = y, color = "blue") +
          geom_text(data = z,
                    size = 3,
                    aes(label = message)) +
          xlab("Year") +
          scale_x_continuous(breaks = x_axis_params,
                             limits = range(x_axis_params)) +
          theme(strip.background = element_blank(),
                strip.text = element_text(hjust = 0),
                strip.clip = "off",
                axis.title.y = element_blank(),
                panel.grid = element_blank(),
                panel.background = element_blank())
      },
      x = all.list, y = curv.list,
      SIMPLIFY = FALSE)

      out.list <- list(detr[, !c(colnames(detr) %in% "year")],
                       curv[, !c(colnames(detr) %in% "year")],
                       trans[, !c(colnames(detr) %in% "year")],
                       messages,
                       plot.list)

    } else { # i.e., detrending == FALSE

      # convert all rwl format data.frames to long format, and then rbind them together.
      rw <- as.data.frame(rwl)
      rw[, "year"] <- rownames(rw) |> as.numeric()
      long.rw <- tidyr::pivot_longer(rw, cols = -year,
                                     names_to = "series" , values_to = "value")
      long.rw$type <- "rw"

      trans <- as.data.frame(trans)
      trans[, "year"] <- rownames(trans) |> as.numeric()
      long.trans <- tidyr::pivot_longer(trans, cols = -year,
                                        names_to = "series" , values_to = "value")
      long.trans$type <- "pwr.t"


      all.df <- rbind(long.rw, long.trans) |>
        na.omit() |> as.data.frame()

      all.df$year <- as.numeric(all.df$year)
      all.df$label <- ifelse(all.df$type %in% "rw",
                             paste("Original ring widths - series ID:", all.df$series),
                             "Transformed series")

      # merge the message data with the rest of the data.
      mess.df <- data.frame(series = names(messages),
                            message = messages)

      all.df <- merge(all.df, mess.df, by = c("series"))

      all.list <- split(all.df, f = as.factor(all.df$series))

      # Make the plots
      plot.list <- lapply(X = all.list,
                          FUN = \(x) {

                            max.trans <- x[x[, "type"] %in% "pwr.t", "value"] |> max(na.rm = TRUE)
                            this.order <- c(unique(x[, "label"])[startsWith(unique(x[, "label"]),
                                                                            prefix = "Or")],
                                            unique(x[, "label"])[startsWith(unique(x[, "label"]),
                                                                            prefix = "Tr")])
                            x[, "type"] <- factor(x[, "label"],
                                                  levels = this.order)

                            z <- x[startsWith(x[, "label"], prefix = "Tr"),][1,]
                            z[, "value"] <- max.trans
                            z[, "year"] <- median(x[, "year"])

                            x_axis_params <- seq(min(na.omit(x[, "year"])),
                                                 max(na.omit(x[, "year"])),
                                                 length.out = 5) |>
                              round(digits = -1)

                            ggplot2::ggplot(x, aes(year, value)) +
                              geom_line(linewidth = 0.25) +
                              facet_wrap(~type,
                                         ncol = 1,
                                         scales = "free_y") +
                              geom_text(data = z,
                                        size = 2.5,
                                        aes(label = message)) +
                              xlab("Year") +
                              scale_x_continuous(breaks = x_axis_params,
                                                 limits = range(x_axis_params)) +
                              theme(strip.background = element_blank(),
                                    strip.text = element_text(hjust = 0),
                                    strip.clip = "off",
                                    axis.title.y = element_blank(),
                                    panel.grid = element_blank(),
                                    panel.background = element_blank())
                          })


      out.list <- list(trans[, !c(colnames(trans) %in% "year")], messages, plot.list)

    }
    out.list

  } ## End of function


