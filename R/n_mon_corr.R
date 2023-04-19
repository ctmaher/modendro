#' Growth-climate correlations via moving windows
#'
#' @description
#' Function to take a tree ring chronology and a monthly climate variable to compute
#' monthly aggregate correlations for every combination of consecutive months for calendar years and water years (Oct-Sept)
#' The function computes the correlations for a specified number lagged years (i.e., years before ring formation).
#' This function is designed for exploratory data analysis on growth-climate relationships in trees.
#' 3 years is the default and is a sensible maximum for most analyses.
#'
#' @details
#' This
#'
#'
#' @param chrono a chron object (such as that produced by dplR's `chron()`).
#' @param clim a `data.frame` with at least 3 columns: year, month (numeric), and a climate variable.
#' @param var character vector - the colname of the climate variable of interest in the `clim` data.frame.
#' @param chrono.col character vector - the colname of the chronology series (default is "xxxstd", which is the defualt produced by dplR's `chron()`).
#' @param agg.fun character vector specifying the function to use for aggregating monthly climate combinations.
#' Options are "mean" or "sum", e.g., for temperature or precipitation data, respectively. Default is "mean".
#' @param max.lag numeric vector specifying how many years of lag to calculate calculations for. Default is 3 years.
#' @param corr.method character vector specifying which correlation method to use. Passes to `cor.test()`.
#' @param chrono.name character vector - the name of your chronology. This is used in the title of your plot.
#' @param plots logical vector indicating whether or not to produce plots.
#'
#' @return A 1-2 element list containing data.frames of the correlation results and the default plots of the same data.
#'
#' @import plyr
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' Will add some later
#'



n_mon_corr <- function(chrono = NULL, clim = NULL, var = NULL,
                       chrono.col = "xxxstd", agg.fun = "mean",
                       max.lag = 3, corr.method = "pearson",
                       chrono.name = NULL, plots = TRUE){

  if (chrono.col %in% colnames(chrono)){

  } else {
    stop("Please supply the chrono.col argument with the correct name
         of the growth variable in the chronology data.frame")
  }


  # make sure that "year" columns are labelled as such
  colnames(clim)[which((substr(colnames(clim), start = 1, stop = 2)
                        %in% c("Ye","ye")) == T)] <- "year"
  colnames(chrono)[which((substr(colnames(chrono), start = 1, stop = 2)
                          %in% c("Ye","ye")) == T)] <- "year"
  # same for month
  colnames(clim)[which((substr(colnames(clim), start = 1, stop = 2)
                        %in% c("Mo","mo")) == T)] <- "month"

  # make sure year and month are integers
  clim$year <- as.integer(clim$year)
  chrono$year <- as.integer(chrono$year)
  clim$month <- as.integer(clim$month)

  # add water year for Oct-Dec
  clim$water.year <- ifelse(clim$month %in% c(10,11,12), clim$year + 1, clim$year)

  # Find min and max years for the correlations - for the plots later
  min.y <- max(min(chrono$year), min(clim$year))
  max.y <- min(max(chrono$year), max(clim$year))

  # Create a vector of all possible combinations of months
  # Regular calendar year first
  mos.mat <- expand.grid(1:12,1:12)
  mos <- apply(mos.mat, MARGIN = 1, FUN = function(x) seq(from = x[1], length.out = x[2]))
  unreal <- lapply(mos, FUN = function(x) x < 13)
  unreal2 <- lapply(unreal, FUN = function(x) length(x[x == F]))
  mos <- mos[which(unreal2 < 1)]

  # Now water year months
  wy <- c(10:12,1:9)
  wy.list1 <- vector("list", 9) # 9 positions x 3 starting points
  new.pos <- 4:12
  for (i in seq_along(new.pos)){
    wy.list1[[i]] <- wy[1:new.pos[i]]
  }
  wy.list2 <- vector("list", 9) # 9 positions x 3 starting points
  for (i in seq_along(new.pos)){
    wy.list2[[i]] <- wy[2:new.pos[i]]
  }
  wy.list3 <- vector("list", 9) # 9 positions x 3 starting points
  for (i in seq_along(new.pos)){
    wy.list3[[i]] <- wy[3:new.pos[i]]
  }
  # Lastly cover the single months and 10:12
  wy.list0 <- list(10,11,12,c(10,11),c(11,12),c(10,11,12))
  # put them together
  wy.mos <- c(wy.list0,wy.list1,wy.list2,wy.list3)


  # How to build in the annual lags?
  # Perhaps I can hold the data.frame of results in a list, with the length of the list being equal to
  # the max lag?
  lag.seq <- 0:max.lag
  #lag.list <- vector("list", length(lag.seq))
  #for (y in seq_along(lag.seq)){
  lag.list <- ldply(lag.seq, .fun = function(l){
    # Challenge: vectorize the code below. The order of month groupings doesn't really matter,
    # So I could make this a lot faster by not doing a loop.

    # calendar year
    cal.cor.res <- ldply(mos, .fun = function(x){
      # Aggregate the variable of interest for the given month sequence
      clim.mo <- aggregate(formula(paste(var, "year", sep = "~")), data = clim[clim$month %in% x,],
                           FUN = function(z){ifelse(agg.fun %in% "mean", mean(z), sum(z))})
      # attach the chronology to the climate data - enter the current lag before this step
      chrono.new <- chrono
      chrono.new$year <- chrono.new$year - l
      clim.mo <- merge(clim.mo, chrono.new, by.x = "year", by.y = "year")
      # Run the correlation test
      # Need a reliable way to select the treering chronology variable
      ct <- cor.test(clim.mo[,var], clim.mo[, chrono.col], method = corr.method)
      # put the results together in a data.frame
      if (corr.method %in% "pearson") {
        result <- data.frame(months = ifelse(length(x) > 1, deparse(x), paste(x)),
                             coef = round(ct$estimate[[1]],3),
                             p = ct$p.value[[1]],
                             ci.lo = ct$conf.int[1],
                             ci.hi = ct$conf.int[2])
      } else {
        result <- data.frame(months = ifelse(length(x) > 1, deparse(x), paste(x)),
                             coef = round(ct$estimate[[1]],3),
                             p = ct$p.value[[1]])
      }
      # return the result
      result
    })
    cal.cor.res$sig <- ifelse(cal.cor.res$p > 0.05, "",
                              ifelse(cal.cor.res$p <= 0.05 & cal.cor.res$p >= 0.01, "*",
                                     ifelse(cal.cor.res$p <= 0.001, "***","**")))
    cal.cor.res$seas.win <- "Calendar year"


    # water year
    wat.cor.res <- ldply(wy.mos, .fun = function(x){
      # Aggregate the variable of interest for the given month sequence
      clim.mo <- aggregate(formula(paste(var, "water.year", sep = "~")), data = clim[clim$month %in% x,],
                           FUN = function(z){ifelse(agg.fun %in% "mean", mean(z), sum(z))})
      colnames(clim.mo)[1] <- "year"
      # attach the chronology to the climate data - enter the current lag before this step
      chrono.new <- chrono
      chrono.new$year <- chrono.new$year - l
      clim.mo <- merge(clim.mo, chrono.new, by.x = "year", by.y = "year")
      # Run the correlation test
      ct <- cor.test(clim.mo[,var], clim.mo[, chrono.col], method = corr.method)
      # put the results together in a data.frame
      if (corr.method %in% "pearson") {
        result <- data.frame(months = ifelse(length(x) > 1, paste(x[1], x[length(x)], sep = ":"), paste(x)),
                             coef = round(ct$estimate[[1]],3),
                             p = ct$p.value[[1]],
                             ci.lo = ct$conf.int[1],
                             ci.hi = ct$conf.int[2])
      } else {
        result <- data.frame(months = ifelse(length(x) > 1, paste(x[1], x[length(x)], sep = ":"), paste(x)),
                             coef = round(ct$estimate[[1]],3),
                             p = ct$p.value[[1]])
      }
      # return the result
      result
    })
    wat.cor.res$sig <- ifelse(wat.cor.res$p > 0.05, "",
                              ifelse(wat.cor.res$p <= 0.05 & wat.cor.res$p >= 0.01, "*",
                                     ifelse(wat.cor.res$p <= 0.001, "***","**")))
    wat.cor.res$seas.win <- "Water year"

    cor.results <- rbind(cal.cor.res, wat.cor.res)
    cor.results$lag <- ifelse(l == 0, paste(l),paste0("-",l))
    cor.results
  })

  lag.list$lag <- factor(lag.list$lag, levels = lag.seq*-1)
  lag.list <- lag.list[order(lag.list$coef, decreasing = T),] # sort by correlation coef
  # There may be better ways to organize this that allow a breakdown of which monthly aggregates
  # drive the strong correlations. Right now it is a lot to look at.


  if (plots == TRUE) {

    plot.df <- lag.list
    for (i in 1:nrow(plot.df)){
      plot.df$start.mo[i] <- as.numeric(strsplit(plot.df$months, ":")[[i]][1])
    }
    for (i in 1:nrow(plot.df)){
      plot.df$end.mo[i] <- as.numeric(ifelse(length(strsplit(plot.df$months, ":")[[i]]) == 2,
                                             strsplit(plot.df$months, ":")[[i]][2],
                                             strsplit(plot.df$months, ":")[[i]][1]))
    }

    plot.df$start.mo <- ifelse(plot.df$seas.win %in% "Water year" & plot.df$start.mo %in% 12, 0, plot.df$start.mo)
    plot.df$start.mo <- ifelse(plot.df$seas.win %in% "Water year" & plot.df$start.mo %in% 11, -1, plot.df$start.mo)
    plot.df$start.mo <- ifelse(plot.df$seas.win %in% "Water year" & plot.df$start.mo %in% 10, -2, plot.df$start.mo)
    plot.df$end.mo <- ifelse(plot.df$seas.win %in% "Water year" & plot.df$end.mo %in% 12, 0, plot.df$end.mo)
    plot.df$end.mo <- ifelse(plot.df$seas.win %in% "Water year" & plot.df$end.mo %in% 11, -1, plot.df$end.mo)
    plot.df$end.mo <- ifelse(plot.df$seas.win %in% "Water year" & plot.df$end.mo %in% 10, -2, plot.df$end.mo)

    plot.df$simp.sig <- ifelse(plot.df$sig %in% c("*","**","***"), "Sig.","Not sig.")
    mo.xaxis <- c(10:12,1:12)
    title <- ifelse(is.null(chrono.name), paste0("Chronology correlations with ",var),
                    paste0(chrono.name, " chronology correlations with ",var))

    out.plot <- ggplot(plot.df, aes(start.mo, coef, color = simp.sig)) +
      geom_point(shape = 5, size = 1.5) +
      geom_point(data = plot.df, aes(end.mo, coef, color = simp.sig), shape = 124, size = 2) +
      scale_x_continuous(breaks = c(-2:12), labels = mo.xaxis) +
      scale_y_continuous(breaks = seq(-1,1, by = 0.1)) +
      geom_segment(aes(xend = end.mo, yend = coef)) + xlab("Month") +
      ylab(paste0("Correlation coefficient\n(",corr.method,")")) +
      scale_color_manual(name = "", values = c("grey80","black")) +
      theme_bw() + facet_wrap(~lag, ncol = 1) +
      ggtitle(label = title,
              subtitle = paste0("Monthly climate ", agg.fun, "s with annual lag 0:", max.lag," (", min.y,"-",max.y,")"))

    out.list <- list(lag.list, out.plot)
    return(out.list)
  } else {
    return(lag.list)
  }

} # End of function
