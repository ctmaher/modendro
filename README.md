
<!-- README.md is generated from README.Rmd. Please edit that file -->

# modendro: more dendro (tree-ring) functions

<!-- badges: start -->
<!-- badges: end -->

modendro is a collection of functions to do a variety of tree-ring
analyses that are not available elsewhere or are improvements on
existing analyses: various utility functions for wrangling tree-ring
data, power transformation & detrending, flexible climate-growth
correlations to identify the strongest seasonal climate signals, and
disturbance detection & removal.

My inclusion of an analysis method in modendro does not imply broad
endorsement of the method. My goal is to make available some of these
methods that I’ve found useful in the past but had to write my own code
for. Enjoy and I hope some of these are useful for you too!

## Installation

You can install the development version of modendro from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ctmaher/modendro")
```

I recommend you do this often, as modendro is under active development.

## Example: convert rwl-format data.frames to “long” format - and back

The rwl-format data.frame (rownames are years, colnames are series
names) is the standard for the important dplR package, and modendro
inherits this. However, this format is inconvenient for many other
standard operations in R - a “long” format would be much more useful.
Enter modendro’s `rwl_longer()` function.

``` r
library(modendro)
#> Registered S3 method overwritten by 'quantmod':
#>   method            from
#>   as.zoo.data.frame zoo
data("PerkinsSwetnam96")
head(PerkinsSwetnam96)[,1:5] # typical rwl-format
#>     TWP05 TWP07 TWP08 TWP13 TWP14
#> 726    NA    NA    NA    NA    NA
#> 727    NA    NA    NA    NA    NA
#> 728    NA    NA    NA    NA    NA
#> 729    NA    NA    NA    NA    NA
#> 730    NA    NA    NA    NA    NA
#> 731    NA    NA    NA    NA    NA


PS.long <- rwl_longer(rwl = PerkinsSwetnam96,
                      series.name = "series", # name the series IDs whatever you want
                      dat.name = "rw.mm", # same for the measurement data
                      trim = TRUE, # trim off the NAs before and after a series
                      new.val.internal.na = NULL) # leave NAs internal to the measurement series
# Note that we could replace internal NAs here if there were any (e.g., with 0s)

head(PS.long) # familiar format
#>             year series rw.mm
#> RRR05.58876 1319  RRR05  0.42
#> RRR05.58877 1320  RRR05  0.46
#> RRR05.58878 1321  RRR05  0.34
#> RRR05.58879 1322  RRR05  0.32
#> RRR05.58880 1323  RRR05  0.42
#> RRR05.58881 1324  RRR05  0.33
```

Once our tree-ring data is in the familiar (to me anyway) long format,
we do all kinds of data wrangling, analyses and plotting that we
normally do in R. For example, we might want to merge the tree-ring data
with site IDs so that we can add site as a random effect in a model or
for plotting.

``` r

data("PSgroupIDs")

PS.long.sites <- merge(PS.long, PSgroupIDs, by = "series")

library(ggplot2)

ggplot(PS.long.sites, aes(year, rw.mm, group = series)) +
  geom_line(alpha = 0.5) +
  facet_wrap( ~ site, ncol = 1)
```

<img src="man/figures/README-example 1.2-1.png" width="100%" />

Another useful application of `rwl_longer()` is when we have several
rwl-format files that we want to bind together. This is not intuitive
with the rwl-format - the number of rows may not line up and the columns
won’t either. This is simple with the long format.

We’ll borrow some data from the dplR package for this.

We may want to convert this joined data back into a single rwl-format
data.frame, e.g., for export or for unifying analyses in modendro or
dplR. We can do this with modendro’s `longer_rwl()`.

``` r

library(dplR)
#> This is dplR version 1.7.7.
#> dplR is part of openDendro https://opendendro.org.
#> New users can visit https://opendendro.github.io/dplR-workshop/ to get started.
data("ca533")

ca533.long <- rwl_longer(ca533,
                         series.name = "series",
                         dat.name = "rw.mm")

comb.long <- rbind(PS.long, ca533.long)

# note that longer_rwl() doesn't care if you have other variables in your long data.frame
# We'll add some here to illustrate that
comb.long$foo <- "bar"
comb.long$bar <- "foo"

comb.rwl <- longer_rwl(df = comb.long,
                       series.name = "series",
                       dat.name = "rw.mm")

head(comb.rwl)[,1:5]
#>     RRR05 RRR06 RRR07 RRR15 RRR19
#> 626    NA    NA    NA    NA    NA
#> 627    NA    NA    NA    NA    NA
#> 628    NA    NA    NA    NA    NA
#> 629    NA    NA    NA    NA    NA
#> 630    NA    NA    NA    NA    NA
#> 631    NA    NA    NA    NA    NA
```

## Example: power transformation & detrending ala Cook & Peters (1997)

Cook & Peters describe a method in their 1997 paper in The Holocene for
removing age/size trends from tree-ring series in a way that minimizes
bias that can result from using more traditional methods. Specifically,
C&P’s method involves power transforming raw ring widths to homogenize
variance, fitting a trend, and then *subtracting* the trend line to get
power transformed detrended residual series.

In modendro, the `cp_detrend()` function does this. I should say that
you can do this with a string of dplR functions, but `cp_detrend()`
contains the entire operation in a single function. The trend fitting
and subtraction step is done using dplR’s `detrend()` function, allowing
all the detrend curve options availble. Additionally, `cp_detrend()` is
much more transparent and returns information about the power
transformation, including the power used. Another reason for this
detailed output is that this is the process used in the modendro
function `ci_detect()`, which I’ll demonstrate below. The `cp_detrend()`
process can be visualized with the `plot_cp_detrend()` function.

``` r

PS.cp <- cp_detrend(PerkinsSwetnam96, detrend.method = "ModNegExp")
names(PS.cp) # output is a list
#> [1] "Resid. detrended series" "Detrend curves"         
#> [3] "Transformed ring widths" "Transformation metadata"
#> [5] "Detrending metadata"     "Raw ring widths"

PS.cp[["Transformation metadata"]][["RRR19"]]
#>       series optimal.pwr            action
#> RRR19  RRR19   0.2308274 Power transformed

PS.cp.plots <- plot_cp_detrend(PS.cp)

PS.cp.plots[["RRR19"]]
```

<img src="man/figures/README-example 2.1-1.png" width="100%" />

## Example: flexible growth-climate relationships

`n_mon_corr()` is a modendro function that is for exploratory data
analysis of growth-climate relationships using tree-rings and monthly
climate data. The concept is relatively simple, in practice it is not so
simple!

The basic operation is to take monthly climate data and compute every
possible contiguous monthly aggregate window given a user-specified n
months (2-12 month aggregates possible) and for a set of lag years. The
monthly aggregate windows are continuous across (lagged) annual
boundaries. Then, `n_mon_corr()` computes correlations with each of
these monthly aggregates and every single tree-ring series in your rwl.
As such, this is a “brute-force” type of operation and can be somewhat
slow, depending on the parameters.

`n_mon_corr()` works in both hemispheres - if hemisphere = “S” is
selected, then a +1 “lag” is added to whatever lags you specify (e.g.,
-2, -1, 0, +1). This is because the convention of assigning years to
tree-rings in the southern hemisphere is to use the year that growth
starts in. Since growth will usually continue through the Gregorian new
year, the +1 is there to accommodate the full range of possible monthly
climate aggregates into the growing season.

I have built in operations like “prewhitening” so that the both the
climate aggregates (generated inside the function) and the tree-ring
data get treated the same way. In most tree-ring analyses, prewhitening
means residuals from an autoregressive model of the ring-width series,
the goal of which is to remove autocorrelation (includes trends) and
leave only the high-frequency variability that resembles white noise. In
practice, I have found that the simple `ar()` is not always adequate to
model all of the autocorrelation in tree-ring or climate series. A more
complex time-series model does better, like an ARIMA model.
`n_mon_corr()` uses `auto.arima()` from the `forecast` package to find a
best-fit ARIMA model, and then extracts the residuals. ARIMA models
don’t model the variability, so all of the high-frequency component is
retained, including any changes in variability over time (i.e.,
volatility). All this is to say that `n_mon_corr()` does a more honest
job of removing autocorrelation than traditional tree-ring approaches.
We want this before doing cross-correlations between two time series!

The default method for correlations is a Spearman rank correlation -
this is to allow for the detection of possibly non-linear associations
between tree-rings and climate. For significance testing of the
individual correlation analyses, I use a method described by Lun et
al. (2022), Journal of Applied Statistics, 50(14) that adjusts the
significance assuming that there is autocorrelation in the time-series.
For “prewhitened” series, this will probably be a bit conservative, but
is useful if you want to run both “prewhitened” and not series to try to
understand the influence that trends may have on the relationships - as
such it makes sense that both types of series are treated with the same
significance test. Hence the default. Kendall rank correlations are also
possible and use this adjusted significance test. *Pearson correlations
are an option, but there is no adjustment for autocorrelation* in
Pearson. I don’t recommend it.

If you turn off “prewhitening”, be aware of the effects of
autocorrelated time-series in cross-correlations. This is a problem for
the *correlation coefficients*, not the significance testing per se. See
Yule (1926), Journal of the Royal Statistical Society 89(1). The take
away is that for cross-correlations between highly autocorrelated time
series (i.e., low-frequency amplitudes are much greater than
high-frequency amplitudes), correlation coefficients are nearly always
very high, even though they may be meaningless.

The `n_mon_corr()` documentation has more details on the various
arguments.

``` r
# Bring in some climate data associated with the tree ring data we loaded earlier
data("idPRISM")
head(idPRISM)
#>           Date year month    PPT.mm     Tavg.C
#> 1   1895-01-15 1895     1 219.04667  -8.733333
#> 106 1903-02-15 1903     2  15.02333 -12.866667
#> 211 1911-03-15 1911     3  50.58667  -4.066667
#> 316 1919-04-15 1919     4  87.31667   0.800000
#> 421 1927-05-15 1927     5  98.99000   2.400000
#> 526 1935-06-15 1935     6  33.36333   8.300000

PS.corr <- n_mon_corr(rwl = PerkinsSwetnam96, 
                      clim = idPRISM, 
                      clim.var = "Tavg.C",
                      common.years = 1895:1991,
                      agg.fun = "mean",
                      max.win = 6,
                      win.align = "left",
                      max.lag = 2,
                      hemisphere = "N",
                      prewhiten = TRUE,
                      corr.method = "spearman")
#> The following tree-ring series have < 4 years overlap with clim data
#>     and will be removed from rwl:
#> TWP05, TWP13, TWP20, TWP21, SDP36, SDP26, SDP28, SDM10, UPS12, UPS16, UPS31, UPS35, UPS04, UPSM1, UPM5, UPM14, UPM25, RRR19, RRR26, RRR28
#> The following tree-ring series have < 25 years overlap with clim data.
#>     Interpret correlations cautiously.
#> UPM3, UPM12, UPM16, UPM09

names(PS.corr)
#> [1] "Correlation results"             "Climate data (prewhitened)"     
#> [3] "Climate data (raw)"              "Ring-width series (prewhitened)"
```

Note the warnings about series that don’t have enough overlap with the
climate data.

The output of `n_mon_corr()` includes the individual correlation
results, the monthly aggreated climate data (prewhitened & raw), and the
prewhitened ring widths. This allows you to inspect each of these
elements.

As with other `modendro` functions, there is a separate plotting
function. This is designed to deal with multiple (a lot of) ring width
series, so keep that in mind.

The first plot shows the percentage of significant correlations for each
month and each window length, the second shows the mean correlation
coefficients of the same. These are two ways to find the strongest
signals across lags and moving window lengths. Lags are indicated by
labeled rectangles at the bottom of the plots (-2, -1, 0).

I include the aggregated data used to make the plots as well, for some
added transparency.

``` r

PS.corr.plots <- plot_n_mon_corr(PS.corr)

names(PS.corr.plots)
#> [1] "Percent sig. corr. plot" "Mean corr. coef. plot"  
#> [3] "Aggregated data"

PS.corr.plots[["Percent sig. corr. plot"]]
```

<img src="man/figures/README-example 3.2-1.png" width="100%" />

``` r
PS.corr.plots[["Mean corr. coef. plot"]]
```

<img src="man/figures/README-example 3.2-2.png" width="100%" />

``` r
PS.corr.plots[["Aggregated data"]] |> head()
#>    month win.len lag  dir prop.sig   mean.coef
#> 3      1       1  -2 Neg.  0.00000 -0.06144863
#> 4      1       1  -2 Pos. 10.81081  0.15089999
#> 9      1       2  -2 Neg.  0.00000 -0.03508997
#> 10     1       2  -2 Pos. 21.62162  0.21733757
#> 15     1       3  -2 Neg.  0.00000 -0.03507990
#> 16     1       3  -2 Pos. 21.62162  0.22710420
```

## Example: detection and removal of disturbances in tree-ring series

Detection of release or suppression events…
