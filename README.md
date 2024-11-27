
<!-- README.md is generated from README.Rmd. Please edit that file -->

# modendro: more dendro (tree-ring) functions

<!-- badges: start -->
<!-- badges: end -->

modendro is a collection of functions to do a variety of tree-ring
analyses that are not available elsewhere or are improvements on
existing analyses: power transformation & detrending, disturbance
detection, flexible climate-growth correlations to identify the
strongest seasonal climate signals, & various utility functions for
wrangling tree-ring data. The functions in modendro are designed to work
with the common tree-ring measurement formats that have years as
rownames and series names as colnames (i.e., the resulting format from
the `read.` functions in the dplR package). Within R, this format is
just a particularly-structured data.frame. The modendro functions don’t
“know” what file format the measurements originated from, and the only
things that matter are the rownames and the colnames

## Installation

You can install the development version of modendro from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ctmaher/modendro")
```

## Example: convert rwl-format data.frames to “long” format - and back

The rwl-format data.frame is the standard for the important dplR
package, and modendro inherits this. However, this format is
inconvenient for many other standard operations in R - a “long” format
would be much more useful. Enter modendro’s `rwl_longer()` function

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

<img src="man/figures/README-example 2-1.png" width="100%" />

More to come…
