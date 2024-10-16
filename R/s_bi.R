#' Calculate the biweight scale
#'
#' @description
#' This function calculates the biweight scale (s_bi) - a robust variance estimator based on the
#' biweight robust mean.
#' From Hoaglin et al. 1983 p.417, eq. 3.
#' This is to support the cp_detrend function, following the methods of Cook & Peters 1997.
#'
#' @param x a numeric vector
#'
#'
#' @return a 2-column, 1-row data.frame with 2 estimates of s_bi: the main one from the Hoaglin et
#' al book & the Kafadar (1979) version that the authors also mention.
#'
#' @references
#' Hoaglin, D. C., F. Mosteller, and J. W. Tukey. 1983. Understanding robust and exploratory data
#' analysis. New York: Wiley.
#'
#' Cook, E. R., and Peters, K. (1997) Calculating unbiased tree-ring indices for the study of
#' climatic and environmental change.
#' \emph{The Holocene}, \strong{7}(3), 361-370.
#'
#' @importFrom DescTools TukeyBiweight
#' @importFrom DescTools MAD
#' @export
#'
#' @examples
#' s_bi(rnorm(100))
#'

s_bi <- function(x) {

  stopifnot("x is not a numeric vector" =
              is.numeric(x) | is.double(x))

  x_tbrm <- DescTools::TukeyBiweight(x)
  u <- (x - x_tbrm)/(9*DescTools::MAD(x - x_tbrm))
  only_these <- which(abs(u) < 1)
  u <- u[only_these]
  x <- x[only_these]
  num1 <- length(x)^0.5
  num2 <- sum((((x - x_tbrm)^2) * (1 - u^2)^4))^0.5
  den <- sum((1 - u^2) * (1 - 5*(u^2)))

  s_bi <- (num1 * num2) /
    abs(den)

  s_bi_kafadar <- (num1 * num2) /
    ((abs(den)*(-1 + abs(den)))^0.5)

  data.frame(s_bi = s_bi, s_bi_kafadar = s_bi_kafadar)

}
