# Need to generate some fake tree-ring data for testing
rwl <- matrix(nrow = 50, ncol = 25)
rwl <- apply(rwl, MARGIN = 2, FUN = \(x) runif(length(x), 0.1, 2)) |>
  as.data.frame()
colnames(rwl) <- paste0("t", 1:25)
rownames(rwl) <- 1951:2000
class(rwl) <- c("rwl","data.frame")

# d_detrend <- function(data = NULL,
#                       win.len = 15,
#                       pgc.thresh = 50,
#                       d.detrend.method = "AgeDepSpline",
#                       detrend.method = "AgeDepSpline",
#                       nyrs = c(10, 30),
#                       event.type = "release")


test_that("function returns vector", {
  expect_vector(
    d_detrend(data = rwl,
              win.len = 15,
              pgc.thresh = 50,
              d.detrend.method = "AgeDepSpline",
              detrend.method = "AgeDepSpline",
              nyrs = c(10, 30),
              event.type = "release")
  )
})



# Basic error catching

test_that("function returns error for invalid data input", {
  expect_error(
    d_detrend(data = "road to nowhere",
             win.len = 15,
             pgc.thresh = 50,
             d.detrend.method = "AgeDepSpline",
             detrend.method = "AgeDepSpline",
             nyrs = c(10, 30),
             event.type = "release")
  )
})


test_that("function returns error for missing data input", {
  expect_error(
    d_detrend(data = NULL,
              win.len = 15,
              pgc.thresh = 50,
              d.detrend.method = "AgeDepSpline",
              detrend.method = "AgeDepSpline",
              nyrs = c(10, 30),
              event.type = "release")
  )
})



test_that("function returns error for invalid win.len input", {
  expect_error(
    d_detrend(data = rwl,
              win.len = "15",
              pgc.thresh = 50,
              d.detrend.method = "AgeDepSpline",
              detrend.method = "AgeDepSpline",
              nyrs = c(10, 30),
              event.type = "release")
  )
})


test_that("function returns error for invalid pgc.thresh input", {
  expect_error(
    d_detrend(data = rwl,
              win.len = 15,
              pgc.thresh = "50",
              d.detrend.method = "AgeDepSpline",
              detrend.method = "AgeDepSpline",
              nyrs = c(10, 30),
              event.type = "release")
  )
})


test_that("function returns error for invalid d.detrend.method input", {
  expect_error(
    d_detrend(data = rwl,
              win.len = 15,
              pgc.thresh = 50,
              d.detrend.method = "The best spline in the world",
              detrend.method = "AgeDepSpline",
              nyrs = c(10, 30),
              event.type = "release")
  )
})


test_that("function returns error for invalid detrend.method input", {
  expect_error(
    d_detrend(data = rwl,
              win.len = 15,
              pgc.thresh = 50,
              d.detrend.method = "AgeDepSpline",
              detrend.method = 1,
              nyrs = c(10, 30),
              event.type = "release")
  )
})


test_that("function returns error for invalid nyrs input", {
  expect_error(
    d_detrend(data = rwl,
              win.len = 15,
              pgc.thresh = 50,
              d.detrend.method = "AgeDepSpline",
              detrend.method = "AgeDepSpline",
              nyrs = 1,
              event.type = "release")
  )
})


test_that("function returns error for invalid event.type input", {
  expect_error(
    d_detrend(data = rwl,
              win.len = 15,
              pgc.thresh = 50,
              d.detrend.method = "AgeDepSpline",
              detrend.method = "AgeDepSpline",
              nyrs = c(10, 30),
              event.type = "hurricane")
  )
})


