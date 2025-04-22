


# Need to generate some fake tree-ring data for testing
rwl <- matrix(nrow = 50, ncol = 25)
rwl <- apply(rwl, MARGIN = 2, FUN = \(x) runif(length(x), 0.1, 2)) |>
  as.data.frame()
colnames(rwl) <- paste0("t", 1:25)
rownames(rwl) <- 1951:2000
class(rwl) <- c("rwl","data.frame")

# xd_check <- function(data = NULL, # the data you are checking. long format or rwl
#                      ref = NULL, # a single reference data.frame or a list of references.
#                      std.method = "p2yrsL", # standardization method
#                      max.offset = 10, # max offset in years to check
#                      p.thresh = 0.05,
#                      out.format = "rwl"


test_that("function returns vector", {
  expect_vector(
    xd_check(data = rwl,
             ref = NULL,
             std.method = "p2yrsL",
             max.offset = 10,
             p.thresh = 0.05,
             out.format = "rwl")
  )
})



# Basic error catching

test_that("function returns error for invalid data input", {
  expect_error(
    xd_check(data = "road to nowhere",
             ref = NULL,
             std.method = "p2yrsL",
             max.offset = 10,
             p.thresh = 0.05,
             out.format = "rwl")
  )
})


test_that("function returns error for missing data input", {
  expect_error(
    xd_check(data = NULL,
             ref = NULL,
             std.method = "p2yrsL",
             max.offset = 10,
             p.thresh = 0.05,
             out.format = "rwl")
  )
})




test_that("function returns error for missing data input", {
  expect_error(
    xd_check(data = rwl,
             ref = NULL,
             std.method = "p2yrsL",
             max.offset = 10,
             p.thresh = 0.05,
             out.format = "rwl")
  )
})
