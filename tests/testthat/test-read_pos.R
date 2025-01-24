test_that("function returns errors for invalid path input", {
  expect_error(
  read_pos(path = "road to nowhere")
  )
})


test_that("function returns errors for invalid path input", {
  expect_error(
    read_pos(path = 42)
  )
})


test_that("function returns errors for invalid path input", {
  expect_error(
    read_pos()
  )
})


test_that("function returns errors for invalid default.OD input", {
  expect_error(
    read_pos(path = system.file("extdata", package = "modendro"), default.OD = "queseyo")
  )
})


test_that("function returns errors for invalid default.OD input", {
  expect_error(
    read_pos(path = system.file("extdata", package = "modendro"), default.OD = c(1919,2021))
  )
})
