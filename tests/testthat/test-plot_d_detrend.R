

test_that("function returns error for invalid data input", {
  expect_error(
    plot_d_detrend(x = "road to nowhere")
  )
})
