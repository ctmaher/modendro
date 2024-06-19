
####
test_that("Thows errors when fed gobbleygook arguments", {

  expect_error(
    ITRDB_search(species = "WHAT", # Species
                       lon.range = c(-165, -140),
                       lat.range = c(66, 70),
                       simplify = TRUE)
  )

  expect_error(
    ITRDB_search(species = "PCGL",
                       lon.range = c(-16500, -140), # lon.range
                       lat.range = c(66, 70),
                       simplify = TRUE)
  )

  expect_error(
  ITRDB_search(species = "PCGL",
                       lon.range = c(-165, -1400), # lon.range
                       lat.range = c(66, 70),
                       simplify = TRUE)
  )

  expect_error(
  ITRDB_search(species = "PCGL",
                       lon.range = c(-165, -140),
                       lat.range = c(660, 70), # lat.range
                       simplify = TRUE)
  )

  expect_error(
  ITRDB_search(species = "PCGL",
                       lon.range = c(-165, -140),
                       lat.range = c(66, 700), # lat.range
                       simplify = TRUE)
  )

  expect_error(
    ITRDB_search(species = "PCGL",
                 lon.range = c(-165, -140),
                 lat.range = c(66, 70),
                 limit = "No limits bruh", # limit
                 simplify = TRUE)
  )

  expect_error(
    ITRDB_search(species = "PCGL",
                 lon.range = c(-165, -140),
                 lat.range = c(66, 70),
                 limit = 20,
                 simplify = "KISS") # simplify
  )

})

####
test_that("Essential function works - and returns a vector", {
expect_vector(
  ITRDB_search(species = "PCGL",
                           lon.range = c(-160, -150),
                           lat.range = c(66, 70),
                           limit = 5,
                           simplify = TRUE)
  )
})

####
test_that("Species outside of geographic range just returns an error", {
  expect_error(
    ITRDB_search(species = "PILO",
                             lon.range = c(-160, -150),
                             lat.range = c(66, 70),
                             limit = 5,
                             simplify = TRUE)
    )
})

####
test_that("Including species inside & outside of geographic range still returns a vector", {
  expect_vector(
    ITRDB_search(species = c("PILO","PCGL"),
                 lon.range = c(-160, -150),
                 lat.range = c(66, 70),
                 limit = 5,
                 simplify = TRUE)
  )
})

####
test_that("A big search still returns a vector", {
  expect_vector(
    ITRDB_search(species = c("PICO","PISY","PIRE","PISI"), # Widespread pine species
                 lon.range = c(-180, 180), # the whole globe
                 lat.range = c(-80, 80), # all of the areas with trees, and then some
                 limit = 1000,
                 simplify = TRUE)
  )
})
