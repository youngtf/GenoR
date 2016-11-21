library(GenoR)
context("test the behavious about bin_list")

# test_that("", {                      ## the description of the section
#   expect_equal(var_1, var_2)
#   expect_match()
#   expect_output()
#   expect_message()
#   expect_warning()
#   expect_error()
#   expect_is()
# })

# test_that("genome-wide bin_list creation", {
#   
# })



test_that("test the checking of bin_list object", {
  # bin_list = genome_wide_bin_list(sorted_map)
  # expect_is(check_bin_list(bin_list, sorted_map),class = "data.frame")
  bin_list = data.frame(Chr = "N", start_pos = 1, end_pos = 35)
  expect_error(check_bin_list(bin_list, sorted_map))
  bin_list = data.frame(Chr = "1", start_pos = 1000000000, end_pos = 35)
  expect_error(check_bin_list(bin_list, sorted_map))
  bin_list = data.frame(Chr = "1", start_pos = 100, end_pos = 35)
  expect_error(check_bin_list(bin_list, sorted_map))
  bin_list = data.frame(Chr = "1", start_pos = 1, end_pos = 350000000)
  expect_warning(check_bin_list(bin_list, sorted_map))
})
