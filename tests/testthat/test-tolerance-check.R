test_that("Check whether the tolerance check to determine the end of the MM algorithm works", {
  expect_true(check_has_converged(old_value = 1, new_value = 1.00001, tolerance = 0.1))
  expect_false(check_has_converged(old_value = 1, new_value = 2, tolerance = 0.1))
})