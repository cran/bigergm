test_that("Relabeling to always provide the canonical labels works", {
  block <- factor(c(1,4,4,3,2)) 
  expect_equal(relabel(block = block,n_blocks = 4), factor(c(1,2,2,3,4)))
})
