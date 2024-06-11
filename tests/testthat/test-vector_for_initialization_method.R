test_that("Providing a vector as the initialization works", {
  data(toyNet)
  
  # Specify the model that you would like to estimate.
  model_formula <- toyNet ~ edges + nodematch("x") + nodematch("y") + triangle
  # Estimate the model
  bigergm_res <- bigergm(
    object = model_formula,
    # The model you would like to estimate
    n_blocks = 4,
    # The number of blocks
    n_MM_step_max = 1,
    # The maximum number of MM algorithm steps
    estimate_parameters = TRUE,
    # Perform parameter estimation after the block recovery step
    clustering_with_features = TRUE,
    # Indicate that clustering must take into account nodematch on characteristics
    check_block_membership = FALSE)

  res <- bigergm(
    object = model_formula,
    # The model you would like to estimate
    n_blocks = 4,
    # The number of blocks
    n_MM_step_max = 1,
    # The maximum number of MM algorithm steps
    estimate_parameters = TRUE,
    # Perform parameter estimation after the block recovery step
    clustering_with_features = TRUE,
    # Indicate that clustering must take into account nodematch on characteristics
    check_block_membership = FALSE, 
    initialization = bigergm_res$block)

  expect_equal(res$initial_block, bigergm_res$block, check.attributes = FALSE)
  
  res <- bigergm(
    object = model_formula,
    # The model you would like to estimate
    n_block = 4,
    # The number of blocks
    n_MM_step_max = 1,
    # The maximum number of MM algorithm steps
    estimate_parameters = TRUE,
    # Perform parameter estimation after the block recovery step
    clustering_with_features = TRUE,
    # Indicate that clustering must take into account nodematch on characteristics
    check_block_membership = FALSE, 
    initialization = letters[1:4][bigergm_res$block])

  expect_equal(res$initial_block, bigergm_res$block, check.attributes = FALSE)
})
