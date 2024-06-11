test_that("control.ergm settings can be passed to the within estimation from hergm function", {
  # Define some settings for testing purposes
  test_burnin <- 9797
  test_interval <- 3434
  test_method <- 'MCMLE'

  bigergm_formula <- g ~ edges + transitiveties + nodematch("x")

  n_nodes <- 100
  n_blocks <- 2

  nodes_data <- tibble::tibble(
    node_id = 1:n_nodes,
    x = sample(1:2, size = n_nodes, replace = T),
    block = sample(1:n_blocks, size = n_nodes, replace = T)
  )

  g <- network::network.initialize(n = n_nodes)
  g %v% "block" <- nodes_data$block
  g %v% "x" <- nodes_data$x
  
  
  network::set.vertex.attribute(g, "x", nodes_data$x)
  list_feature_matrices <- get_features(g, bigergm_formula)

  coef_between_block <- c(-3, 1)
  coef_within_block <- c(-2, 0.1, 0.5)

  sim_ergm_control <- ergm::control.simulate.formula(
    MCMC.burnin = 4000000,
    MCMC.interval = 200000
  )

  g <- bigergm::simulate_bigergm(
    formula = bigergm_formula,
    coef_between = coef_between_block,
    coef_within = coef_within_block,
    control_within =  sim_ergm_control
  )

  bigergm_res <- bigergm::bigergm(
    g ~ edges + nodematch("x") + transitiveties,
    n_blocks = n_blocks,
    n_em_step_max = 10,
    estimate_parameters = T,
    clustering_with_features = T,
    method_within = 'MLE',
    control_within = ergm::control.ergm(
      MCMC.burnin = test_burnin,
      MCMC.interval = test_interval,
      main.method = test_method, 
      MCMLE.effectiveSize = NULL
    )
  )
  
  used_control <- bigergm_res$est_within$control

  expect_equal(used_control$MCMC.burnin, test_burnin)
  expect_equal(used_control$MCMC.interval, test_interval)
  expect_equal(used_control$main.method, test_method)
})
