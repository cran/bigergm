test_that("Setting a checkpoint for MM iterations works", {
  skip_on_cran()

  set.seed(334)
  # Simulate a network to work with in this unit test.
  # Number of nodes
  N <- 50
  # Number of blocks
  K <- 2
  # Block memberships (same block size)
  memb <- rep(1:K, each = N / K)
  # Covariates
  x <- sample(1:10, size = N, replace = TRUE)
  y <- sample(1:10, size = N, replace = TRUE)

  # Within-block parameters: edges, nodematch("x"), nodematch("y")
  list_within_params <- c(-1, 1, 1)
  # Between-block parameters: edges, nodematch("x"), nodematch("y")
  list_between_params <- c(-3.5, 0.5, 0.5)

  formula <- g ~ edges + nodematch("x") + nodematch("y") 
  g <- network::network.initialize(N, directed = FALSE)
  network::set.vertex.attribute(g, "x", x)
  network::set.vertex.attribute(g, "y", y)
  network::set.vertex.attribute(g, "block", memb)


  g_sim <-
    simulate_bigergm(
      formula = formula,
      coef_within = list_within_params,
      coef_between = list_between_params,
      control_within = ergm::control.simulate.formula(MCMC.burnin = 1000, MCMC.interval = 10),
      seed = 1,
      nsim = 1,
      output = "network"
    )

  ############# 1. Clustering with features ##############################
  # Conduct clustering at once
  initial_weight <- 1000

  cluster_with_feature <-
    bigergm::bigergm(g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,
      n_blocks = K,
      estimate_parameters = TRUE,
      verbose = 0,tol_MM_step = 0.000000001,
      n_MM_step_max = 10,
      initialization = "spectral",
      infomap_python = FALSE,
      clustering_with_features = TRUE,
      check_alpha_update = TRUE,
      compute_pi = TRUE,
      check_lower_bound = TRUE,
      check_block_membership = TRUE,
      weight_for_initialization = initial_weight,
      seed = 334
    )

  
  # Conduct clustering in two steps
  first_step <-
    bigergm::bigergm(
      object = g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,
      n_blocks = K,tol_MM_step = 0.000000001,
      estimate_parameters = TRUE,
      verbose = 0,
      n_MM_step_max = 7,
      initialization = "spectral",
      infomap_python = FALSE,
      clustering_with_features = TRUE,
      check_block_membership = TRUE,
      weight_for_initialization = initial_weight,
      check_alpha_update = TRUE,
      seed = 334
    )

  second_step <-
    bigergm::bigergm(
      object = first_step,tol_MM_step = 0.000000001,
      n_MM_step_max = 3
    )

  
  
  # Check if the calculated lower bounds are identical (both length and values)
  expect_equal(cluster_with_feature$MM_lower_bound, second_step$MM_lower_bound)

  # Check if block memberships are identical over iterations (both length and Yule's coefficient)
  expect_true(all(unlist(purrr::map2(cluster_with_feature$MM_list_z, second_step$MM_list_z, yule) == rep(1, length(second_step$MM_list_z)))))

  # The block should be the same at the end of the second estimation as the one obtained after running 10 iterations of the MM algorithm
  expect_equal(yule(cluster_with_feature$block, second_step$block), 1)

  # Check if alphas are identical over iterations
  for (i in 1:length(cluster_with_feature$MM_list_alpha)) {
    expect_equal(cluster_with_feature$MM_list_alpha[[i]], second_step$MM_list_alpha[[i]], check.attribute = FALSE, tolerance = 1e-2)
  }
  # Check if the alpha after the second checkpoint is the same as the alpha when performing estimation with 10 MM iterations without a checkpoint.
  expect_equal(cluster_with_feature$alpha, second_step$alpha)

  # Check if estimated coefficients with and without a checkpoint are identical.
  expect_equal(coef(cluster_with_feature$est_between), coef(second_step$est_between), tolerance = 1e-10)
  expect_equal(coef(cluster_with_feature$est_within), coef(second_step$est_within), tolerance = 1e-10)


  ############# 2. Clustering without features ##############################
  # Conduct clustering at once
  initial_weight <- 1000

  cluster_without_feature <-
    bigergm::bigergm(g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,
      n_blocks = K,
      estimate_parameters = TRUE,
      verbose = 0,
      n_MM_step_max = 5,
      initialization = "spectral",
      infomap_python = FALSE,
      clustering_with_features = FALSE,
      check_alpha_update = TRUE,
      compute_pi = TRUE,tol_MM_step = 0.000000001,
      check_lower_bound = TRUE,
      check_block_membership = TRUE,
      weight_for_initialization = initial_weight,
      seed = 334
    )

  # Conduct clustering in two steps
  first_step_without_feature <-
    bigergm::bigergm(g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,
      n_blocks = K,
      estimate_parameters = TRUE,
      verbose = 0,tol_MM_step = 0.000000001,
      n_MM_step_max = 3,
      initialization = "spectral",
      infomap_python = FALSE,
      clustering_with_features = FALSE,
      check_alpha_update = TRUE,
      compute_pi = TRUE,
      check_lower_bound = TRUE,
      check_block_membership = TRUE,
      weight_for_initialization = initial_weight,
      seed = 334
    )
  second_step_without_feature <-
    bigergm::bigergm(
      object = first_step_without_feature,
      n_MM_step_max = 2,tol_MM_step = 0.000000001
    )

  # Check if the calculated lower bounds are identical (both length and values)
  expect_equal(cluster_without_feature$EM_lower_bound, second_step_without_feature$EM_lower_bound)

  # Check if block memberships are identical over iterations (both length and Yule's coefficient)
  expect_true(all(unlist(purrr::map2(cluster_without_feature$EM_list_z, second_step_without_feature$EM_list_z, yule) ==
    rep(1, length(second_step_without_feature$EM_list_z)))))

  # The block should be the same at the end of the second estimation as the one obtained after running 10 iterations of the MM algorithm
  expect_equal(yule(cluster_without_feature$block, second_step_without_feature$block), 1)

  # Check if alphas are identical over iterations
  for (i in 1:6) {
    expect_equal(cluster_without_feature$MM_list_alpha[[i]], second_step_without_feature$MM_list_alpha[[i]], check.attribute = FALSE, tolerance = 1e-2)
  }
  # Check if the alpha after the second checkpoint is the same as the alpha when performing estimation with 10 MM iterations without a checkpoint.
  expect_equal(cluster_without_feature$alpha, second_step_without_feature$alpha)

  # Check if estimated coefficients with and without a checkpoint are identical.
  expect_equal(coef(cluster_without_feature$est_between), coef(second_step_without_feature$est_between), tolerance = 1e-10)
  expect_equal(coef(cluster_without_feature$est_within), coef(second_step_without_feature$est_within), tolerance = 1e-10)
})
