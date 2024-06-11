test_that("starting EM iterations and parameter estimation from a given vector of block memberships works", {
  skip_on_cran()
  set.seed(334)
  # Simulate a network to work with in this unit test.
  # Number of nodes
  N <- 1000
  # Number of blocks
  K <- 20
  # Block memberships (same block size)
  memb <- rep(1:K, each = N / K)
  # Covariates
  x <- sample(1:10, size = N, replace = TRUE)
  y <- sample(1:10, size = N, replace = TRUE)

  # Within-block parameters: edges, nodematch("x"), nodematch("y"), triangle
  list_within_params <- c(-3, 1, 1, 0.25)
  # Between-block parameters: edges, nodematch("x"), nodematch("y")
  list_between_params <- c(-8.5, 0.5, 0.5)


  vertex_id <- 1:N

  df <- tibble::tibble(
    id = vertex_id,
    memb = memb,
    x = x,
    y = y
  )
  set.seed(334)
  g <- network::network(x = N,  directed = FALSE)
  g %v% "block" <- df$memb
  g %v% "x" <- df$x
  g %v% "y" <- df$y
  formula <- g ~ edges + nodematch("x") + nodematch("y") + triangle
  set.seed(334)
  g_sim <-
    simulate_bigergm(
      formula = formula,
      coef_within = list_within_params,
      coef_between = list_between_params,
      control_within = ergm::control.simulate.formula(MCMC.burnin = 1000000, MCMC.interval = 1000),
      seed = 1,
      nsim = 1,
      output = "network"
    )
  # Conduct clustering
  cluster_with_feature <-
    bigergm::bigergm(g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,
      n_blocks = K,
      estimate_parameters = FALSE,
      verbose = 0,
      n_MM_step_max =  3,
      initializatiom = 3,
      infomap_python = FALSE,
      clustering_with_features = TRUE,
      check_alpha_update = TRUE,
      compute_pi = TRUE,
      check_lower_bound = TRUE,
      check_blocks = TRUE,seed = 123
    )

  # Check if starting from the previously estimated block memberships works.
  expect_error(result <-
                 bigergm::bigergm(g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,
                                  initialization = cluster_with_feature$block,
      n_MM_step_max = 2,
      estimate_parameters = FALSE,
      seed = 123
    ), NA)

  # Check if starting from block memberships initialized Python's infomap works.
  expect_error(result2 <-
    bigergm::bigergm(g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,
                     initialization =  system.file("extdata", "initialized_cluster_data_by_infomap.clu", package = "bigergm"),
      n_MM_step_max= 1,
      estimate_parameters = FALSE,
      verbose = 1,seed = 123,
    ), NA)

  
  # Check if starting parameter estimation from a given vector of block memberships works.
  expect_error(result3 <-
    bigergm::bigergm(g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,blocks = result$block,
      verbose = 1,
      seed = 123
    ), NA)

  # Check if not specifying n_blocks when initialized_cluster_data and block_membership are null yields an error.
  expect_error(result4 <-
    bigergm::bigergm(g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,
      verbose = 1,
      seed = 123
    ))

  })
