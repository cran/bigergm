test_that("combining within- and between-block edgelists while removing duplicated edges works", {
  edgelist_within <- matrix(c(1, 1, 2, 3, 2, 3, 4, 5), nrow = 4, ncol = 2)
  edgelist_between <- matrix(c(1, 2, 3, 4, 3, 5), nrow = 3, ncol = 2)

  # Create a true edgelist
  true_edgelist <-
    data.frame(
      tail = c(1, 1, 2, 3, 1, 2, 3),
      head = c(2, 3, 4, 5, 4, 3, 5)
    ) %>%
    dplyr::distinct(tail, head) %>%
    dplyr::arrange(tail) %>%
    as.matrix()

    # Get edgelist from R function
  edgelist <- combine_within_between_edges(edgelist_within, edgelist_between, 
                                           use_fast_between_simulation = FALSE,
                                           old_vertex_names = c(1,2,3,4,5))

  # Check if it works
  expect_equal(true_edgelist, unclass(edgelist), check.attributes = FALSE)
})

test_that("correctly attaching vertex ids, block memberships, and vertex features to the simulated network works", {
  set.seed(1)
  # Prepare ingredients for simulating a network
  N <- 1000
  K <- 10

  list_within_params <- c(-2, 1, 1, 0.76, 0.08)
  list_between_params <- c(-7, 2, 2)
  formula <- g ~ edges + nodematch("x") + nodematch("y") + triangle + kstar(2)

  memb <- sample(1:K, size = N, replace = TRUE)
  vertex_id <- as.character(11:(11 + N - 1))

  x <- sample(1:20, size = N, replace = TRUE)
  y <- sample(LETTERS, size = N, replace = TRUE)


  df <- tibble::tibble(
    id = vertex_id,
    memb = memb,
    x = x,
    y = y
  )
  g <- network::network.initialize(n = N, directed = FALSE)
  g %v% "vertex.names" <- df$id
  g %v% "block" <- df$memb
  g %v% "x" <- df$x
  g %v% "y" <- df$y
  # Simulate a network
  g_sim <-
    simulate_bigergm(
      formula = formula,
      coef_within = list_within_params,
      coef_between = list_between_params,
      control_within = ergm::control.simulate.formula(),
      seed = 1,
      nsim = 1
    )
  
  
  # Check if vertex ids, block memberships, and vertex features are stored correctly in the network object.
  expect_equal(df$id, network::get.vertex.attribute(g_sim, "vertex.names"), check.attributes = FALSE)
  expect_equal(df$memb, network::get.vertex.attribute(g_sim, "block"), check.attributes = FALSE)
  expect_equal(df$x, network::get.vertex.attribute(g_sim, "x"), check.attributes = FALSE)
  expect_equal(df$y, network::get.vertex.attribute(g_sim, "y"), check.attributes = FALSE)
})

test_that("simulating a network from a given network works", {
  set.seed(1)
  # Prepare ingredients for simulating a network
  N <- 100
  K <- 3

  list_within_params <- c(-3, 1, 1, 0.76, 0.08)
  list_between_params <- c(-5, 2, 2)
  formula <- g ~ edges + nodematch("x") + nodematch("y") + triangle + kstar(2)

  memb <- sample(1:K, size = N, replace = TRUE)
  vertex_id <- as.character(11:(11 + N - 1))

  x <- sample(1:20, size = N, replace = TRUE)
  y <- sample(LETTERS, size = N, replace = TRUE)

  df <- tibble::tibble(
    id = vertex_id,
    memb = memb,
    x = x,
    y = y
  )
  g <- network::network.initialize(n = N, directed = FALSE)
  g %v% "vertex.names" <- df$id
  g %v% "block" <- df$memb
  g %v% "x" <- df$x
  g %v% "y" <- df$y
  
  
  # Simulate a network
  g_sim <-
    simulate_bigergm(
      formula = formula,
      coef_within = list_within_params,
      coef_between = list_between_params,
      control_within = ergm::control.simulate.formula(
        MCMC.burnin = 100,
        MCMC.interval = 10
      ),
      seed = 1,
      nsim = 1,
      only_within = TRUE
    )
  # Simulate a network from a given edgelist
  g2 <-
    simulate_bigergm(
      formula = formula,network = g_sim,
      coef_within = list_within_params,
      coef_between = list_between_params,
      control_within = ergm::control.simulate.formula(
        MCMC.burnin = 0,
        MCMC.interval = 1
      ),
      seed = 1,
      nsim = 1, 
      only_within = TRUE
    )
  
  expect_match(class(g2), "network")

  g2 <- network::as.edgelist(g2)
  g_sim <- network::as.edgelist(g_sim)
  # Check if the network is correctly generated
  expect_equal(nrow(g_sim), nrow(g2))
  expect_true(all(g_sim == g2))
})

test_that("simulating a network with a seed works", {
  set.seed(1)
  # Prepare ingredients for simulating a network
  N <- 100
  K <- 3
  
  list_within_params <- c(-3, 1, 1, 0.76, 0.08)
  list_between_params <- c(-5, 2, 2)
  formula <- g ~ edges + nodematch("x") + nodematch("y") + triangle + kstar(2)
  
  memb <- sample(1:K, size = N, replace = TRUE)
  vertex_id <- as.character(11:(11 + N - 1))
  
  x <- sample(1:20, size = N, replace = TRUE)
  y <- sample(LETTERS, size = N, replace = TRUE)
  
  df <- tibble::tibble(
    id = vertex_id,
    memb = memb,
    x = x,
    y = y
  )
  g <- network::network.initialize(n = N, directed = FALSE)
  g %v% "vertex.names" <- df$id
  g %v% "block" <- df$memb
  g %v% "x" <- df$x
  g %v% "y" <- df$y
  
  
  # Simulate a network only within-block
  g_sim <-
    simulate_bigergm(
      formula = formula,
      coef_within = list_within_params,
      coef_between = list_between_params,
      seed = 1,
      nsim = 1,
      only_within = TRUE
    )
  # Simulate a network from a given edgelist
  g_sim2 <-
    simulate_bigergm(
      formula = formula,
      coef_within = list_within_params,
      coef_between = list_between_params,
      seed = 1,
      nsim = 1, 
      only_within = TRUE
    )
  g_sim <- network::as.edgelist(g_sim)
  g_sim2 <- network::as.edgelist(g_sim2)
  expect_equal(g_sim, g_sim2)
  
  # Simulate whole network 
  g_sim <-
    simulate_bigergm(
      formula = formula,
      coef_within = list_within_params,
      coef_between = list_between_params,
      seed = 123,
      nsim = 1
    )
  # Simulate a network from a given edgelist
  g_sim2 <-
    simulate_bigergm(
      formula = formula,
      coef_within = list_within_params,
      coef_between = list_between_params,
      seed = 123,
      nsim = 1
    )
  g_sim <- network::as.edgelist(g_sim)
  g_sim2 <- network::as.edgelist(g_sim2)
  expect_equal(g_sim, g_sim2)

  
  # Simulate two networks and check the  
  g_sim <-
    simulate_bigergm(
      formula = formula,
      coef_within = list_within_params,
      coef_between = list_between_params,
      seed = 123,
      nsim = 2, 
      output = "stats"
    )
  # Simulate a network from a given edgelist
  g_sim2 <-
    simulate_bigergm(
      formula = formula,
      coef_within = list_within_params,
      coef_between = list_between_params,
      seed = 123,
      nsim = 2, 
      output = "stats"
    )
  expect_equal(g_sim$between_network, g_sim2$between_network)
  expect_equal(g_sim$within_network, g_sim2$within_network)
})

test_that("The simulation begins from the observed network given in the formula", {
  set.seed(1)
  # Prepare ingredients for simulating a network
  N <- 100
  K <- 10

  list_within_params <- c(-3, 1, 1, 0.76, 0.08)
  list_between_params <- c(-5, 2, 2)
  formula <- g ~ edges + nodematch("x") + nodematch("y") + triangle + kstar(2)

  memb <- sample(1:K, size = N, replace = TRUE)
  vertex_id <- as.character(11:(11 + N - 1))

  x <- sample(1:20, size = N, replace = TRUE)
  y <- sample(LETTERS, size = N, replace = TRUE)

  # Start with an epmty network 
  g <- network::network.initialize(n = N, directed = FALSE)
  g %v% "vertex.names" <- vertex_id
  g %v% "block" <-memb
  g %v% "x" <- x
  g %v% "y" <- y
  # Simulate a network
  g_sim <-
    simulate_bigergm(
      formula = formula,
      coef_within = list_within_params,
      coef_between = list_between_params,
      control_within = ergm::control.simulate.formula(
        MCMC.burnin = 0,
        MCMC.interval = 1
      ),
      seed = 1,
      nsim = 1, 
      only_within = TRUE
    )

  expect_match(class(g_sim), "network")

  g_sim <- network::as.edgelist(g_sim)

  # Check if the network is correctly generated
  expect_equal(nrow(g_sim), 0)
})

test_that("generating multiple networks works", {
  set.seed(1)
  # Prepare ingredients for simulating a network
  N <- 100
  K <- 10
  
  list_within_params <- c(-3, 1, 1, 0.76, 0.08)
  list_between_params <- c(-5, 2, 2)
  formula <- g ~ edges + nodematch("x") + nodematch("y") + triangle + kstar(2)
  
  memb <- sample(1:K, size = N, replace = TRUE)
  vertex_id <- as.character(11:(11 + N - 1))
  
  x <- sample(1:20, size = N, replace = TRUE)
  y <- sample(LETTERS, size = N, replace = TRUE)
  
  # Start with an epmty network 
  g <- network::network.initialize(n = N, directed = FALSE)
  g %v% "vertex.names" <- vertex_id
  g %v% "block" <-memb
  g %v% "x" <- x
  g %v% "y" <- y
  # Simulate a network
  g_sim <-
    simulate_bigergm(
      formula = formula,
      coef_within = list_within_params,
      coef_between = list_between_params,
      seed = 1,
      nsim = 3
    )
  
  expect_equal(length(g_sim), 3)
  expect_match(class(g_sim[[1]]), "network")
  expect_match(class(g_sim[[2]]), "network")
  expect_match(class(g_sim[[3]]), "network")
})

test_that("generating stats and avoiding generating within-block links while simulating a between-block network works", {
  set.seed(1)
  # Prepare ingredients for simulating a network
  N <- 100
  K <- 10
  
  list_within_params <- c(-3, 1, 1, 0.76, 0.08)
  list_between_params <- c(-5, 2, 2)
  formula <- g ~ edges + nodematch("x") + nodematch("y") + triangle + kstar(2)
  
  memb <- sample(1:K, size = N, replace = TRUE)
  vertex_id <- as.character(11:(11 + N - 1))
  
  x <- sample(1:20, size = N, replace = TRUE)
  y <- sample(LETTERS, size = N, replace = TRUE)
  
  # Start with an epmty network 
  g <- network::network.initialize(n = N, directed = FALSE)
  g %v% "vertex.names" <- vertex_id
  g %v% "block" <-memb
  g %v% "x" <- x
  g %v% "y" <- y

  # Simulate a network
  g_sim <-
    simulate_bigergm(
      formula = formula,
      coef_within = list_within_params,
      coef_between = list_between_params,
      seed = 1,
      nsim = 3, output = "stats"
    )

  # Check whether network stats are generated.
  expect_equal(c("edges", "nodematch.x", "nodematch.y", "triangle", "kstar2"), colnames(g_sim[[1]]))
  expect_equal(c("edges", "nodematch.x", "nodematch.y"), colnames(g_sim[[2]]))
})

test_that("Simulating networks using formula without externality terms works", {
  set.seed(1)
  # Prepare ingredients for simulating a network
  N <- 100
  K <- 10
  
  list_within_params <- c(-3,0.76, 0.08)
  list_between_params <- c(-5)
  formula <- g ~ edges + triangle + kstar(2)
  
  memb <- sample(1:K, size = N, replace = TRUE)
  vertex_id <- as.character(11:(11 + N - 1))

  # Start with an epmty network 
  g <- network::network.initialize(n = N, directed = FALSE)
  g %v% "vertex.names" <- vertex_id
  g %v% "block" <-memb
  
  expect_error( simulate_bigergm(
    formula = formula,
    coef_within = list_within_params,
    coef_between = list_between_params,
    seed= 1,
    nsim = 3
  ), NA)
})

test_that("Simulating networks using formula without covariates works", {
  set.seed(1)
  # Prepare ingredients for simulating a network
  N <- 100
  K <- 10
  
  list_within_params <- c(-3)
  list_between_params <- c(-5)
  formula <- g ~ edges 
  
  memb <- sample(1:K, size = N, replace = TRUE)
  vertex_id <- as.character(11:(11 + N - 1))
  
  # Start with an epmty network 
  g <- network::network.initialize(n = N, directed = FALSE)
  g %v% "vertex.names" <- vertex_id
  g %v% "block" <-memb
  
  expect_error( simulate_bigergm(
    formula = formula,
    coef_within = list_within_params,
    coef_between = list_between_params,
    seed= 1,
    nsim = 3
  ), NA)
})

test_that("Setting output = 'edgelist' returns a single edgelist for the full
          simulated network", {
  net <- bigergm::toyNet
  # Shuffle the vertex names because we want to check whether the IDs are in the correct order
  network::network.vertex.names(net) <- as.character(sample(401:600, size = 200, replace = FALSE))
  nodes_data <- intergraph::asDF(net)$vertexes %>%
    dplyr::select(vertex.names, x, y, block) %>%
    dplyr::rename(id = vertex.names)

  model_formula <-
    net ~ edges + nodematch('x') + nodematch('y') + triangle + kstar(2)

  
  bigergm_res <- bigergm::bigergm(
    object = model_formula,
    n_blocks = 4,n_MM_step_max = 10,
    clustering_with_features = TRUE,
    verbose = 1
  )

  simulate_output <- function(output){
    bigergm::simulate_bigergm(
      formula = model_formula,
      seed = 44,
      coef_within = bigergm_res$est_within$coefficients,
      coef_between = bigergm_res$est_between$coefficients,
      control_within = ergm::control.simulate.formula(
        MCMC.burnin = 100000,
        MCMC.interval = 500
      ),
      nsim = 1,
      output = output
    )
  }

  sim_net <- simulate_output('network')
  # The class of the returned object should be `network`
  expect_s3_class(sim_net, 'network')
  # Convert the network to an edgelist for comparison
  sim_net_edgelist <- network::as.edgelist(sim_net)
  sim_edgelist <- simulate_output('edgelist')
  # The class of the returned object should be `edgelist`
  expect_s3_class(sim_edgelist, 'edgelist')

  # The list of edges should be the same
  expect_true(all(sim_net_edgelist, sim_edgelist))

  # Both should have the same number of nodes
  expect_true(attr(sim_net_edgelist, 'n') == attr(sim_edgelist, 'n'))

  expect_true(all(attr(sim_net_edgelist, 'vnames') == attr(sim_edgelist, 'vnames')))
})

test_that("Setting output = 'edgelist' returns a list of full-network edgelists", {
  net <- bigergm::toyNet
  # Shuffle the vertex names because we want to check whether the IDs are in the correct order
  network::network.vertex.names(net) <- as.character(sample(401:600, size = 200, replace = FALSE))
  nodes_data <- intergraph::asDF(net)$vertexes %>%
    dplyr::select(vertex.names, x, y, block) %>%
    dplyr::rename(id = vertex.names)

  model_formula <-
    net ~ edges + nodematch('x') + nodematch('y') + triangle + kstar(2)

  bigergm_res <- bigergm::bigergm(
    object = model_formula,
    n_blocks = 4,n_MM_step_max = 10,
    clustering_with_features = TRUE,
    verbose = 1
  )

  simulate_output <- function(output){
    bigergm::simulate_bigergm(
      formula = model_formula,
      seed = 44,
      coef_within = bigergm_res$est_within$coefficients,
      coef_between = bigergm_res$est_between$coefficients,
      control_within = ergm::control.simulate.formula(
        MCMC.burnin = 100000,
        MCMC.interval = 500
      ),
      nsim = 3,
      output = output
    )
  }

  sim_nets <- simulate_output('network')
  # The returned object should be a list of networks
  expect_equal(class(sim_nets), 'list')

  for(i in 1:length(sim_nets)){
    expect_s3_class(sim_nets[[i]], 'network')
  }

  # Convert the network to an edgelist for comparison
  sim_nets_edgelists <- sim_nets %>%
    purrr::map(network::as.edgelist)

  
  sim_edgelists <- simulate_output('edgelist')
  # The class of the returned object should be `list`
  expect_equal(class(sim_edgelists), 'list')

  for(i in 1:length(sim_edgelists)){
    expect_s3_class(sim_edgelists[[i]], 'edgelist')
  }

  expect_equal(length(sim_nets_edgelists), 3)
  expect_equal(length(sim_edgelists), 3)

  # The list of edges should be the same
  for (i in 1:3){
    expect_true(all(sim_nets_edgelists[[i]], sim_edgelists[[i]]))
    # Both should have the same number of nodes
    expect_true(attr(sim_nets_edgelists[[i]], 'n') == attr(sim_edgelists[[i]], 'n'))
    expect_true(all(attr(sim_nets_edgelists[[i]], 'vnames') == attr(sim_edgelists[[i]], 'vnames')))
  }
})

