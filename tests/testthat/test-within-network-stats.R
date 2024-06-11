test_that("generating multiple within-block networks works", {
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
  
  # Obtain the stats
  within_sim_stats <- simulate_bigergm(
    formula = formula,
    control_within = ergm::control.simulate.formula(),
    seed = 1,
    nsim = 3, only_within = TRUE, output = "stats",
    coef_between = list_between_params,
    coef_within = list_within_params
  )

  expected_terms <- statnet.common::list_rhs.formula(formula)

  expect_equal(nrow(within_sim_stats$within_network), 3)
  expect_equal(length(expected_terms), length(names(within_sim_stats$within_network)))

})

test_that("simulating a network from a given edgelist works", {
  set.seed(1)
  # Prepare ingredients for simulating a network
  N <- 1000
  K <- 10

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
      control_within = ergm::control.simulate.formula(),
      seed = 1,
      nsim = 1, only_within = TRUE, output = "network",
      coef_between = list_between_params,
      coef_within = list_within_params
    )

  
  formula <- g_sim ~ edges + nodematch("x") + nodematch("y") + triangle + kstar(2)
  
  # Simulate a within-block network from a given edgelist
  g2 <-  simulate_bigergm(
    formula = formula,
    control_within =  ergm::control.simulate.formula(
      MCMC.burnin = 0,
      MCMC.interval = 1
    ),
    seed = 1,
    nsim = 1, only_within = TRUE, output = "network",
    coef_between = list_between_params,
    coef_within = list_within_params
  ) 
  

  expect_match(class(g2), "network")

  g2 <- network::as.edgelist(g2)
  g_sim <- network::as.edgelist(g_sim)
  
  # Check if the network is correctly generated
  expect_equal(nrow(g_sim), nrow(g2))
  expect_true(all(g_sim == g2))
})

test_that("The within-simulation begins from an empty network by default", {
  set.seed(1)
  # Prepare ingredients for simulating a network
  N <- 1000
  K <- 10
  
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
      control_within =  ergm::control.simulate.formula(
        MCMC.burnin = 0,
        MCMC.interval = 1
      ),
      seed = 1,
      nsim = 1, only_within = TRUE, output = "network",
      coef_between = list_between_params,
      coef_within = list_within_params
    )
  
  expect_match(class(g_sim), "network")
  g_sim <- network::as.edgelist(g_sim)

  # Check if the network is correctly generated
  expect_equal(nrow(g_sim), 0)
})
