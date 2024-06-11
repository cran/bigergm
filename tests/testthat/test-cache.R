set.seed(334)

# Simulate a random network for testing
simulate_network <- function(){
  N <- 50
  K <- 5
  memb <- rep(1:K, each = N / K)
  x <- sample(1:10, size = N, replace = TRUE)
  y <- sample(1:10, size = N, replace = TRUE)
  list_within_params <- c(-1, 1, 1, 0.5)
  list_between_params <- c(-3.5, 0.5, 0.5)

  formula <- g ~ edges + nodematch("x") + nodematch("y") + triangle

  vertex_id <- 1:N

  df <- tibble::tibble(
    id = vertex_id,
    memb = memb,
    x = x,
    y = y
  )
  g <- network::network.initialize(n = N, directed = FALSE)
  network::set.vertex.attribute(g, "block", df$memb)
  network::set.vertex.attribute(g, "x", df$x)
  network::set.vertex.attribute(g, "y", df$y)
  
  simulate_bigergm(
      formula = formula,
      coef_within = list_within_params,
      coef_between = list_between_params,
      control_within = ergm::control.simulate.formula(MCMC.burnin = 1000000, MCMC.interval = 1000),
      seed = 1,
      nsim = 1,
      output = "network"
    )
}

# Function that checks that the number of files in the directory is the expected number
check_files <- function(dir, expected_number){
  expect_equal(length(list.files(dir, '*.rds')), expected_number)
}

cleanup <- function(){
  do.call(file.remove, list(list.files(tempdir(), '*.rds', full.names = TRUE)))
}

test_that('Estimation with a disk cache stores data in the correct directory', {
  skip_on_cran()
  on.exit(cleanup())

  # Use this directory for caching in disk
  dir <- tempdir()

  # Simulate a network
  g_1 <- simulate_network()

  # There should be no cached files in the directory
  check_files(dir, 0)

  # Perform estimation
  bigergm::bigergm(
    object = g_1 ~ edges + nodematch("x"),
    n_blocks =  3,
    n_MM_step_max = 3,
    initialization = "infomap",
    clustering_with_features = TRUE,
    verbose=2,
    cache = cachem::cache_disk(dir)
  )

  # The estimation should have stored one RDS object in the cache directory.
  check_files(dir, 1)

  bigergm::bigergm(
    object = g_1 ~ edges + nodematch("x"),
    n_blocks = 3,n_MM_step_max = 3,
    initialization = "infomap",
    clustering_with_features = TRUE,
    verbose=2,
    cache = cachem::cache_disk(dir)
  )

  # Running again the estimation on the same network should reuse the previously stored RDS object
  # and not store a new one.
  check_files(dir, 1)

  # Generate a different network
  g_2 <- simulate_network()

  # Perform estimation on the new network.
  bigergm::bigergm(
    object = g_2 ~ edges + nodematch("x"),
    n_blocks = 3,n_MM_step_max = 3,
    initialization = "infomap",
    clustering_with_features = TRUE,
    verbose=2,
    cache = cachem::cache_disk(dir)
  )

  # The network has changed, so the previously cached RDS is not reused and a new cache file is generated.
  check_files(dir, 2)

})
