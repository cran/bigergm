test_that("Check if algorithm automatically checks if covariates are specified", {
  # Simulate an ERGM with local dependence and two clusters and sampled covariates x and y
  set.seed(123) # Set seed
  directed <- F
  n_actors <- 25
  empty_network <- network.initialize(directed = directed, n = n_actors)
  # The subnetworks are sampled separately
  network_11 <- simulate_formula(empty_network~edges + gwesp(decay = log(2), fixed = T), coef = c(-2.5,0.25))
  network_22 <- simulate_formula(empty_network~edges  + gwesp(decay = log(2), fixed = T), coef = c(-2.5,0.25))
  network_12 <- simulate_formula(empty_network~edges, coef = c(-5))
  # And then stitched together
  network = rbind(cbind(as.matrix(network_11),as.matrix(network_12)), 
                  cbind(t(as.matrix(network_12)),as.matrix(network_22)))
  # Now set the vertex.names before to be unique
  rownames(network) = 1:nrow(network)
  colnames(network) = 1:nrow(network)
  # Set up network object
  tmp  <- network(directed = directed, network)
  model_formula <- tmp~edges+ gwesp(decay = log(2), fixed = T)
  
  bigergm_res_covariate <- expect_warning(bigergm(verbose = F,object = model_formula, # The model you would like to estimate
                       n_blocks = 2, # The number of blocks
                       n_MM_step_max =2, # The maximum number of MM algorithm steps
                       estimate_parameters = TRUE, # Perform parameter estimation after the block recovery step
                       clustering_with_features = TRUE, initialization ="infomap",
                       check_blocks = TRUE, seed_infomap = 123))
  
  bigergm_res_no_covariate <- expect_no_warning(bigergm(verbose = F,object = model_formula, # The model you would like to estimate
                                                    n_blocks = 2, # The number of blocks
                                                    n_MM_step_max =2, # The maximum number of EM algorithm steps
                                                    estimate_parameters = TRUE, # Perform parameter estimation after the block recovery step
                                                    clustering_with_features = FALSE, initialization ="infomap",
                                                    check_blocks = TRUE, seed_infomap = 123))

  expect_equal(bigergm_res_covariate$est_within$coefficients, bigergm_res_no_covariate$est_within$coefficients)
  
  
  
  })
