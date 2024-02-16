test_that("Check if algorithm automatically checks if covariates are specified", {
  # Simulate an ERGM with local dependence and two clusters and sampled covariates x and y
  set.seed(123) # Set seed
  directed <- F
  n_actors <- 25
  empty_network <- network(directed = directed, x = n_actors)
  # The subnetworks are sampled separately
  network_11 <- simulate(empty_network~edges + gwesp(decay = log(2), fixed = T), coef = c(-2.5,0.25))
  network_22 <- simulate(empty_network~edges  + gwesp(decay = log(2), fixed = T), coef = c(-2.5,0.25))
  network_12 <- simulate(empty_network~edges, coef = c(-5))
  # And then stitched together
  network = rbind(cbind(as.matrix(network_11),as.matrix(network_12)), 
                  cbind(t(as.matrix(network_12)),as.matrix(network_22)))
  # Now set the vertex.names before to be unique
  rownames(network) = 1:nrow(network)
  colnames(network) = 1:nrow(network)
  # Set up network object
  tmp  <- network(directed = directed, network)
  model_formula <- tmp~edges+ gwesp(decay = log(2), fixed = T)
  
  hergm_res_covariate <- expect_warning(hergm(verbose = F,object = model_formula, # The model you would like to estimate
                       n_clusters = 2, # The number of blocks
                       n_MM_step_max =2, # The maximum number of EM algorithm steps
                       estimate_parameters = TRUE, # Perform parameter estimation after the block recovery step
                       clustering_with_features = TRUE, initialization_method =1,
                       check_block_membership = TRUE, seed_infomap = 123))
  hergm_res_no_covariate <- expect_no_warning(hergm(verbose = F,object = model_formula, # The model you would like to estimate
                                                    n_clusters = 2, # The number of blocks
                                                    n_MM_step_max =2, # The maximum number of EM algorithm steps
                                                    estimate_parameters = TRUE, # Perform parameter estimation after the block recovery step
                                                    clustering_with_features = FALSE, initialization_method =1,
                                                    check_block_membership = TRUE, seed_infomap = 123))

  expect_equal(hergm_res_covariate$est_within$coefficients, hergm_res_no_covariate$est_within$coefficients)
  
  
  
  })