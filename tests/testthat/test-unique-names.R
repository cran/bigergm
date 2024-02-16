test_that("Check if algorithm works even if vertex.names are nonunique", {
  # Simulate an ERGM with local dependence and two clusters and sampled covariates x and y
  set.seed(123) # Set seed
  directed <- F
  n_actors <- 300
  empty_network <- network(directed = directed, x = n_actors)
  x_info <- sample(x = c(1,2,3),size = n_actors, replace = T)
  y_info <- sample(x = c(1,2,3),size = n_actors, replace = T)
  empty_network %v% "y" = y_info
  empty_network %v% "x" = x_info
  # The subnetworks are sampled separately
  network_11 <- simulate(empty_network~edges +nodematch("x") + nodematch("y")+ gwesp(decay = log(2), fixed = T), coef = c(-3.5,0.5,-0.5,0.25))
  network_22 <- simulate(empty_network~edges +nodematch("x") + nodematch("y") + gwesp(decay = log(2), fixed = T), coef = c(-3.5,0.5,-0.5,0.25))
  network_12 <- simulate(empty_network~edges +nodematch("x") + nodematch("y"), coef = c(-7,0.5,0.5))
  # And then stitched together
  network = rbind(cbind(as.matrix(network_11),as.matrix(network_12)), 
                  cbind(t(as.matrix(network_12)),as.matrix(network_22)))
  
  # Start with the version where there are noninuque vertex.names 
  tmp  <- network(directed = directed, network)
  tmp %v% "y" <- c(y_info, y_info)
  tmp %v% "x" <- c(x_info, x_info)
  model_formula <- tmp~edges +nodematch("x") + nodematch("y")+ gwesp(decay = log(2), fixed = T)
  hergm_res_nonunique <- expect_warning(hergm(verbose = F,object = model_formula, # The model you would like to estiamte
                       n_clusters = 2, # The number of blocks
                       n_MM_step_max =50, # The maximum number of EM algorithm steps
                       estimate_parameters = T, # Perform parameter estimation after the block recovery step
                       clustering_with_features = T, initialization_method =1,
                       check_block_membership = TRUE))
  
  # Now set the vertex.names before to be unique
  rownames(network) = 1:nrow(network)
  colnames(network) = 1:nrow(network)
  tmp  <- network(directed = directed, network)
  tmp %v% "y" <- c(y_info, y_info)
  tmp %v% "x" <- c(x_info, x_info)
  model_formula <- tmp~edges +nodematch("x") + nodematch("y")+ gwesp(decay = log(2), fixed = T)
  hergm_res_unique <- expect_no_warning(hergm(verbose = F,object = model_formula, # The model you would like to estiamte
                                              n_clusters = 2, # The number of blocks
                                              n_MM_step_max =50, # The maximum number of EM algorithm steps
                                              estimate_parameters = T, # Perform parameter estimation after the block recovery step
                                              clustering_with_features = T, initialization_method =1,
                                              check_block_membership = TRUE))

  expect_equal(hergm_res_nonunique$partition, hergm_res_unique$partition)
  

  })