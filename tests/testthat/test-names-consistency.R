test_that("Check if algorithm works even if vertex.names are nonunique", {
  # Simulate an ERGM with local dependence and two clusters and sampled covariates x and y
  set.seed(123) # Set seed
  directed <- TRUE
  n_actors <- 300
  empty_network <- network(directed = directed, x = n_actors)
  x_info <- sample(x = c(1,2,3),size = n_actors, replace = TRUE)
  y_info <- sample(x = c(1,2,3),size = n_actors, replace = TRUE)
  empty_network %v% "y" = y_info
  empty_network %v% "x" = x_info
  # The subnetworks are sampled separately
  network_11 <- simulate(empty_network~edges +nodematch("x") + nodematch("y")+ gwesp(decay = log(2), fixed = TRUE), coef = c(-3.5,0.5,-0.5,0.25))
  network_22 <- simulate(empty_network~edges +nodematch("x") + nodematch("y") + gwesp(decay = log(2), fixed = TRUE), coef = c(-3.5,0.5,-0.5,0.25))
  network_12 <- simulate(empty_network~edges +nodematch("x") + nodematch("y"), coef = c(-7,0.5,0.5))
  network_21 <- simulate(empty_network~edges +nodematch("x") + nodematch("y"), coef = c(-7,0.5,0.5))

  # And then stitched together
  network = rbind(cbind(as.matrix(network_11),as.matrix(network_12)),
                  cbind(t(as.matrix(network_21)),as.matrix(network_22)))

  # Start with the version where there are noninuque vertex.names
  tmp  <- network(directed = directed, network)
  tmp %v% "y" <- c(y_info, y_info)
  tmp %v% "x" <- c(x_info, x_info)
  formula <- tmp~ edges + nodematch("x") + nodematch("y")+ gwesp(decay = log(2), fixed = T)
  
  bigergm_res_repeated <- bigergm(verbose = F,object = formula, # The model you would like to estimate
                                   n_blocks = 2, # The number of blocks
                                   n_MM_step_max =50, # The maximum number of EM algorithm steps
                                   estimate_parameters =F, # Perform parameter estimation after the block recovery step
                                   clustering_with_features = T, initialization ="infomap",
                                   check_blocks = TRUE, add_intercepts = F)
  
  expect_equal(tmp%v%"vertex.names", 
               bigergm_res_repeated$checkpoint$network %v% "vertex.names")
  # Now set the vertex.names before to be unique
  rownames(network) = 1:nrow(network)
  colnames(network) = 1:nrow(network)
  tmp  <- network(directed = directed, network)
  tmp %v% "y" <- c(y_info, y_info)
  tmp %v% "x" <- c(x_info, x_info)
  model_formula <- tmp~edges +nodematch("x") + nodematch("y")+ gwesp(decay = log(2), fixed = T)
  tmp%v%"vertex.names" <- letters[sample(1:26, 2*n_actors, replace = T)]
  
  bigergm_res_letter <- expect_no_warning(bigergm(verbose = F,object = model_formula,
                                              n_blocks = 2, # The number of blocks
                                              n_MM_step_max =50, # The maximum number of EM algorithm steps
                                              estimate_parameters = F, # Perform parameter estimation after the block recovery step
                                              clustering_with_features = T, initialization ="infomap",
                                              check_blocks = TRUE))
  expect_equal(tmp%v%"vertex.names", 
               bigergm_res_letter$checkpoint$network %v% "vertex.names")
  # Check if the blocks are the same
  expect_equal(bigergm_res_letter$block, bigergm_res_repeated$block)
  })