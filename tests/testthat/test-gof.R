get_dummy_net <- function(n_nodes, n_blocks,  seed = 1234) {
  bigergm_formula <- g ~ edges  + nodematch("x")

  nodes_data <- tibble::tibble(
    node_id = 1:n_nodes,
    x = sample(1:2, size = n_nodes, replace = T),
    block = sample(1:n_blocks, size = n_nodes, replace = T)
  )

  g <- network::network.initialize(n = n_nodes,directed = FALSE)
  g%v% "block" <- nodes_data$block
  network::set.vertex.attribute(g, "x", nodes_data$x)
  preprocessed_features <- get_features(g, bigergm_formula)

  coef_between_block <- c(-3, 1)
  coef_within_block <- c(-2, 0.1)

  sim_control_within <- ergm::control.simulate.formula(
    MCMC.burnin = 400,
    MCMC.interval = 200
  )

  g <- bigergm::simulate_bigergm(
    formula = bigergm_formula,
    coef_between = coef_between_block,
    coef_within = coef_within_block,
    control_within = sim_control_within
  )

  
  bigergm_res <- bigergm::bigergm(
    g ~ edges + nodematch("x"),
    n_blocks = n_blocks,
    n_MM_step_max = 2,
    estimate_parameters = T,
    use_infomap_python = F, 
    clustering_with_features = T
  )

  list(
    bigergm_res = bigergm_res,
    g = g,
    nodes_data = nodes_data,
    K = n_blocks,
    preprocessed_features = preprocessed_features,
    vertex_id_var = "node_id",
    block_id_var = "block",
    control_within = sim_control_within
  )
}
sim <- get_dummy_net(n_nodes = 50,n_blocks =  2, seed = 1234)


test_that("Returned GOF dataframe has the correct fields", {
  g <- sim$g
  test_gof_res <- gof(
    sim$bigergm_res,
    control_within = sim$control_within,
    nsim = 3, 
    compute_geodesic_distance = FALSE
  )
  
  
  for (stat_type in c("original", "simulated")) {
    stats <- test_gof_res[[stat_type]]
    expect_false(is.null(stats))
    for (stat in c("network_stats", "degree_dist", "esp_dist")) {
      expect_false(is.null(stats[[stat]]))
    }
    expect_true(is.null(stats[["geodesic_dist"]]))
  }
})

test_that("GOF network stats have the right fields and terms", {
  g <- sim$g

  test_gof_res <- gof(
    sim$bigergm_res,
    control_within = sim$control_within,
    nsim = 3
  )

  expected_terms <- ergm::ergm_model(sim$bigergm_res$est_within$formula)$terms %>%
    purrr::map(function(t) {
      `$`(t, name)
    })

  for (stat_type in c("original", "simulated")) {
    stat_type_df <- test_gof_res[[stat_type]]
    actual_terms <- colnames(stat_type_df$network_stats)

    actual_terms[stringr::str_detect(actual_terms, "nsim", negate = TRUE)] %>%
      setdiff(c("value", "stat")) %>%
      length() %>%
      expect_equal(0)

    stat_type_df$network_stats$stat %>%
      unique() %>%
      stringr::str_replace("[.].*", "") %>%
      setdiff(expected_terms) %>%
      length() %>%
      expect_equal(0)
  }

})

test_that("GOF degree stats have the right fields and terms", {
  g <- sim$g

  test_gof_res <- gof(
    sim$bigergm_res,
    control_within = sim$control_within,
    nsim = 3
  )

  for (stat_type in c("original", "simulated")) {
    stat_type_df <- test_gof_res[[stat_type]]
    actual_terms <- colnames(stat_type_df$degree_dist)

    actual_terms[stringr::str_detect(actual_terms, "nsim", negate = TRUE)] %>%
      setdiff(c("degree", "share")) %>%
      length() %>%
      expect_equal(0)

    expect_lte(max(stat_type_df$degree_dist$degree), g$gal$n)
    expect(
      min(stat_type_df$degree_dist$share) >= 0 && max(stat_type_df$degree_dist$share) <= 1,
      failure_message = "Some degree shares are out of bounds"
    )
  }
})

test_that("GOF esp stats have the right fields and terms", {
  g <- sim$g
  
  test_gof_res <-  gof(
    sim$bigergm_res,
    control_within = sim$control_within,
    nsim = 3
  )

  
  for (stat_type in c("original", "simulated")) {
    stat_type_df <- test_gof_res[[stat_type]]
    actual_terms <- colnames(stat_type_df$esp_dist)
    actual_terms[stringr::str_detect(actual_terms, "nsim", negate = TRUE)] %>%
      setdiff(c("label", "esp")) %>%
      length() %>%
      expect_equal(0)

    expect_lte(max(stat_type_df$esp_dist$label), min(g$gal$n, 10))
    expect(
      min(stat_type_df$esp_dist$esp) >= 0 && max(stat_type_df$esp_dist$esp) <= (g$gal$n^2),
      failure_message = "Some esp counts are out of bounds."
    )
  }
})

test_that("GOF geodesic distance is returned when requested", {
  g <- sim$g

  test_gof_res <-  gof(
    sim$bigergm_res,
    control_within = sim$control_within,
    nsim = 3, compute_geodesic_distance = TRUE
  )

  for (stat_type in c("original", "simulated")) {
    stat_type_df <- test_gof_res[[stat_type]]
    actual_terms <- colnames(stat_type_df$geodesic_dist)

    actual_terms[stringr::str_detect(actual_terms, "nsim", negate = TRUE)] %>%
      setdiff(c("dist", "pairs")) %>%
      length() %>%
      expect_equal(0)

    # Some of the distances will be Inf, and that's ok (that's how ergm returns them).
    non_inf <- stat_type_df$geodesic_dist$dist[!is.infinite(stat_type_df$geodesic_dist$dist)]
    expect_lte(max(non_inf), g$gal$n)
    expect(
      (min(stat_type_df$geodesic_dist$pairs) >= 0) && (max(stat_type_df$geodesic_dist$pairs) <= (g$gal$n^2)),
      failure_message = "Some geodesic distance pairs are out of bounds."
    )
  }

})

test_that("Return GOF statistics including only within-block connections", {
  g <- sim$g
  test_gof_res <-gof(
    sim$bigergm_res,
    control_within = sim$control_within,
    nsim = 3, type = "within", 
    compute_geodesic_distance = FALSE
  )

    # check that the network stats belong to the within-block sub network only
  edgelist <- network::as.edgelist(g) %>% as.data.frame
  colnames(edgelist) <- c('src', 'dst')
  nodes_with_blocks <- data.frame(id = 1:length(network::network.vertex.names(g)), block=sim$bigergm_res$block)
  actual_within_conns <- edgelist %>%
    dplyr::left_join(nodes_with_blocks, by = c('src' = 'id')) %>%
    dplyr::left_join(nodes_with_blocks, by = c('dst' = 'id'), suffix=c('.src', '.dst')) %>%
    dplyr::filter(block.src == block.dst) %>%
    nrow

  within_conns_from_gof <- (test_gof_res$original$network_stats %>% dplyr::filter(stat == 'edges'))[, 2]

  expect_equal(within_conns_from_gof, actual_within_conns)

  
  for (stat_type in c("original", "simulated")) {
    stats <- test_gof_res[[stat_type]]
    expect_false(is.null(stats))
    for (stat in c("network_stats", "degree_dist", "esp_dist")) {
      expect_false(is.null(stats[[stat]]))
    }
    expect_true(is.null(stats[["geodesic_dist"]]))
  }
  })

test_that("Within-connections GOF can be started from the observed network", {
  g <- sim$g

  control_within <- ergm::control.simulate.formula(
    MCMC.burnin = 0,
    MCMC.interval = 1
  )
  test_gof_res <-gof(
    sim$bigergm_res,
    control_within = control_within,
    start_from_observed = TRUE,
    nsim = 2, type = "within"
  )

  first_simulation_stats <- test_gof_res$simulated$network_stats %>%
    dplyr::filter(nsim == 1) %>%
    dplyr::select(-nsim)

  original_network_stats <- test_gof_res$original$network_stats
  expect_equal(original_network_stats,first_simulation_stats)

  })

test_that("Full GOF can be started from the observed network", {
  sim <- get_dummy_net(100, 4)
  g <- sim$g

  control_within <- ergm::control.simulate.formula(
    MCMC.burnin = 0,
    MCMC.interval = 1
  )

  test_gof_res <- gof(
    sim$bigergm_res, 
    type = 'full',
    control_within = control_within,
    nsim = 2,
    start_from_observed = TRUE
  )

  first_simulation_stats <-test_gof_res$simulated$network_stats %>%
    dplyr::filter(nsim == 1)

  # If it starts from the observed network, the stats should not be zero
  expect_true(all(first_simulation_stats['value'] > 0))
})

