# Extracts the degree distribution from a network and returns it as a data frame.
# @param net a statnet network object
# @return a data frame
to_degree_dist_df <- function(net) {
  degree_dist <- igraph::degree.distribution(intergraph::asIgraph(net))
  data.frame(
    degree = 0:(length(degree_dist) - 1),
    share = degree_dist
  )
}

# Extracts the geodesic distance distribution from a network and returns it as a dataframe.
# @param net a statnet network object
# @return a data frame
to_geodesic_dist_df <- function(net) {
  dist <- ergm::ergm.geodistdist(net)
  labels <- as.numeric(names(dist))
  names(dist) <- NULL
  geodesic_distances_df <- data.frame(
    dist = labels,
    pairs = dist[[1]]
  )
  geodesic_distances_df
}

# Extracts the edgewise shared partners distribution (undirected).
# @param net a statnet network object
# @return a data frame
to_edgewise_shared_partners_df <- function(net) {
  esp_dist <- summary(net ~ esp(1:min(net$gal$n - 2, 10)))
  labels <- names(esp_dist) %>%
    purrr::map(function(lab) {
      stringr::str_replace(lab, "esp", "")
    }) %>%
    as.numeric()

  names(esp_dist) <- NULL

  data.frame(
    label = labels,
    esp = esp_dist
  )
}

# Swaps the network on the lhs of a formula for a new one with the given environment
# @param new_net A network object to be inserted into the lhs of the formula
# @param net_formula The target formula
# @param env The environment to assign to the formula
# @return A new formula with the lhs swapped
swap_formula_network <- function(new_net, net_formula, env) {
  rhs <- as.character(net_formula)[3]
  as.formula(paste(deparse(substitute(new_net)), "~", rhs), env = env)
}

# Separates a formula into its between and within components. The between component excludes
# terms which introduce dyadic dependence.
# @param target_formula a target formula
# @return a list containing the between and within formulas
separate_formulas <- function(target_formula) {
  str_net <- as.character(target_formula)[2]
  net <- get(str_net, envir = environment(target_formula))
  terms <- ergm::ergm_model(target_formula)$terms
  varnames <- list_rhs.formula(target_formula) %>% as.character()
  dep_terms <-
    terms %>% purrr::map(function(t) {
      dep <- t$dependence
      is_dep <- is.null(dep) || dep
    }) %>% unlist()
  between_rhs <- varnames[!dep_terms]
  within_rhs <- varnames

  between_formula <- paste(str_net, "~", paste(between_rhs, collapse = " + "))
  within_formula <- paste(str_net, "~", paste(within_rhs, collapse = " + "))

  list(
    between = formula(between_formula, env = environment(target_formula)),
    within = formula(within_formula, env = environment(target_formula))
  )
}

# Gets the GOF stats for a formula
# If a network is passed, that one is used to obtain the network statistics,
# otherwise the netwok in the formula is used.
# @param sim_formula a formula
# @param net a statnet network object
# @param sim_number the ID of the current simulation
# @param compute_geodesic_distance if TRUE, includes the geodesic distance in the result object
# @return a list with the goodness-of-fit statistics
get_gof_stats <- function(sim_formula, net = NULL, sim_number = NULL, compute_geodesic_distance = FALSE) {
  stats_formula <- sim_formula

  if (!is.null(net)) {
    stats_formula <- swap_formula_network(net, stats_formula, environment())
  }

  network_stats <- summary(stats_formula)
  network_stats <- data.frame(stat = names(network_stats), value = network_stats)
  rownames(network_stats) <- NULL

  formula_net <- get(as.character(stats_formula)[2], envir = environment(stats_formula))
  degree_dist <- to_degree_dist_df(formula_net)
  esp_dist <- to_edgewise_shared_partners_df(formula_net)

  if (compute_geodesic_distance == TRUE) {
    geodesic_dist <- to_geodesic_dist_df(formula_net)
  } else {
    geodesic_dist <- NULL
  }

  stats <- list(
    network_stats = network_stats,
    degree_dist = degree_dist,
    esp_dist = esp_dist,
    geodesic_dist = geodesic_dist
  )

  if (!is.null(sim_number)) {
    stats$network_stats$n_sim <- sim_number
    stats$degree_dist$n_sim <- sim_number
    stats$esp_dist$n_sim <- sim_number

    if (!is.null(stats$geodesic_dist)) {
      stats$geodesic_dist$n_sim <- sim_number
    }
  }

  stats
}


#' Goodness of fit statistics for HERGM
#' @param net the target network
#' @param data_for_simulation a dataframe with node-level covariates
#' @param list_feature_matrices a list of feature adjacency matrices
#' @param colname_vertex_id the name of the column that contains the node id
#' @param colname_block_membership the name o the column that contains the block affiliation of each node
#' @param bigergm_results a bigergm results object
#' @param type the type of evaluation to perform. Can take the values `full` or `within`. `full` performs the evaluation on all edges, and `within` only considers within-block edges.
#' @param ergm_control MCMC parameters as an instance of ergm.control
#' @param seed the seed to be passed to simulate_hergm
#' @param n_sim the number of simulations to employ for calculating goodness of fit
#' @param prevent_duplicate see `simulate_hergm`
#' @param compute_geodesic_distance if `TRUE`, the distribution of geodesic distances is also computed (considerably increases computation time on large networks. `FALSE` by default.)
#' @param start_from_observed if `TRUE`, MCMC uses the observed network as a starting point
#' @param ... Additional arguments, to be passed to lower-level functions
#' @return \code{\link{gof_bigergm}} returns a list with two entries. 
#' The first entry 'original' is another list of the network stats, degree distribution, edgewise-shared partner distribution, and geodesic distance distribution (if \code{compute_geodesic_distance = TRUE}) of the observed network. 
#' The second entry is called 'simulated' is also list compiling the network stats, degree distribution, edgewise-shared partner distribution, and geodesic distance distribution (if \code{compute_geodesic_distance = TRUE}) of all simulated networks.  
#' @examples
#' data(toyNet)
#' \donttest{
#' # Specify the model that you would like to estimate.
#' model_formula <- toyNet ~ edges + nodematch("x") + nodematch("y") + triangle
#' # Estimate the model
#' nodes_data <- data.frame(
#'   node_id = 1:toyNet$gal$n,
#'   x = toyNet %v% "x",
#'   y = toyNet %v% "y",
#'   block = toyNet %v% "block"
#' )
#' list_feature_matrices <- bigergm::get_list_sparse_feature_adjmat(toyNet, model_formula)
#' estimate <- hergm(model_formula,n_clusters = 4)
#' gof_res <- bigergm::gof_bigergm(
#'   toyNet,
#'   list_feature_matrices = list_feature_matrices,
#'   data_for_simulation = nodes_data,
#'   colname_vertex_id = "node_id",
#'   colname_block_membership = "block",
#'   bigergm_results = estimate,
#'   n_sim = 100
#' )
#' }
#' 
#' @export
gof_bigergm <- function(net,
                           data_for_simulation,
                           list_feature_matrices,
                           colname_vertex_id,
                           colname_block_membership,
                           bigergm_results,
                           type = 'full',
                           ergm_control = ergm::control.simulate.formula(),
                           seed = NULL,
                           n_sim = 1,
                           prevent_duplicate = TRUE,
                           compute_geodesic_distance = FALSE,
                           start_from_observed = FALSE,
                           ...) {
  # Setup
  gof_formula <- swap_formula_network(net, bigergm_results$est_within$formula, environment())
  coef_within_block <- coef(bigergm_results$est_within)
  coef_between_block <- coef(bigergm_results$est_between)

  # Validate the simulation type
  allowed_type_values <- c('full', 'within')
  if (!type %in% allowed_type_values){
    stop("The `type` argument must be any of 'full' or 'within'")
  }

  seed_edgelist = NULL

  if (type == 'full'){
    original_stats <- get_gof_stats(gof_formula, compute_geodesic_distance = compute_geodesic_distance)

    if(start_from_observed){
      seed_edgelist <- network::as.edgelist(net)
    }

  } else {
    sorted_dataframe <- sort_block_membership(data_for_simulation, colname_vertex_id, colname_block_membership)
    seed_edgelist_within <- arrange_edgelist(network::as.edgelist(net), sorted_dataframe)$edgelist_within
    within_network <- generate_seed_network(gof_formula, sorted_dataframe, edgelist = seed_edgelist_within, directed = FALSE)

    original_stats <- get_gof_stats(gof_formula, net = within_network, compute_geodesic_distance = compute_geodesic_distance)

    if(start_from_observed){
      seed_edgelist <- network::as.edgelist(within_network)
    }
  }

  # Simulate the first network by initializing it from zero. The burnin here is the one set by the user.
  if (type == 'full'){
    base_network <- simulate_hergm(
      formula_for_simulation = gof_formula,
      data_for_simulation = data_for_simulation,
      colname_vertex_id = colname_vertex_id,
      colname_block_membership = colname_block_membership,
      seed_edgelist = seed_edgelist,
      coef_within_block = coef_within_block,
      coef_between_block = coef_between_block,
      ergm_control = ergm_control,
      seed_for_within = seed,
      seed_for_between = seed,
      directed = FALSE,
      n_sim = 1,
      output = "network",
      prevent_duplicate = prevent_duplicate,
      list_feature_matrices = list_feature_matrices,
      use_fast_between_simulation = TRUE
    )
  } else {
    base_network <- bigergm::simulate_hergm_within(
      formula_for_simulation = gof_formula,
      data_for_simulation = data_for_simulation,
      colname_vertex_id = colname_vertex_id,
      colname_block_membership = colname_block_membership,
      seed_edgelist = seed_edgelist,
      coef_within_block = coef_within_block,
      output = 'network',
      ergm_control = ergm_control,
      seed = seed,
      n_sim = 1
    )
  }

  # Get the statistics for the first network
  sim_stats <- get_gof_stats(gof_formula, base_network, 1, compute_geodesic_distance = compute_geodesic_distance)
  results <- list(
    original = list(
      network_stats = original_stats$network_stats,
      degree_dist = original_stats$degree_dist,
      esp_dist = original_stats$esp_dist,
      geodesic_dist = original_stats$geodesic_dist
    ),
    simulated = list(
      network_stats = sim_stats$network_stats,
      degree_dist = sim_stats$degree_dist,
      esp_dist = sim_stats$esp_dist,
      geodesic_dist = sim_stats$geodesic_dist
    )
  )

  effective_nsim <- n_sim - 1

  if (effective_nsim > 0) {
    # Now replace the burnin with the interval and simulate networks one by one.
    ergm_control$MCMC.burnin <- ergm_control$MCMC.interval

    for (i in 1:effective_nsim) {
      if ((i + 1) %% 50 == 0) {
        message(paste("Simulation:", i + 1))
      }

      if(type == 'full'){
        base_network <- simulate_hergm(
          formula_for_simulation = gof_formula,
          list_feature_matrices = list_feature_matrices,
          data_for_simulation = data_for_simulation,
          colname_vertex_id = colname_vertex_id,
          colname_block_membership = colname_block_membership,
          seed_edgelist = network::as.edgelist(base_network),
          coef_within_block = coef_within_block,
          coef_between_block = coef_between_block,
          ergm_control = ergm_control,
          seed = seed,
          directed = FALSE,
          n_sim = 1,
          output = "network",
          prevent_duplicate = prevent_duplicate,
          use_fast_between_simulation = TRUE,
          ...
        )
      } else {
        base_network <- bigergm::simulate_hergm_within(
          formula_for_simulation = gof_formula,
          data_for_simulation = data_for_simulation,
          colname_vertex_id = colname_vertex_id,
          colname_block_membership = colname_block_membership,
          seed_edgelist = network::as.edgelist(base_network),
          coef_within_block = coef_within_block,
          output = 'network',
          ergm_control = ergm_control,
          seed = seed,
          n_sim = 1
        )
      }

      sim_stats <- get_gof_stats(gof_formula, base_network, i + 1, compute_geodesic_distance = compute_geodesic_distance)
      results$simulated$network_stats <- rbind(results$simulated$network_stats, sim_stats$network_stats)
      results$simulated$degree_dist <- rbind(results$simulated$degree_dist, sim_stats$degree_dist)
      results$simulated$esp_dist <- rbind(results$simulated$esp_dist, sim_stats$esp_dist)
      if (
        !(is.null(results$simulated$geodesic_dist)) &
          !(is.null(sim_stats$geodesic_dist))
      ) {
        results$simulated$geodesic_dist <- rbind(results$simulated$geodesic_dist, sim_stats$geodesic_dist)
      }
    }

    message("Simulation Finished")
  }

  return(results)
}
