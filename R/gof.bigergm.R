# Extracts the degree distribution from a network and returns it as a data frame.
# @param net a statnet network object
# @return a data frame
to_degree_dist_df <- function(net) {
  if(net$gal$directed){
    degrees_in_tmp <- unlist(lapply(net$iel, length))
    degrees_out_tmp <- unlist(lapply(net$oel, length))
    
    data.frame(
      degree = c(0:max(degrees_in_tmp), 0:max(degrees_out_tmp)),
      share = c(as.vector(table(factor(degrees_in_tmp, levels = 0:max(degrees_in_tmp)))/length(degrees_in_tmp)), 
                as.vector(table(factor(degrees_out_tmp, levels = 0:max(degrees_out_tmp)))/length(degrees_out_tmp))),
      type = c(rep("in", max(degrees_in_tmp)+1), rep("out", max(degrees_out_tmp)+1))
    )
  } else {
    degrees_tmp <- unlist(lapply(net$oel, length)) + unlist(lapply(net$iel, length))
    data.frame(
      degree = 0:max(degrees_tmp),
      share = as.vector(table(factor(degrees_tmp, levels = 0:max(degrees_tmp)))/length(degrees_tmp))
    )
  }
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
    pairs = dist
  )
  geodesic_distances_df[geodesic_distances_df$pairs != 0,]
}

# Extracts the edgewise shared partners distribution (undirected).
# @param net a statnet network object
# @return a data frame
to_edgewise_shared_partners_df <- function(net, to = NULL) {
  if(is.null(to)){
    to <- min(net$gal$n - 2, 10)
  }
  esp_dist <- summary(net ~ esp(1:to))
  labels <- as.numeric(gsub("[^0-9]+", "", names(esp_dist)))
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
separate_formulas <- function(target_formula, network = ergm.getnetwork(target_formula)) {
  str_net <- as.character(target_formula)[2]
  net <- get(str_net, envir = environment(target_formula))
  terms <- ergm::ergm_model(target_formula,nw = network)$terms
  varnames <- list_rhs.formula(target_formula) %>% as.character()
  dep_terms <-
    terms %>% purrr::map(function(t) {
      dep <- t$dependence
      is_dep <- is.null(dep) || dep
    }) %>% unlist()
  between_rhs <- varnames[!dep_terms]
  between_rhs <- between_rhs[!is.na(between_rhs)]
  if(length(between_rhs) == 0){
    between_rhs <- "edges"
  }
  within_rhs <- varnames

  between_formula <- paste(str_net, "~", paste(between_rhs, collapse = " + "))
  within_formula <- paste(str_net, "~", paste(within_rhs, collapse = " + "))

  list(
    between = formula(between_formula, env = environment(target_formula)),
    within = formula(within_formula, env = environment(target_formula))
  )
}

# Function to normalize values in a dataset by the observed values
normalize_values <- function(data, obs_data) {
  # Find unique stats in the dataset
  stats <- unique(data$stat)
  # Loop through each stat and normalize its values
  for (stat in stats) {
    # Subset the data for the current stat
    stat_data <- data[data$stat == stat, ]
    normalized_values <- (stat_data$value - obs_data$value[obs_data$stat == stat])/sqrt(var(stat_data$value))
    
    # Update the original dataset with normalized values
    data[data$stat == stat, "value"] <- normalized_values
  }
  
  return(data)
}


# Gets the GOF stats for a formula
# If a network is passed, that one is used to obtain the network statistics,
# otherwise the netwok in the formula is used.
# @param sim_formula a formula
# @param net a statnet network object of the entire network 
# @param net_within a statnet network object of the within-block network
# @param sim_number the ID of the current simulation
# @param compute_geodesic_distance if TRUE, includes the geodesic distance in the result object
# @param to numeric value indicating to which level should the edgewise shared partners be computed
# @return a list with the goodness-of-fit statistics
get_gof_stats <- function(sim_formula, net = NULL,net_within = NULL, sim_number = NULL, compute_geodesic_distance = FALSE, to = NULL) {
  stats_formula <- sim_formula
  if (!is.null(net_within)) {
    stats_formula <- swap_formula_network(net_within, stats_formula, environment())
  }
  network_stats <- summary(stats_formula)
  network_stats <- data.frame(stat = names(network_stats), value = network_stats)
  rownames(network_stats) <- NULL

  # formula_net <- get(as.character(stats_formula)[2], envir = environment(stats_formula))
  degree_dist <- to_degree_dist_df(net)
  esp_dist <- to_edgewise_shared_partners_df(net,to = to)

  if (compute_geodesic_distance == TRUE) {
    geodesic_dist <- to_geodesic_dist_df(net)
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
    stats$network_stats$nsim <- sim_number
    stats$degree_dist$nsim <- sim_number
    stats$esp_dist$nsim <- sim_number

    if (!is.null(stats$geodesic_dist)) {
      stats$geodesic_dist$nsim <- sim_number
    }
  }
  stats
}


#' Conduct Goodness-of-Fit Diagnostics on a Exponential Family Random Graph
#' Model for big networks 
#' 
#' @description
#' A sample of graphs is randomly drawn from the specified model.  The first
#' argument is typically the output of a call to \code{\link{bigergm}} and the
#' model used for that call is the one fit.
#'
#' By default, the sample consists of 100 simulated networks, but this sample
#' size (and many other settings) can be changed using the \code{ergm_control} 
#' argument described above.
#' @param object An \code{\link{bigergm}} object.
#' @param ... Additional arguments, to be passed to \code{\link{simulate_bigergm}}, 
#' which, in turn, passes the information to \code{\link[ergm]{simulate_formula}}.
#' See documentation for \code{\link{bigergm}}.
#' @examples
#' data(toyNet)
#' \donttest{
#' # Specify the model that you would like to estimate.
#'data(toyNet)
#'# Specify the model that you would like to estimate.
#'model_formula <- toyNet ~ edges + nodematch("x") + nodematch("y") + triangle
#'estimate <- bigergm(model_formula,n_blocks = 4)
#'gof_res <- gof(estimate,
#'nsim = 100
#')
#'plot(gof_res)
#' }
#' @param type the type of evaluation to perform. Can take the values `full` or `within`. `full` performs the evaluation on all edges, and `within` only considers within-block edges.
#' @param control_within MCMC parameters as an instance of `control.simulate.formula` to be used for the within-block simulations.
#' @param seed the seed to be passed to simulate_bigergm. If `NULL`, a random seed is used.
#' @param nsim the number of simulations to employ for calculating goodness of fit, default is 100.
#' @param compute_geodesic_distance if `TRUE`, the distribution of geodesic distances is also computed (considerably increases computation time on large networks. `FALSE` by default.)
#' @param start_from_observed if `TRUE`, MCMC uses the observed network as a starting point. If `FALSE`, MCMC starts from a random network.
#' @param simulate_sbm if `TRUE`, the between-block connections are simulated from the estimated stochastic block model from the first stage not the estimated ERGM.
#' @return \code{\link{gof.bigergm}} returns a list with two entries. 
#' The first entry 'original' is another list of the network stats, degree distribution, edgewise-shared partner distribution, and geodesic distance distribution (if \code{compute_geodesic_distance = TRUE}) of the observed network. 
#' The second entry is called 'simulated' is also list compiling the network stats, degree distribution, edgewise-shared partner distribution, and geodesic distance distribution (if \code{compute_geodesic_distance = TRUE}) of all simulated networks.  
#' @export
gof.bigergm <- function(object, ...,  
                        type = 'full',
                        control_within = ergm::control.simulate.formula(),
                        seed = NULL,
                        nsim = 100,
                        compute_geodesic_distance = TRUE,
                        start_from_observed = TRUE, 
                        simulate_sbm = FALSE) {
  
  if(is.null(object$checkpoint)){
    stop("The object must have a checkpoint to be used in the gof function")
  }
  # Setup
  gof_formula <- object$formula
  coef_within_block <- coef(object$est_within)
  coef_between_block <- coef(object$est_between)
  # Validate the simulation type
  allowed_type_values <- c('full', 'within')
  if (!type %in% allowed_type_values){
    stop("The `type` argument must be any of 'full' or 'within'")
  }
  
  if (type == 'full'){
    if(start_from_observed){
      net <- object$checkpoint$network
      net_within <- get_within_networks(net, net %v% "block")
    } else {
      net_tmp <- object$checkpoint$network
      net <- network.initialize(n = net_tmp$gal$n, directed = net_tmp$gal$directed)
      # Copy all vertex.attributes from net_tmp to net
      net <- copy_vertex_attributes(net = net, net_tmp = net_tmp)
      net_within <- get_within_networks(net, net %v% "block")
      rm(net_tmp)
    }
    original_stats <- get_gof_stats(object$formula, net = object$checkpoint$network,
                                    net_within =  net_within, 
                                    compute_geodesic_distance = compute_geodesic_distance)
    
    
  } else {
    net <- object$checkpoint$network
    net_within <- get_within_networks(net, net %v% "block")
    net <- get_within_networks(net, net %v% "block",combined_networks = FALSE)
    original_stats <- get_gof_stats(gof_formula, net = net,
                                    net_within = net_within, 
                                    compute_geodesic_distance = compute_geodesic_distance)
    if(!start_from_observed){
      net <- copy_vertex_attributes(net = network.initialize(n = net_within$gal$n,
                                                             directed = net_within$gal$directed), 
                                    net_tmp = net_within)
      net_within <- get_within_networks(net, net %v% "block")
    }
  }
  # gof_formula <- swap_formula_network(net_formula = gof_formula, new_net = net, env = environment())
  # Simulate all networks
  if (type == 'full'){
    simulated_networks <- simulate_bigergm(
      formula = gof_formula,
      coef_within = coef_within_block,
      coef_between = coef_between_block,
      network = net, 
      control_within = control_within,
      seed = seed,
      nsim = nsim,
      output = "network", 
      within_formula = object$est_within$within_formula,
      add_intercepts = object$checkpoint$add_intercepts,
      clustering_with_features = object$checkpoint$clustering_with_features, 
      sbm_pi = object$sbm_pi, 
      simulate_sbm = simulate_sbm, ...
    )
  } else {
    simulated_networks <- simulate_bigergm(
      formula = gof_formula,
      coef_within = coef_within_block,
      coef_between = coef_between_block,
      control_within = control_within,
      network = net, 
      seed = seed,
      nsim = nsim,
      output = "network", 
      only_within = TRUE, 
      within_formula = object$est_within$within_formula,
      add_intercepts = object$checkpoint$add_intercepts,
      clustering_with_features = object$checkpoint$clustering_with_features, 
      sbm_pi = object$sbm_pi, 
      simulate_sbm = simulate_sbm, 
      return_within = TRUE, ...
    )
  }

  if(control_within$parallel>1){
    simulated <- parallel::mcmapply(FUN = function(x,y, compute_geodesic_distance, original_stats) {
      get_gof_stats(gof_formula, net = x,net_within =  get_within_networks(x, x %v% "block"),
                    sim_number = y, 
                    compute_geodesic_distance = compute_geodesic_distance,
                    to = original_stats$esp_dist$label[max(which(original_stats$esp_dist$esp != 0))])
    }, x = simulated_networks, y= 1:length(simulated_networks),
    MoreArgs = list(compute_geodesic_distance = compute_geodesic_distance,  original_stats = original_stats),
    mc.cores=control_within$parallel, SIMPLIFY = FALSE)
  } else {
    simulated <-  mapply(FUN = function(x,y, compute_geodesic_distance, original_stats) {
      get_gof_stats(gof_formula, net = x,net_within =  get_within_networks(x, x %v% "block"),
                    sim_number = y, 
                    compute_geodesic_distance = compute_geodesic_distance,
                    to = original_stats$esp_dist$label[max(which(original_stats$esp_dist$esp != 0))])
    }, x = simulated_networks, y= 1:length(simulated_networks), 
    MoreArgs = list(compute_geodesic_distance = compute_geodesic_distance,  
    original_stats = original_stats),
    SIMPLIFY = FALSE)
  }
  
  
  
  # Get the statistics for all simulated networks

  geodesic_dist <- lapply(simulated, function(x) x$geodesic_dist)
  network_stats <- lapply(simulated, function(x) x$network_stats)
  network_stats <- do.call(rbind, network_stats)
  network_stats_normalized <- normalize_values(data = network_stats, obs_data = original_stats$network_stats)
  degree_dist <- lapply(simulated, function(x) x$degree_dist)
  esp_dist <- lapply(simulated, function(x) x$esp_dist)
  
  
  results <- list(
    original = list(
      network_stats = original_stats$network_stats,
      degree_dist = original_stats$degree_dist,
      esp_dist = original_stats$esp_dist,
      geodesic_dist = original_stats$geodesic_dist
    ),
    simulated = list(
      network_stats = network_stats,
      network_stats_normalized = network_stats_normalized,
      degree_dist = do.call(rbind, degree_dist),
      esp_dist = do.call(rbind, esp_dist),
      geodesic_dist = do.call(rbind, geodesic_dist)
    )
  )
  return(structure(results, class = "gof.bigergm"))
}






#' @export
plot.gof.bigergm <- function(x, ...) {
  if(!is.null(x$simulated$degree_dist$type)){
    
    # In degree
    original_indegree_distribution <- x$original$degree_dist[x$original$degree_dist$type == "in",]
    simulated_indegree_distribution <- x$simulated$degree_dist[x$simulated$degree_dist$type == "in",]
    levels =sort(unique(c(original_indegree_distribution$degree,simulated_indegree_distribution$degree)))
    if(length(levels[!levels %in% original_indegree_distribution$degree])>0){
      original_indegree_distribution <- rbind(original_indegree_distribution, 
                                      cbind(degree = levels[!levels %in% original_indegree_distribution$degree], 
                                            share = rep(0,length(levels[!levels %in% original_indegree_distribution$degree])), 
                                            type = "in"))
      original_indegree_distribution <- original_indegree_distribution[order(as.numeric(original_indegree_distribution$degree)),]
    }
    
    if(length(levels[!levels %in% simulated_indegree_distribution$degree])>0){
      simulated_indegree_distribution <- rbind(simulated_indegree_distribution, 
                                       cbind(degree = levels[!levels %in% simulated_indegree_distribution$degree], 
                                             share = rep(0,length(levels[!levels %in% simulated_indegree_distribution$degree])), 
                                             type = "in",
                                             nsim = 1))
      simulated_indegree_distribution <- simulated_indegree_distribution[order(as.numeric(simulated_indegree_distribution$degree)),]
    }
    simulated_indegree_distribution$degree <- factor(simulated_indegree_distribution$degree, levels =levels)
    original_indegree_distribution$degree <- factor(original_indegree_distribution$degree, levels =levels)
    simulated_indegree_distribution$share <- as.numeric(simulated_indegree_distribution$share)
    original_indegree_distribution$share <- as.numeric(original_indegree_distribution$share)

    boxplot(simulated_indegree_distribution$share~
              simulated_indegree_distribution$degree,
            xlab = "In-Degree", ylab = "Share")
    lines(original_indegree_distribution$share, col = "red")
    # Out degree
    
    original_outdegree_distribution <- x$original$degree_dist[x$original$degree_dist$type == "out",]
    simulated_outdegree_distribution <- x$simulated$degree_dist[x$simulated$degree_dist$type == "out",]
    levels =sort(unique(c(original_outdegree_distribution$degree,simulated_outdegree_distribution$degree)))
    
    if(length(levels[!levels %in% original_outdegree_distribution$degree])>0){
      original_outdegree_distribution <- rbind(original_outdegree_distribution, 
                                               cbind(degree = levels[!levels %in% original_outdegree_distribution$degree], 
                                                     share = rep(0,length(levels[!levels %in% original_outdegree_distribution$degree])), 
                                                     type = "out"))
      original_outdegree_distribution <- original_outdegree_distribution[order(as.numeric(original_outdegree_distribution$degree)),]
    }
    
    if(length(levels[!levels %in% simulated_outdegree_distribution$degree])>0){
      simulated_outdegree_distribution <- rbind(simulated_outdegree_distribution, 
                                                cbind(degree = levels[!levels %in% simulated_outdegree_distribution$degree], 
                                                      share = rep(0,length(levels[!levels %in% simulated_outdegree_distribution$degree])), 
                                                      type = "out",
                                                      nsim = 1))
      simulated_outdegree_distribution <- simulated_outdegree_distribution[order(as.numeric(simulated_outdegree_distribution$degree)),]
    }
    simulated_outdegree_distribution$degree <- factor(simulated_outdegree_distribution$degree, levels =levels)
    original_outdegree_distribution$degree <- factor(original_outdegree_distribution$degree, levels =levels)
    simulated_outdegree_distribution$share <- as.numeric(simulated_outdegree_distribution$share)
    original_outdegree_distribution$share <- as.numeric(original_outdegree_distribution$share)
    boxplot(simulated_outdegree_distribution$share~
              simulated_outdegree_distribution$degree,
            xlab = "Out-Degree", ylab = "Share")
    lines(original_outdegree_distribution$share, col = "red")
  } else {
    levels =sort(unique(c(x$original$degree_dist$degree,x$simulated$degree_dist$degree)))
    if(length(levels[!levels %in% x$original$degree_dist$degree])>0){
      x$original$degree_dist <- rbind(x$original$degree_dist, 
                                        cbind(degree = levels[!levels %in% x$original$degree_dist$degree], 
                                              share = rep(0,length(levels[!levels %in% x$original$degree_dist$degree]))))
      x$original$degree_dist <- x$original$degree_dist[order(x$original$degree_dist$degree),]
    }
    
    if(length(levels[!levels %in% x$simulated$degree_dist$degree])>0){
      x$simulated$degree_dist <- rbind(x$simulated$degree_dist, 
                                      cbind(degree = levels[!levels %in% x$simulated$degree_dist$degree], 
                                            share = rep(0,length(levels[!levels %in% x$simulated$degree_dist$degree])), 
                                            nsim = 1))
      x$simulated$degree_dist <- x$simulated$degree_dist[order(x$simulated$degree_dist$degree),]
    }
    x$simulated$degree_dist$degree <- factor(x$simulated$degree_dist$degree, levels =levels)
    x$original$degree_dist$degree <- factor(x$original$degree_dist$degree, levels =levels)
    boxplot(x$simulated$degree_dist$share~x$simulated$degree_dist$degree, xlab = "Degree", ylab = "Share")
    lines(x$original$degree_dist$share, col = "red")
  }
  # ESP distribution
  levels =sort(unique(c(x$original$esp_dist$label,x$simulated$esp_dist$label)))
  if(length(levels[!levels %in% x$original$esp_dist$label])>0){
    x$original$esp_dist <- rbind(x$original$esp_dist, 
                                      cbind(label = levels[!levels %in% x$original$esp_dist$label], 
                                            esp = rep(0,length(levels[!levels %in% x$original$esp_dist$label]))))
    x$original$esp_dist <- x$original$esp_dist[order(as.numeric(x$original$esp_dist$label)),]
  }
  
  if(length(levels[!levels %in% x$simulated$esp_dist$label])>0){
    x$simulated$esp_dist <- rbind(x$simulated$esp_dist, 
                                 cbind(label = levels[!levels %in% x$simulated$esp_dist$label], 
                                       esp = rep(0,length(levels[!levels %in% x$simulated$esp_dist$label])),
                                       nsim = 0))
    x$simulated$esp_dist <- x$simulated$esp_dist[order(as.numeric(x$simulated$esp_dist$label)),]
  }
  x$simulated$esp_dist$label <- factor(x$simulated$esp_dist$label, levels =levels)
  x$original$esp_dist$label <- factor(x$original$esp_dist$label, levels =levels)
  boxplot(x$simulated$esp_dist$esp ~x$simulated$esp_dist$label, xlab = "Edgewise-shared Partner", ylab = "Number")
  lines(x$original$esp_dist$esp, col = "red")
  
  if(!is.null(x$original$geodesic_dist)){
    levels =sort(unique(c(x$original$geodesic_dist$dist,x$simulated$geodesic_dist$dist)))
    if(length(levels[!levels %in% x$original$geodesic_dist$dist])>0){
      x$original$geodesic_dist <- rbind(x$original$geodesic_dist, 
                                        cbind(dist = levels[!levels %in% x$original$geodesic_dist$dist], 
                                              pairs = rep(0,length(levels[!levels %in% x$original$geodesic_dist$dist]))))
      x$original$geodesic_dist <- x$original$geodesic_dist[order(x$original$geodesic_dist$dist),]
    }
    
    if(length(levels[!levels %in% x$simulated$geodesic_dist$dist])>0){
      x$simulated$geodesic_dist <- rbind(x$simulated$geodesic_dist, 
                                        cbind(dist = levels[!levels %in% x$simulated$geodesic_dist$dist], 
                                              pairs = rep(0,length(levels[!levels %in% x$simulated$geodesic_dist$dist])), 
                                              nsim = 1))
      x$simulated$geodesic_dist <- x$simulated$geodesic_dist[order(x$simulated$geodesic_dist$dist),]
    }
    x$simulated$geodesic_dist$dist <- factor(x$simulated$geodesic_dist$dist, levels =levels)
    x$original$geodesic_dist$dist <- factor(x$original$geodesic_dist$dist, levels =levels)
    boxplot(x$simulated$geodesic_dist$pairs ~x$simulated$geodesic_dist$dist, xlab = "Geodesic Distance", ylab = "Number")
    lines(x$original$geodesic_dist, col = "red")
  }
  boxplot(x$simulated$network_stats_normalized$value~x$simulated$network_stats_normalized$stat,
          xlab = "Statistics", ylab = "Normalized simulated statistics")
  abline(a = 0, b = 0, col = "red")
}

