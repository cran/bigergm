#' Simulate networks under Exponential Random Graph Models (ERGMs) under local dependence
#' @description
#' This function simulates networks under Exponential Random Graph Models (ERGMs) with local dependence.
#' There is also an option to simulate only within-block networks and a S3 method for the class `bigergm`.
#' 
#' @param formula  An \R \code{\link{formula}} object of the form
#' \code{y ~ <model terms>}, where \code{y} is a
#' \code{\link[network]{network}} object. 
#' The network object must contain block information as a vertex attribute with the name 'block'.
#' For the details on the possible \code{<model terms>}, see
#' \code{\link[ergm]{ergmTerm}} and Morris, Handcock and Hunter (2008).
#' All terms that induce dependence are excluded from the between block model, while the within block model includes all terms.
#' The \code{\link[ergm.multi]{L-ergmTerm}} is supported to enable size-dependent coefficients for the within-blocks model. 
#' Note, however, that for size-dependent parameters of terms that are included in the between-blocks model, 
#' the intercept in the linear model provided to \code{\link[ergm.multi]{L-ergmTerm}} should not include the intercept. 
#' See the second example of \code{\link{bigergm}} for a demonstration. 
#' @param coef_within a vector of within-block parameters. The order of the parameters should match that of the formula.
#' @param coef_between a vector of between-block parameters. The order of the parameters should match that of the formula without externality terms.
#' @param control_within auxiliary function as user interface for fine-tuning ERGM simulation for within-block networks.
#' @param seed seed value (integer) for network simulation.
#' @param nsim number of networks generated.
#' @param network a network object to be used as a seed network for the simulation (if none is provided, the network on the lhs of the `formula` is used).
#' @param only_within If this is TRUE, only within-block networks are simulated.
#' @param output Normally character, one of "network" (default), "stats", "edgelist", to determine the output format.
#' @param verbose If this is TRUE/1, the program will print out additional information about the progress of simulation.
#' @param ... Additional arguments, passed to \code{\link[ergm]{simulate_formula}}.
#'
#' @references 
#' Morris M, Handcock MS, Hunter DR (2008). Specification of Exponential-Family Random Graph Models: Terms and Computational Aspects.
#' Journal of Statistical Software, 24.
#' 
#' @examples
#' data(toyNet)
#' # Specify the model that you would like to estimate.
#' model_formula <- toyNet ~ edges + nodematch("x") + nodematch("y") + triangle
#' # Simulate network stats
#' sim_stats <- bigergm::simulate_bigergm(
#' formula = model_formula,
#'   # Formula for the model
#' coef_between = c(-4.5,0.8, 0.4),
#' # The coefficients for the between connections
#' coef_within = c(-1.7,0.5,0.6,0.15),
#' # The coefficients for the within connections
#' nsim = 10,
#' # Number of simulations to return
#' output = "stats",
#' # Type of output
#' )
#'
#' @return Simulated networks, the output form depends on the parameter \code{output} 
#' (default is a list of networks).
#' @export
#' @importFrom Matrix sparseMatrix
#' @importFrom stats as.formula coef formula kmeans
#' @importFrom methods as
simulate_bigergm <- function(formula,
                           coef_within,
                           coef_between,
                           network = ergm.getnetwork(formula), 
                           control_within = ergm::control.simulate.formula(),
                           only_within = FALSE,
                           seed = NULL,
                           nsim = 1,
                           output = "network",
                           verbose = 0,
                           ...) {
  # Save the ... parameters in a list
  list_info <- list(...)
  # For the entire simulation process the vertex.names are changed to be integers running from 1 to n 
  old_vertex_names <- network::network.vertex.names(network)
  network::network.vertex.names(network) <- 1:length(old_vertex_names)
  
  
  if(NA %in% get.vertex.attribute(network,attrname = "block")){
    stop("The network object does not contain block information. Please add block information to the network object as a vertex attribute with the name 'block'.")
  }

  # If output == "network":
  if (!is.function(output) && output == "network") {
    temp <- "edgelist"
  } else {
    temp <- output
  }

  # Create formula for simulating within- and between-block networks.
  ## Extract the RHS of the given formula
  formula_rhs <- as.character(formula)[3]
  ## Create a formula for simulation
  ### For within-block network
  formula_for_simulation_within <- as.formula(glue::glue("network ~ {formula_rhs}"))

  ### For between-block network
  # Exclude all N() terms for the between-network model
  ### Here, the LHS network object can be arbitrary as long as the RHS is correct.
  network_blocked <- get_within_networks(network, network%v%"block")
  formula_for_simulation_between <- separate_formulas(formula_for_simulation_within, network = network_blocked)[[1]]
  term_names <- rownames(attr(terms(formula_for_simulation_between),"factors"))
  if(sum(grepl("N\\(",term_names))>0){
    formula_for_simulation_between <- update(formula_for_simulation_between, 
                                            new = paste("~. -",paste(term_names[grepl("N\\(",term_names)], 
                                                                     collapse = "-"))) 
  }
  # Exclude the offset terms
  formula_for_simulation_between <- no.offset(formula_for_simulation_between)
  if("within_formula" %in% names(list_info)){
    formula_for_simulation_within <- list_info$within_formula
  }
  # Send a message about which value of a coefficient is attached to which term.
  # This is to make sure that `coef_within_block` and `coef_between_block` are correctly specified as the user intends.
  attach_terms_to_coef <- function(formula, coef) {
    terms <- as.character(formula)[3]
    terms <- unlist(stringr::str_split(string = terms, pattern = " \\+ "))
    terms <- stringr::str_replace_all(terms, "\"", "'")
    names(coef) <- terms
    return(coef)
  }
  coefs <- purrr::map2(
    list(formula_for_simulation_within, formula_for_simulation_between),
    list(coef_within, coef_between),
    attach_terms_to_coef
  )

  # Print out the specified coefficients
  if (verbose > 0) {
    message("Specified coefficients for the within-block model:")
    print(coefs[[1]])

    message("Specified coefficients for the between-block model:")
    print(coefs[[2]])
  }

  ############################################################################
  ############# Draw within-block connections ################################
  ############################################################################

  if (verbose > 0) {
    message("Simulating within-block networks.")
  }
  
  edgelist_within <- draw_within_block_connection(
    seed_network = network_blocked,
    formula_for_simulation = formula_for_simulation_within,
    coef_within_block = coef_within,
    ergm_control = control_within,
    output = temp,
    seed = seed,
    nsim = nsim,
    verbose = verbose, 
    ...
  )
  
  ############################################################################
  ############# Draw between-block connections ###############################
  ############################################################################

  if (verbose > 0) {
    message("Simulating between-block networks.")
  }
  edgelist_between <-
    draw_between_block_connection(
      network = network,
      formula_for_simulation = formula_for_simulation_between,
      coef_between_block = coef_between,
      seed = seed,
      nsim = nsim,
      verbose = verbose,
      draw_networks = !only_within, ...
    )
  res <- list()
  
  # if("return_within" %in% names(list_info)){
  #   if(list_info$"return_within"){
  #     browser()  
  #     if (nsim == 1) {
  #       
  #     }
  #     res <- append(res,edgelist_within)
  #     generate_network_for_output(edgelist = edgelist, 
  #                                 sorted_dataframe = edgelist_between$node_data,
  #                                 old_vertex_names)
  #     
  #   }
  # }
  
  if((nsim>1)& (output == "edgelist")){
    # Shorten if number of effective samples is smaller than nsim
    edgelist_between$output <- edgelist_between$output[1:length(edgelist_within)]
  }
  if (output %in% c('network', 'edgelist')) {
    if (nsim == 1) {
      # Combine within- and between-block edges while removing duplicated links.
      edgelist <- combine_within_between_edges(edgelist_within, edgelist_between$output, TRUE, old_vertex_names)
      if (output == 'edgelist'){
        simulation_output <- edgelist
      } else {
        # Create a final network
        simulation_output <- generate_network_for_output(edgelist = edgelist, 
                                                         sorted_dataframe = edgelist_between$node_data, 
                                                         old_vertex_names)
      }
      # Return the output
      if (verbose > 0) {
        message("One entire network has been generated.")
      }
      attr(simulation_output,"vnames") <- old_vertex_names
      return(simulation_output)
    } else {
      # Combine within- and between-block edges while removing duplicated links.
      fn_combine <- function(edgelist_within, edgelist_between) {
        combine_within_between_edges(edgelist_within, edgelist_between, TRUE, old_vertex_names)
      }
      edgelist <- purrr::map2(edgelist_within, edgelist_between$output, fn_combine)
      # Create a final network
      if (output == 'edgelist'){
        simulation_output <- edgelist
      } else {
        # browser()
        fn_final_network <- function(edgelist) {
          return(generate_network_for_output(edgelist = edgelist, 
                                             sorted_dataframe = edgelist_between$node_data,old_vertex_names = old_vertex_names))
        }
        simulation_output <- purrr::map(edgelist, fn_final_network)
      }

      # Return the output
      if (verbose > 0) {
        message(glue::glue("{nsim} entire networks have been generated."))
      }
      res <- append(res,simulation_output)
    }
  } else if (output == "stats") {
    output_within <- data.frame(edgelist_within)
    if(only_within){
      output_between <- NULL
    } else {
      output_between <- get_between_stats(edgelist_between$output,
                                          between_formula = formula_for_simulation_between)
    }
    res <- append(res,list(
      within_network = output_within,
      between_network = output_between
    ))
  }
  # In this case the result depends on the function set as output
  else {
    cat("here")
    res <- append(res,list(
      within_network = edgelist_within,
      between_network = edgelist_between$output
    ))
  }
  return(res)
  }
#' Simulate networks under Exponential Random Graph Models (ERGMs) under local dependence
#' @description
#' This function simulates networks under the Exponential Random Graph Model (ERGM) 
#' with local dependence with all parameters set according to the estimated model (`object`).
#' See \code{\link{simulate_bigergm}} for details of the simulation process 
#' @param object an object of class `bigergm`
#' @param seed seed value (integer) for network simulation.
#' @param verbose If this is TRUE/1, the program will print out additional information about the progress of simulation.
#' @param output Normally character, one of "network" (default), "stats", "edgelist", to determine the output of the function.
#' @param control_within \code{\link[ergm]{control.simulate.formula}} object for fine-tuning ERGM simulation of within-block networks.
#' @param nsim number of networks to be randomly drawn from the given distribution on the set of all networks.
#' @param ... Additional arguments, passed to \code{\link[ergm]{simulate_formula}}.
#'
#' @return Simulated networks, the output form depends on the parameter \code{output} 
#' (default is a list of networks).
#' @export
simulate.bigergm <- function (object, nsim = 1, seed = NULL, ..., output = "network",
                             control_within = ergm::control.simulate.formula(), verbose = 0) {
  if (is.null(seed)) {
    seed <- sample.int(10^5, 1)
  }
  if (is.null(object$control)) {
    object$control <- ergm::control.simulate.formula()
  }
  if(is.null(object$checkpoint)){
    stop("The object must have a checkpoint to be used in the gof function")
  }
  simulate_bigergm(
    formula = object$formula,
    network = object$checkpoint$network,
    coef_within = object$est_within$coefficients,
    coef_between = object$est_between$coefficients,
    control_within = control_within,
    seed = seed,
    nsim = nsim,
    output = output,
    verbose = verbose,
    within_formula = object$est_within$within_formula,
    add_intercepts = object$checkpoint$add_intercepts,
    clustering_with_features = object$checkpoint$clustering_with_features, 
    sbm_pi = object$sbm_pi
  )
}
# Extract Covariate Names
# Extracts the names of covariates used in the formula
# @param formula_for_simulation the formula to check
# @return A list of covariate names (can be empty)
#' @importFrom magrittr %>%
extract_covariate_names <- function(formula_for_simulation){
  covars_pattern = '"[^"]*"'
  as.character(formula_for_simulation)[3] %>%
    stringr::str_extract_all(pattern = covars_pattern) %>%
    unlist %>%
    stringr::str_remove_all(pattern = '\"')
}



# Draw within-block connections
# @param seed_network a seed network from which a network will be simulated.
# @param formula_for_simulation formula for simulating a network
# @param coef_within_block a vector of within-block parameters. The order of the parameters should match that of the formula.
# @param ergm_control auxiliary function as user interface for fine-tuning ERGM simulation
# @param output Normally character, one of "network" (default), "stats", "edgelist", to determine the output format.
# @param seed seed value (integer) for the random number generator.
# @param nsim Number of networks to be randomly drawn from the given distribution on the set of all networks.
# @param verbose If this is TRUE/1, the program will print out additionalinformation about the progress of simulation.
# @param ... Additional arguments, to be passed to lower-level functions
# @export
# @return Simulated within-block connections, the output form depends on the parameter \code{output}.
draw_within_block_connection <- function(seed_network,
                                         formula_for_simulation,
                                         coef_within_block,
                                         ergm_control,
                                         output = "network",
                                         seed,
                                         nsim,
                                         verbose,
                                         ...) {
  # Extract the RHS of the formula
  # formula_rhs <- as.character(formula_for_simulation)[3]
  # This part is just important if the network was already blocked before 
  # attributes_tmp <- network::list.vertex.attributes(seed_network)
  # for(i in c(".NetworkID", ".NetworkName")){
  #   if(i %in% attributes_tmp){
  #     delete.vertex.attribute(seed_network, i)
  #   }
  # }
  # Split network up into block
  if("combined_networks" %in% class(seed_network)){
    seed_network_blocked <- seed_network
  } else {
    seed_network_blocked <- get_within_networks(seed_network, seed_network%v%"block")
  }
  # Simulate within-block links 
  # The suppressWarnings is necessary because the simulate_formula function throws a warning if some of the ... parameters are not used.
  suppressWarnings(
    within_conn <- ergm::simulate_formula(
      object = formula_for_simulation,
      basis = seed_network_blocked,
      nsim = nsim,
      coef = coef_within_block,
      seed = seed,
      control = ergm_control,
      output = output,
      verbose = verbose, 
      ...
    )  
  )
  
  # Cut out empty networks (we might have to change this at some point)
  if((nsim>1)& (output == "edgelist")){
    within_conn <- within_conn[unlist(lapply(within_conn, function(x) !is.null(dim(x))))]  
  }

  
  if(nsim>1){
    # Match the vertex.names of the seed network to the blocked network
    matching_tmp <- match(seed_network_blocked%v% "vertex.names", 1:length(seed_network_blocked%v% "vertex.names"))
    if (output == "edgelist") {
      within_conn <- lapply(within_conn, function(x) {
        x$.head <- matching_tmp[x$.head]
        x$.tail <- matching_tmp[x$.tail]
        attr(x,"vnames") <- seed_network%v%"vertex.names"
        return(x)
      })
    }
    } else {
    # Match the vertex.names of the seed network to the blocked network
    matching_tmp <- match(seed_network_blocked%v% "vertex.names", 1:length(seed_network_blocked%v% "vertex.names"))
    if (output == "edgelist") {
      within_conn$.head <- matching_tmp[within_conn$.head]
      within_conn$.tail <- matching_tmp[within_conn$.tail]
      attr(within_conn,"vnames") <- seed_network%v%"vertex.names"
    }
  }
  # Return the output
  return(within_conn)
}



# Draw between-block connections.
# @description 
draw_between_block_connection <- function(network, 
                                          formula_for_simulation,
                                          coef_between_block,
                                          seed = NULL,
                                          nsim = 1,
                                          verbose = 0,
                                          ergm_control = ergm::control.simulate.formula(),
                                          draw_networks = TRUE, 
                                          ...) {
  # Set seed if it is not provided
  if(is.null(seed)){
    seed <- sample.int(10^5, 1)
  }
  set.seed(seed)
  # Number of nodes
  N_node <- network$gal$n
  preprocessed_features <- get_features(network, formula_for_simulation)
  dot_list <- list(...)
  
  if(draw_networks){
    # This just changes the type of sparse matrices from "nsCMatrix" to "dsCMatrix"
    preprocessed_features$list_sparse_feature_adjmat <- lapply(preprocessed_features$list_sparse_feature_adjmat, function(x) x*1)
    sim_between <- function(seed) {
      # Simulate one between-block network. The output format is a sparse matrix.
      if("add_intercepts" %in% names(dot_list) && "clustering_with_features" %in% names(dot_list)){
        if(dot_list$add_intercepts && dot_list$clustering_with_features && dot_list$simulate_sbm){
          between_conn <- simulate_between_network_covariates(numOfVertices = N_node,
                                                              list_feature_adjmat = preprocessed_features$list_sparse_feature_adjmat,
                                                              block_membership = as.numeric(network %v% "block"), 
                                                              coef_between = lapply(dot_list$sbm_pi, function(x) as(Class = "dgCMatrix", object = x)),
                                                              directed = network$gal$directed, seed = seed)
        } else if(dot_list$add_intercepts && !dot_list$clustering_with_features && dot_list$simulate_sbm){
          between_conn <- simulate_between_network_no_covariates(numOfVertices = N_node,
                                                                 block_membership = as.numeric(network %v% "block"), 
                                                                 coef_between = as(Class = "dgCMatrix", object = dot_list$sbm_pi),
                                                                 directed = network$gal$directed, seed = seed)
        } else {
          coef_between_block <- check_edge_term(coef_between_block)
          between_conn <- simulate_between_network(
            numOfVertices = N_node, 
            list_feature_adjmat = preprocessed_features$list_sparse_feature_adjmat,
            coef_between = coef_between_block,
            block_membership = network %v% "block",
            directed =  network$gal$directed, seed = seed
          )
        }
      } else {
        coef_between_block <- check_edge_term(coef_between_block)
        # Simulate one between-block network. The output format is a sparse matrix.
        between_conn <- simulate_between_network(
          numOfVertices = N_node, 
          list_feature_adjmat = preprocessed_features$list_sparse_feature_adjmat,
          coef_between = coef_between_block,
          block_membership = network %v% "block",
          directed =  network$gal$directed, seed = seed
        )
      }
      
      # Convert the sparse matrix into an edgelist.
      between_conn <-
        as.data.frame(Matrix::summary(between_conn)) %>%
        dplyr::select("i", "j") %>%
        dplyr::rename(X1 = "i", X2 = "j")
      
      between_conn <-  as.matrix(between_conn)
      attr(between_conn, "n") <- N_node
      attr(between_conn, "vnames") <- network %v% "vertex.names"
      attr(between_conn, "directed") <- network$gal$directed
      attr(between_conn, "bipartite") <- FALSE
      attr(between_conn, "loops") <- FALSE
      attr(between_conn, "class") <- c("edgelist", "matrix")
      # Return the output
      return(between_conn)
    }
    
    
    # If you simulate just one between-block network:
    if (nsim == 1) {
      output <- sim_between(seed)
      if (verbose > 0) {
        message("Simulated one between-block network.")
      }
      # Return the output
      return(list(output = output, node_data = preprocessed_features$node_data))
    }
    
    # If you simulate more than one between-block network:
    if (nsim > 1) {
      i = 0
      output <-
        foreach(i = 1:nsim) %do% {
          net <- sim_between(seed*i)
          if (verbose > 0) {
            message(glue::glue("Finished between-block network simulation {i} of {nsim}."))
          }
          return(net)
        }
      if (verbose > 0) {
        message("Finished simulating between-block networks.")
      }
      return(list(output = output, node_data = preprocessed_features$node_data))
    }
  } else {
    if(nsim == 1){
      return(list(output = data.frame(tail = integer(0), head = integer(0)), 
                  node_data = preprocessed_features$node_data))  
    } else {
      return(list(output = rep(list( data.frame(tail = integer(0), head = integer(0))),nsim), 
                  node_data = preprocessed_features$node_data))
    }
    
    
    
  }
}


# Combine within- and between-block edges while removing duplicated links.
combine_within_between_edges <- function(edgelist_within, edgelist_between, use_fast_between_simulation, old_vertex_names) {
  # Store network information necessary to reconstruct a network object later.
  n <- attr(edgelist_within, "n")
  directed <- attr(edgelist_within, "directed")
  bipartite <- attr(edgelist_within, "directed")
  loops <- attr(edgelist_within, "loops")

  # Convert the edgelists into data frames
  edgelist_within <- data.frame(edgelist_within)
  colnames(edgelist_within) <- c("tail", "head")
  if(length(edgelist_between) == 0){
    edgelist_between <- data.frame(tail = integer(0), head = integer(0))
  } else {
    edgelist_between <- data.frame(edgelist_between)
  }
  edgelist_between <- data.frame(edgelist_between)
  colnames(edgelist_between) <- c("tail", "head")

  # From df_between, remove between-block edges that also appear in df_within.
  # To do so, we use dplyr::anti_join.
  # If use_fast_between_simulation, skip this step.
  if (!use_fast_between_simulation) {
    edgelist_between <-
      dplyr::anti_join(edgelist_between, edgelist_within, by = c("tail", "head"))
  }

  # Bind within- and between-block edges
  output <-
    rbind(edgelist_within, edgelist_between) %>%
    dplyr::arrange(tail) %>%
    as.matrix()

  # Attach network info
  attr(output, "n") <- n
  attr(output, "vnames") <- old_vertex_names
  attr(output, "directed") <- directed
  attr(output, "bipartite") <- bipartite
  attr(output, "loops") <- loops
  class(output) <- c("matrix_edgelist", "edgelist", class(output))

  # Return the output
  return(output)
}


# Create a final network
# @param sorted_dataframe a data frame generated by `sort_block_membership`
# @param edgelist an edgelist that contain both within- and between-block edges without duplication
# @param edgelist an edgelist that contain both within- and between-block edges without duplication
# @param old_vertex_names the original vertex names
#' @importFrom foreach foreach %do%
generate_network_for_output <- function(sorted_dataframe, 
                                        edgelist, old_vertex_names) {
  directed <- attr(edgelist, "directed")
  # Initialize a network object from the edgelist
  if(length(edgelist) == 0){
    edgelist <- data.frame(tail = 1, head = 2)
    g <- network::network(edgelist, directed = directed, 
                          matrix.type = "edgelist", vertex.attr = sorted_dataframe)
    delete.edges(g,1)
  } else{
    g <- network::network(edgelist, directed = directed, 
                          matrix.type = "edgelist", vertex.attr = sorted_dataframe)
  }
  # Restore the original vertex names
  g %v% "vertex.names" <- old_vertex_names
  return(g)
}



# Converts an edgelist into a matrix of sufficient network statistics
# @param net the net to extract the covariates from
# @param edgelist the edgelist
# @param between_formula the formula for the between connections
# @return a matrix of sufficient network statistics
edgelist_to_stats <- function(net, edgelist, between_formula) {
  g_copy <- network::network(edgelist, matrix.type = "edgelist", directed = FALSE, loops = FALSE)
  attrs <- network::list.vertex.attributes(net)
  for (attr in attrs) {
    network::set.vertex.attribute(g_copy, attr, network::get.vertex.attribute(net, attr))
  }
  form <- as.formula(paste("g_copy", "~", as.character(between_formula)[3]))
  summary(form)
}

# Converts a list of edgelists into a data frame of network statistics
# @param edgelists the list of edgelists
# @param between_formula the formula for the between connections
# @return a data frame of sufficient network statistics
get_between_stats <- function(edgelists, between_formula) {
  net <- get(as.character(between_formula)[2], envir = environment(between_formula))
  between_stats <-
    edgelists %>%
    purrr::map(function(el) {
      edgelist_to_stats(net, el, between_formula)
    }) %>%
    purrr::reduce(function(a, b) {
      rbind(a, b)
    })
  rownames(between_stats) <- NULL
  data.frame(between_stats)
}

#' @importFrom magrittr %>% %<>%

