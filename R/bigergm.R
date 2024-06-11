#' bigergm: Exponential-family random graph models for large networks with local dependence
#' @description
#' The function \code{bigergm} estimates and simulates three classes of exponential-family 
#' random graph models for large networks under local dependence: 
#'\enumerate{
#'  \item{The p_1 model of Holland and Leinhardt (1981) in exponential-family form and extensions by Vu, Hunter, and Schweinberger (2013), Schweinberger, Petrescu-Prahova, and Vu (2014), Dahbura et al. (2021), and Fritz et al. (2024) to both directed and undirected random graphs with additional model terms, with and without covariates.}
#'  \item{The stochastic block model of Snijders and Nowicki (1997) and Nowicki and Snijders (2001) in exponential-family form.}
#'  \item{The exponential-family random graph models with local dependence of Schweinberger and Handcock (2015), with and without covariates. The exponential-family random graph models with local dependence replace the long-range dependence of conventional exponential-family random graph models by short-range dependence. Therefore, exponential-family random graph models with local dependence replace the strong dependence of conventional exponential-family random graph models by weak dependence, reducing the problem of model degeneracy (Handcock, 2003; Schweinberger, 2011) and improving goodness-of-fit (Schweinberger and Handcock, 2015). In addition, exponential-family random graph models with local dependence satisfy a weak form of self-consistency in the sense that these models are self-consistent under neighborhood sampling (Schweinberger and Handcock, 2015), which enables consistent estimation of neighborhood-dependent parameters (Schweinberger and Stewart, 2017; Schweinberger, 2017).}
#' }
#' @useDynLib bigergm
#' @importFrom utils tail
#' @importFrom cachem cache_disk
#' @importFrom tidyr expand 
#' @importFrom Rcpp sourceCpp
#' @importFrom stats simulate terms update
#' @importFrom ergm ergm.getnetwork gof control.ergm ergmMPLE
#' @importFrom graphics abline boxplot lines
#' @importFrom network `%v%` list.vertex.attributes get.vertex.attribute as.network delete.edges as.network
#' @importFrom network `%v%<-` delete.vertex.attribute network.initialize set.vertex.attribute
#' @param object An \R \code{\link{formula}} object or \code{\link{bigergm}} class object. 
#' If a formula is given, the function estimates a new model specified by it. 
#' It needs to be of the form
#' \code{y ~ <model terms>}, where \code{y} is a
#' \code{\link[network]{network}} object. 
#' For the details on the possible \code{<model terms>}, see
#' \code{\link[ergm]{ergmTerm}} and Morris, Handcock and Hunter (2008).
#' All terms that induce dependence are excluded from the between block model, while the within block model includes all terms.
#' When you pass a \code{\link{bigergm}} class object to the function, you continue from the previous MM step.
#' Note that the block allocation (which is either provided by parameter \code{blocks} or estimated in the first step) is saved as the vertex.attribute `block` of the network. 
#' This attribute can also be used in the specified formula.
#' The \code{\link[ergm.multi]{L-ergmTerm}} is supported to enable size-dependent coefficients for the within-blocks model. 
#' Note, however, that for size-dependent parameters of terms that are included in the between-blocks model, 
#' the intercept in the linear model provided to \code{\link[ergm.multi]{L-ergmTerm}} should not include the intercept. 
#' See the second example below for a demonstration. 
#' @param add_intercepts Boolean value to indicate whether adequate intercepts 
#' should be added to the provided formula so that the model in the first stage 
#' of the estimation is a nested model of the estimated model in the second stage of the estimation.
#' @param n_blocks The number of blocks. This must be specified by the user.
#' When you pass a \code{\link{bigergm}} class object to the function, you don't have to specify this argument.
#' @param n_cores The number of CPU cores to use.
#' @param blocks The pre-specified block memberships for each node.
#' If `NULL`, the latent community structure is estimated, assuming that the number of communities is `n_blocks`.
#' @param estimate_parameters If `TRUE`, both clustering and parameter estimation are implemented.
#' If `FALSE`, only clustering is executed.
#' @param verbose A logical or an integer: if this is TRUE/1,
#' the program will print out additional information about the progress of estimation and simulation.
#' A higher value yields lower level information.
#' @param n_MM_step_max The maximum number of MM iterations.
#' Currently, no early stopping criteria is introduced. Thus `n_MM_step_max` MM iterations are exactly implemented.
#' @param tol_MM_step Tolerance regarding the relative change of the lower bound of the likelihood 
#' used to decide on the convergence of the clustering step
#' @param initialization How the blocks should be initialized.
#' If `infomap` (the default), `igraph`' or Python's infomap is implemented.
#' If `random`, the initial clusters are randomly uniformally selected.
#' If `spectral`, spectral clustering is conducted.
#' If `walktrap`, the walktrap clustering algorithm as implemented in \code{\link[igraph]{cluster_walktrap}} is conducted.
#' If \code{initialization} 
#' is a vector of integers of the same length as the number of nodes in the provided network (in \code{object}),
#' then the provided vector is used as the initial cluster assignment.  
#' If \code{initialization} is a string relating to a file path, \code{\link{bigergm}} will interpret it as block allocations saved in Python's infomap .clu format under that path.
#' @param use_infomap_python If `TRUE`, the cluster initialization is implemented using Pythons' infomap.
#' @param virtualenv_python Which virtual environment should be used for the infomap algorithm? 
#' @param seed_infomap seed value (integer) for the infomap algorithm, which can be used to initialize the estimation of the blocks.
#' @param weight_for_initialization weight value used for cluster initialization. The higher this value, the more weight is put on the initialized block allocation.
#' @param seed seed value (integer) for the random number generator.
#' @param method_within If "MPLE" (the default), then the maximum pseudolikelihood estimator is implemented when estimating the within-block network model.
#' If "MLE", then an approximate maximum likelihood estimator is conducted. If "CD" (EXPERIMENTAL), the Monte-Carlo contrastive divergence estimate is returned. 
#' @param control_within A list of control parameters for the \code{\link[ergm]{ergm}} function used to estimate the parameters of the within model. See \code{\link[ergm]{control.ergm}} for details.
#' @param only_use_preprocessed If `TRUE`, the function only uses the preprocessed data from a previous fit but does not continue the estimation from its final iteration, instead the estimation is started again from the provided initialization.
#' @param clustering_with_features If `TRUE`, clustering is implemented using the discrete covariates specified in the formula.
#' @param compute_pi If `TRUE`, this function keeps track of pi matrices at each MM iteration.
#' If the network is large, we strongly recommend to set to be `FALSE`.
#' @param check_alpha_update If `TRUE`, this function keeps track of alpha matrices at each MM iteration.
#' If the network is large, we strongly recommend to set to be `FALSE`.
#' @param check_blocks If TRUE, this function keeps track of estimated block memberships at each MM iteration.
#' @param cache a `cachem` cache object used to store intermediate calculations such as eigenvector decomposition results.
#' @param return_checkpoint If `TRUE`, the function returns the checkpoint list. For most applications, this should be set to `TRUE` but if memory space needed by the output is an issue, set to `FALSE`.
#' @param ... Additional arguments, to be passed to lower-level functions (mainly to the \code{\link[ergm]{ergm}} function used for the estimation of within-block connections).
#' @references 
#' Babkin, S., Stewart, J., Long, X., and M. Schweinberger (2020). Large-scale estimation of random graph models with local dependence. Computational Statistics and Data Analysis, 152, 1--19.
#' 
#' Dahbura, J. N. M., Komatsu, S., Nishida, T. and Mele, A. (2021), ‘A structural model of business cards exchange networks’.
#' https://arxiv.org/abs/2105.12704
#' 
#' Fritz C., Georg C., Mele A., and Schweinberger M. (2024). A strategic model of software dependency networks.
#' https://arxiv.org/abs/2402.13375
#' 
#' Handcock, M. S. (2003). Assessing degeneracy in statistical models of social networks. Technical report, Center for Statistics and the Social Sciences, University of Washington, Seattle. \cr 
#' https://csss.uw.edu/Papers/wp39.pdf
#' 
#' Holland, P. W. and S. Leinhardt (1981). An exponential family of probability distributions for directed graphs. Journal of the American Statistical Association, Theory & Methods, 76, 33--65.
#' 
#' Morris M, Handcock MS, Hunter DR (2008). Specification of Exponential-Family Random Graph Models: Terms and Computational Aspects.
#' Journal of Statistical Software, 24.
#' 
#' Nowicki, K. and T. A. B. Snijders (2001). Estimation and prediction for stochastic blockstructures. Journal of the American Statistical Association, Theory & Methods, 96, 1077--1087.
#' 
#' Schweinberger, M. (2011). Instability, sensitivity, and degeneracy of discrete exponential families. Journal of the American Statistical Association, Theory & Methods, 106, 1361--1370.
#' 
#' Schweinberger, M. (2020). Consistent structure estimation of exponential-family random graph models with block structure. Bernoulli, 26, 1205--1233.
#' 
#' Schweinberger, M. and M. S. Handcock (2015). Local dependence in random graph models: characterization, properties, and statistical inference. Journal of the Royal Statistical Society, Series B (Statistical Methodology), 7, 647-676.
#' 
#' Schweinberger, M., Krivitsky, P. N., Butts, C.T. and J. Stewart (2020). Exponential-family models of random graphs: Inference in finite, super, and infinite population scenarios. Statistical Science, 35, 627-662.
#' 
#' Schweinberger, M. and P. Luna (2018). HERGM: Hierarchical exponential-family random graph models. Journal of Statistical Software, 85, 1--39.
#' 
#' Schweinberger, M., Petrescu-Prahova, M. and D. Q. Vu (2014). Disaster response on September 11, 2001 through the lens of statistical network analysis. Social Networks, 37, 42--55.
#' 
#' Schweinberger, M. and J. Stewart (2020). Concentration and consistency results for canonical and curved exponential-family random graphs. The Annals of Statistics, 48, 374--396.
#' 
#' Snijders, T. A. B. and K. Nowicki (1997). Estimation and prediction for stochastic blockmodels for graphs with latent block structure. Journal of Classification, 14, 75--100.
#' 
#' Stewart, J., Schweinberger, M., Bojanowski, M., and M. Morris (2019). Multilevel network data facilitate statistical inference for curved {ERGM}s with geometrically weighted terms. Social Networks, 59, 98--119.
#' 
#' Vu, D. Q., Hunter, D. R. and M. Schweinberger (2013). Model-based clustering of large networks. Annals of Applied Statistics, 7, 1010--1039.
#' @return An object of class 'bigergm' including the results of the fitted model. 
#' These include: 
#'\describe{
#'  \item{call:}{call of the mode}
#'  \item{block:}{vector of the found block of the nodes into cluster}
#'  \item{initial_block:}{vector of the initial block of the nodes into cluster}
#'  \item{sbm_pi:}{Connection probabilities represented as a \code{n_blocks x n_blocks} matrix from the first stage of the estimation between all clusters}
#'  \item{MM_list_z:}{list of cluster allocation for each node and each iteration}
#'  \item{MM_list_alpha:}{list of posterior distributions of cluster allocations for all nodes for each iteration}
#'  \item{MM_change_in_alpha:}{change in 'alpha' for each iteration}
#'  \item{MM_lower_bound:}{ vector of the evidence lower bounds from the MM algorithm}
#'  \item{alpha: }{matrix representing the converged posterior distributions of cluster allocations for all nodes}
#'  \item{counter_e_step:}{ integer number indicating the number of iterations carried out}
#'  \item{adjacency_matrix:}{sparse matrix representing the adjacency matrix used for the estimation}
#'  \item{estimation_status:}{character stating the status of the estimation}
#'  \item{est_within:}{\code{\link[ergm]{ergm}} object of the model for within cluster connections }
#'  \item{est_between:}{\code{\link[ergm]{ergm}} object of the model for between cluster connections}
#'  \item{checkpoint:}{list of information to continue the estimation (only returned if \code{return_checkpoint = TRUE})}
#'  \item{membership_before_kmeans:}{vector of the found blocks of the nodes into cluster before the final check for bad clusters}
#'  \item{estimate_parameters:}{binary value if the parameters in the second step of the algorithm should be estimated or not}
#' }
#' @examples
#' # Load an embedded network object.
#' data(toyNet)
#'
#' # Specify the model that you would like to estimate.
#' model_formula <- toyNet ~ edges + nodematch("x") + nodematch("y") + triangle
#' # Estimate the model
#' bigergm_res <- bigergm(
#'   object = model_formula,
#'   # The model you would like to estimate
#'   n_blocks = 4,
#'   # The number of blocks
#'   n_MM_step_max = 10,
#'   # The maximum number of MM algorithm steps
#'   estimate_parameters = TRUE,
#'   # Perform parameter estimation after the block recovery step
#'   clustering_with_features = TRUE,
#'   # Indicate that clustering must take into account nodematch on characteristics
#'   check_blocks = FALSE)
#'   
#'  # Example with N() operator
#'  
#'  \dontrun{
#' set.seed(1)
#' # Prepare ingredients for simulating a network
#' N <- 500
#' K <- 10
#' 
#' list_within_params <- c(1, 2, 2,-0.5)
#' list_between_params <- c(-8, 0.5, -0.5)
#' formula <- g ~ edges + nodematch("x") + nodematch("y")  + N(~edges,~log(n)-1)
#' 
#' memb <- sample(1:K,prob = c(0.1,0.2,0.05,0.05,0.10,0.1,0.1,0.1,0.1,0.1), 
#'                size = N, replace = TRUE)
#' vertex_id <- as.character(11:(11 + N - 1))
#' 
#' x <- sample(1:2, size = N, replace = TRUE)
#' y <- sample(1:2, size = N, replace = TRUE)
#' 
#' 
#' df <- tibble::tibble(
#'   id = vertex_id,
#'   memb = memb,
#'   x = x,
#'   y = y
#' )
#' g <- network::network.initialize(n = N, directed = FALSE)
#' g %v% "vertex.names" <- df$id
#' g %v% "block" <- df$memb
#' g %v% "x" <- df$x
#' g %v% "y" <- df$y
#' 
#' # Simulate a network
#' g_sim <-
#'   simulate_bigergm(
#'     formula = formula,
#'     coef_within = list_within_params,
#'     coef_between = list_between_params,
#'     nsim = 1, 
#'     control_within = control.simulate.formula(MCMC.burnin = 200000))
#' 
#' estimation <- bigergm(update(formula,new = g_sim~.), n_blocks = 10, 
#'                       verbose = T)
#' summary(estimation)
#' }
#' @export
bigergm <- function(object,
                    add_intercepts = FALSE, 
                    n_blocks = NULL,
                    n_cores = 1,
                    blocks = NULL,
                    estimate_parameters = TRUE,
                    verbose = 0,
                    n_MM_step_max = 100,
                    tol_MM_step = 0.0001,
                    initialization = "infomap",
                    use_infomap_python = FALSE, 
                    virtualenv_python = "r-bigergm",
                    seed_infomap = NULL,
                    weight_for_initialization = 1000,
                    seed = NULL,
                    method_within = "MPLE",
                    control_within = ergm::control.ergm(),
                    clustering_with_features = TRUE,
                    compute_pi = FALSE,
                    check_alpha_update = FALSE,
                    check_blocks = FALSE,
                    cache = NULL,
                    return_checkpoint = TRUE, 
                    only_use_preprocessed = FALSE,
                    ...) {
  ###################################################################################
  ###### Preparations for estimation ################################################
  ###################################################################################
  
  # Add a check to look for the indicator whether the estimation should not be continued but 
  # only the preprocessed data be used 
  if ("bigergm" %in% class(object)) {
    if(only_use_preprocessed){
      n_blocks <- object$checkpoint$n_blocks
      preprocessed_data <- object$checkpoint$preprocessed_data
      object <- object$checkpoint$formula
    }
  }
  # When the given object is formula:
  if ("formula" %in% class(object)) {
    # Check if initialization is a vector of block memberships, a .clu file to be read or whether blocks are known 
    check <- all(stringr::str_detect(initialization, ".clu")) + (length(initialization) >1) + !is.null(blocks)
    # If n_cluster is missing and blocks and initialized_cluster_data are NULL, stop the process.
    if (!check && is.null(n_blocks)){
      stop("\nThe argument 'n_blocks' is missing. Please specify the number of clusters.")
    } 
    
    if(length(initialization) >1){
      # If initialization is a vector of block memberships, use it as the initial block memberships but leave the n_blocks argument if provided as is.
      if (is.null(n_blocks)){
        n_blocks <- length(unique(initialization))  
      }
    } else if (all(stringr::str_detect(initialization, ".clu"))){
      initialized_cluster_data <- readr::read_delim(initialization, delim = " ", skip = 9, col_names = c("node_id", "block", "flow"), col_types = "iid")
      initialization <- as.numeric(dplyr::arrange(initialized_cluster_data, by_group = .data$node_id)$block)
      n_blocks <- length(unique(initialization))
    } 
    if (!is.null(blocks)){
      if (!is.null(n_blocks)){
        warning("The argument 'blocks' is given. The argument 'n_blocks' will be ignored.")
      }
      n_blocks <- length(unique(blocks))
      # Have the blocks be a numeric vector
      blocks <- as.numeric(factor(blocks))
    }
    
    
    # Get the formula
    formula <- object
    
    # Get network object from formula
    network <- ergm::ergm.getnetwork(formula)
    
    if(clustering_with_features) {
      # Check if clustering_with_features is TRUE without any features
      list_varname <- as.character(formula)[3]
      list_varname <- unlist(stringr::str_extract_all(string = list_varname, pattern = '"[^"]*"'))
      list_varname <- stringr::str_remove_all(string = list_varname, pattern = '\"')
      if(length(list_varname) == 0){
        warning("Parameter clustering_with_features was set to TRUE although no covariates are needed for the given formula")
        clustering_with_features <- FALSE
      }  
    }
    if(!exists("preprocessed_data")){
      preprocessed_data <- NULL
    }
    # If preprocessed_data is not NULL, check if the order of features names is the same with that of formula.
    if (!is.null(preprocessed_data)) {
      # If clustering_with_features = FALSE, stop the process.
      if (!clustering_with_features) {
        stop("\nSet clustering_with_features = TRUE since you are going to use vertex features for cluster estimation.")
      }
    }
  }
  
  MM_restart_object <- NULL
  # When the given object has a class of "bigergm":
  if ("bigergm" %in% class(object)) {
    # Inherit the following arguments from the previous estimation.
    # These arguments are so important that they won't be replaced by the arguments given by the user.
    network <- object$checkpoint$network
    formula <- object$checkpoint$formula
    n_blocks <- object$checkpoint$n_blocks
    clustering_with_features <- object$checkpoint$clustering_with_features
    add_intercepts <- object$checkpoint$add_intercepts
    preprocessed_data <- object$checkpoint$preprocessed_data
    # Inherit the following arguments from the previous estimation if not given by the user.
    # If given by the user, discard the inherited argument and use the given one.
    vec_arguments <-
      c(
        "n_cores",
        "estimate_parameters",
        "verbose",
        "n_MM_step_max",
        "seed",
        "method_within",
        "compute_pi",
        "check_alpha_update",
        "check_blocks"
      )
    
    message("Arguments:")
    message(glue::glue("Number of clusters = {n_blocks}"))
    for (i in 1:length(vec_arguments)) {
      if (do.call(missing, list(vec_arguments[[i]])) == TRUE) {
        assign(vec_arguments[[i]], object$checkpoint[[vec_arguments[[i]]]])
      }
      message(glue::glue("{vec_arguments[[i]]} = {get(vec_arguments[[i]])}"))
    }
    
    # Prepare an object to restart the MM with.
    MM_restart_object <-
      list(
        alpha = object$alpha,
        list_alpha = object$MM_list_alpha,
        list_z = object$MM_list_z,
        z_memb_init = object$z_memb_final_before_kmeans,
        change_in_alpha = object$MM_change_in_alpha,
        lower_bound = object$MM_lower_bound,
        counter_e_step = object$counter_e_step,
        adjacency_matrix = object$adjacency_matrix
      )
  }
  
  
  # For the entire estimation process the vertex.names are changed to be integers running from 1 to n 
  old_vertex_names <- network::network.vertex.names(network)
  network::network.vertex.names(network) <- 1:length(old_vertex_names)
  
  # If the blockss of each node are specified in the variable 'blocks' as integers between 1 and n_blocks,
  # the specified block memberships are fixed.
  all_blockss_fixed <- ifelse(is.null(blocks), FALSE, TRUE)
  # Make sure that if all_blockss_fixed == TRUE, estimate_parameters must also be TRUE.
  estimate_parameters <- ifelse(all_blockss_fixed, TRUE, estimate_parameters)
  
  ###################################################################################
  ###### First step: Estimating block memberships ###################################
  ###################################################################################
  
  # Estimate the memberships if they are not specified.
  if (!all_blockss_fixed) {
    set.seed(seed)
    # Estimate the block memberships
    answer <- MM_wrapper(
      network = network,
      formula = formula,
      n_blocks = n_blocks,
      n_MM_step_max = n_MM_step_max,
      tol_MM_step = tol_MM_step, 
      min_size = 2,
      initialization = initialization,
      use_infomap_python = use_infomap_python,
      virtualenv_python = virtualenv_python, 
      seed_infomap = seed_infomap,
      clustering_with_features = clustering_with_features,
      list_multiplied_feature_matrices = preprocessed_data,
      verbose = verbose,
      weight_for_initialization = weight_for_initialization,
      compute_pi = compute_pi,
      check_alpha_update = check_alpha_update,
      check_blocks = check_blocks,
      MM_restart_object = MM_restart_object,
      cache = cache
    )
    
    blocks <- answer$z_memb_final
    initial_block <- answer$z_memb_init
    membership_before_kmeans <- answer$z_memb_final_before_kmeans
    sbm_pi <- answer$Pi
    MM_list_alpha <- answer$list_alpha
    MM_list_z <- answer$list_z
    MM_change_in_alpha <- answer$change_in_alpha
    MM_lower_bound <- answer$lower_bound
    alpha <- answer$alpha
    counter_e_step <- answer$counter_e_step
    adjacency_matrix <- as(answer$adjacency_matrix, "nMatrix")
    preprocessed_data <- answer$list_multiplied_feature_matrices
  }
  else {
    if (verbose > 0) {
      message("\nSkipping Steps 1 and 2: z specified")
    }
    param_MM_wrapper <- NULL
    initial_block <- NULL
    membership_before_kmeans <- NULL
    sbm_pi <- NULL
    MM_list_z <- NULL
    MM_list_alpha <- NULL
    MM_change_in_alpha <- NULL
    MM_lower_bound <- NULL
    alpha <- NULL
    counter_e_step <- NULL
    adjacency_matrix <- NULL
    preprocessed_data <- NULL
  }
  # Relabel the block to be in canonical form
  blocks <- relabel(block = blocks, n_blocks = n_blocks)
  ####################################################################################################
  ###### Second step: Estimating between-block parameters ############################################
  ####################################################################################################
  network %v% "block" <- as.numeric(blocks)
  if (estimate_parameters) {
    if (verbose > 0) {
      message("\nEstimate between-block parameters")
    }
    ## Estimate between-block parameters
    est_between <- est_between(
      formula = formula,
      network = network,
      add_intercepts = add_intercepts, 
      clustering_with_features = clustering_with_features
    )
  } else {
    # If you don't estimate any parameters...
    est_between <- NULL
  }
  
  ####################################################################################################
  ###### Third step: Estimating within-block parameters ##############################################
  ####################################################################################################
  
  # When you estimate with-block parameters, then:
  if (estimate_parameters) {
    if (verbose > 0) {
      message("\n\nStep 3: Estimate parameters conditional on z")
    }
    # Estimate within-block parameters
    if (verbose == 0) {
      suppressMessages(est_within <- est_within(formula =formula,network = network,
                                                block = blocks,number_cores = n_cores,
                                                seed = NULL,method = method_within,
                                                add_intercepts = add_intercepts,
                                                clustering_with_features = clustering_with_features,
                                                control = control_within,
                                                ...))
    } else {
      est_within <- est_within(formula =formula,network = network,
                               block = blocks,number_cores = n_cores,
                               seed = NULL,method = method_within,
                               add_intercepts = add_intercepts,
                               clustering_with_features = clustering_with_features,
                               control = control_within,
                               ...)
    }
    
    
    
    ####################################################################################################
    ###### Store the results ###########################################################################
    ####################################################################################################
    
    estimation_status <- ifelse(est_within$failure, "Estimation failed", "Estimated")
  } else {
    # estimate_parameters = FALSE:
    message("\n")
    estimation_status <- "Not estimated"
    est_within <- NULL
  }
  
  if(return_checkpoint){
    # Change the names back to the old ones
    network%v% "vertex.names" <- old_vertex_names
    # Store the given arguments. These will be inherited for the next MM.
    checkpoint <- list(
      network = network,
      formula = formula,
      n_blocks = n_blocks,
      clustering_with_features = clustering_with_features,
      preprocessed_data = preprocessed_data,
      n_cores = n_cores,
      estimate_parameters = estimate_parameters,
      verbose = verbose,
      n_MM_step_max = n_MM_step_max,
      seed = seed,
      method_within = method_within,
      compute_pi = compute_pi,
      check_alpha_update = check_alpha_update,
      check_blocks = check_blocks, 
      add_intercepts = add_intercepts
    )
  } else{
    checkpoint <- NULL
  }
  
  # Store the results in a list
  output <- list(
    call = sys.calls()[[1]], 
    formula = formula,
    block = blocks,
    initial_block = initial_block,
    sbm_pi = sbm_pi,
    MM_list_z = MM_list_z,
    MM_list_alpha = MM_list_alpha,
    MM_change_in_alpha = MM_change_in_alpha,
    MM_lower_bound = MM_lower_bound,
    alpha = alpha,
    counter_e_step = counter_e_step,
    adjacency_matrix = adjacency_matrix,
    directed = network$gal$directed,
    estimation_status = estimation_status,
    est_within = est_within,
    est_between = est_between,
    checkpoint = checkpoint,
    membership_before_kmeans = membership_before_kmeans, 
    estimate_parameters = estimate_parameters
  )
  
  # Return the output
  return(structure(output, class = "bigergm"))
}


