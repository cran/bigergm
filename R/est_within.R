#' Estimate a within-block network model.
#' 
#' @description
#' Function to estimate the within-block model. Both pseudo-maximum likelihood and monte carlo approximate maximum likelihood estimators are implemented. 
#' 
#' @importFrom parallel mclapply
#' @importFrom magrittr %>%
#' @importFrom ergm.multi Networks
#' @importFrom ergm ergm
#' @importFrom utils combn
#' @importFrom stats binomial glm logLik
#' @importFrom rlang .data 
#' @importFrom foreach foreach %do% %dopar%
#' @param formula  An \R \code{\link{formula}} object of the form
#' \code{y ~ <model terms>}, where \code{y} is a
#' \code{\link[network]{network}} object. 
#' The network object must contain block information as a vertex attribute with the name 'block'.
#' For the details on the possible \code{<model terms>}, see
#' \code{\link[ergm]{ergmTerm}} and Morris, Handcock and Hunter (2008).
#' The \code{\link[ergm.multi]{L-ergmTerm}} is supported to enable size-dependent coefficients.
#' @param network a network object with one vertex attribute called 'block' representing which node belongs to which block
#' @param seed seed value (integer) for the random number generator
#' @param method If "MPLE" (the default), then the maximum pseudolikelihood estimator is returned.
#' If "MLE", then an approximate maximum likelihood estimator is returned.
#' @param add_intercepts Boolean value to indicate whether adequate intercepts 
#' should be added to the provided formula so that the model in the first stage 
#' of the estimation is a nested model of the estimated model in the second stage of the estimation
#' @param clustering_with_features Boolean value to indicate if the clustering 
#' was carried out making use of the covariates or not (only important if \code{add_intercepts = TRUE})
#' @param ... Additional arguments, to be passed to the \code{\link[ergm]{ergm}} function
#' @param return_network Boolean value to indicate if the network object should be returned in the output. 
#' This is needed if the user wants to use, e.g., the \code{\link[ergm]{gof}} function as opposed to the \code{\link{gof.bigergm}} function.
#' @importFrom rlang %||%
#' @references 
#' Morris M, Handcock MS, Hunter DR (2008). Specification of Exponential-Family Random Graph Models: Terms and Computational Aspects.
#' Journal of Statistical Software, 24.
#' @return 'ergm' object of the estimated model.
#' @examples
#' adj <- c(
#' c(0, 1, 0, 0, 1, 0),
#' c(1, 0, 1, 0, 0, 1),
#' c(0, 1, 0, 1, 1, 0),
#' c(0, 0, 1, 0, 1, 1),
#' c(1, 0, 1, 1, 0, 1),
#' c(0, 1, 0, 1, 1, 0)
#' )
#' adj <- matrix(data = adj, nrow = 6, ncol = 6)
#' rownames(adj) <- as.character(1001:1006)
#' colnames(adj) <- as.character(1001:1006)
#' 
#' # Use non-consecutive block names
#' block <- c(70, 70, 70, 70, 95, 95)
#' g <- network::network(adj, matrix.type = "adjacency", directed = FALSE)
#' g %v% "block" <- block
#' g %v% "vertex.names" <- 1:length(g %v% "vertex.names")
#' est <- est_within(
#' formula = g ~ edges,
#'   network = g,
#'   parallel = FALSE,
#'   verbose = 0,
#'   initial_estimate = NULL,
#'   seed = NULL,
#'   method = "MPLE", 
#'   add_intercepts = FALSE,
#'   clustering_with_features = FALSE
#'   )
#' @export
est_within <-
  function(formula,
           network,
           seed = NULL,
           method = "MPLE",
           add_intercepts = TRUE, 
           clustering_with_features = FALSE, 
           return_network = FALSE,
           ...) {
    varargs <-list(...)
    # %||% extracts the value on the left with a default value if null
    control <- varargs$control %||% ergm::control.ergm()
    offset.coef <- varargs$offset.coef %||% NULL
    block <- network %v% "block" 

    network <- get_within_networks(network, block)
    terms <- ergm::ergm_model(formula, nw = network)$terms
    terms <- terms[unlist(lapply(terms, function(x) !is.null(x$coef.names)))]
    varnames <-
      list_rhs.formula(formula) %>%
      as.character()
    dep_terms <-
      terms %>% purrr::map(function(t) {
        dep <- t$dependence
        is_dep <- is.null(dep) || dep || t$offset || grepl("N\\(",t$coef.names)
      }) %>% unlist()
    
    if(add_intercepts){
      if(clustering_with_features){
        paste_tmp <- paste(c("nodematch('block', diff = TRUE)",
                             paste0(varnames[!dep_terms],":nodematch('block', diff = TRUE)"),
                             varnames[dep_terms]), collapse = '+')
        within_formula <- as.formula(glue::glue("network ~ {paste_tmp}"))
      } else {
        within_formula <- as.formula(glue::glue("network ~ nodematch('block', diff = TRUE) + {paste(varnames, collapse = '+')}"))
      }
    } else {
      within_formula <- as.formula(glue::glue("network ~ {paste(varnames, collapse = '+')}"))
    }
    
    model_est <- ergm::ergm(formula = within_formula, estimate = method,
                            control = control,offset.coef = offset.coef)
    
    if(!return_network){
      # Remove unnecessary network objects
      model_est$newnetwork <- NULL
      model_est$network <- NULL
      model_est$within_formula <- within_formula  
    }
    return(model_est)
  }


