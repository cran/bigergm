#' Estimate between-block parameters
#' 
#' @description
#' Function to estimate the between-block model by relying on the maximum likelihood estimator. 
#' 
#' @importFrom magrittr %>%
#' @importFrom RcppArmadillo armadillo_set_seed
#' @importFrom ergm ergm
#' @param formula  An \R \code{\link{formula}} object of the form
#' \code{y ~ <model terms>}, where \code{y} is a
#' \code{\link[network]{network}} object. 
#' The network object must contain block information as a vertex attribute with the name 'block'.
#' For the details on the possible \code{<model terms>}, see
#' \code{\link[ergm]{ergmTerm}} and Morris, Handcock and Hunter (2008).
#' All terms that induce dependence are excluded from the between block model.
#' @param network a network object with one vertex attribute called 'block' representing which node belongs to which block
#' @param add_intercepts Boolean value to indicate whether adequate intercepts 
#' should be added to the provided formula so that the model in the first stage
#'  of the estimation is a nested model of the estimated model in the second stage of the estimation
#' @param clustering_with_features Boolean value to indicate if the clustering 
#' was carried out making use of the covariates or not (only important if \code{add_intercepts = TRUE})

#' @return 'ergm' object of the estimated model.
#' @references 
#' Morris M, Handcock MS, Hunter DR (2008). Specification of Exponential-Family Random Graph Models: Terms and Computational Aspects.
#' Journal of Statistical Software, 24.
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
#' block <- c(50, 70, 95, 50, 95, 70)
#' g <- network::network(adj, matrix.type = "adjacency")
#' g %v% "block" <- block
#' est <- est_between(
#'   formula = g ~ edges,network = g,
#'   add_intercepts = FALSE, clustering_with_features = FALSE
#' )
#' @export
#' @importFrom statnet.common list_rhs.formula 
#' @importFrom reticulate py_module_available
est_between <- function(formula, network,
                                   add_intercepts = TRUE,
                                   clustering_with_features = FALSE) {

  # Exclude all N() terms for the between-network model
  term_names <- rownames(attr(terms(formula),"factors"))
  if(sum(grepl("N\\(",term_names))>0){
    formula <- update(formula, new = paste("~. -",paste(term_names[grepl("N\\(",term_names)], collapse = "-")))   
  }
  # Exclude the offset terms
  formula <- no.offset(formula)
  # Create a formula that contains only dyad-independent terms. i.e. exclude externality terms like triangle.
  terms <- ergm::ergm_model(formula, nw = network)$terms
  block <- network%v% "block"
  varnames <-
    list_rhs.formula(formula) %>%
    as.character()
  dep_terms <-
    terms %>% purrr::map(function(t) {
      dep <- t$dependence
      is_dep <- is.null(dep) || dep
    }) %>% unlist()
  between_rhs <- varnames[!dep_terms]
  if(add_intercepts){
    tmp_mat <- matrix(NA, nrow = length(unique(block)), ncol = length(unique(block)))
    tmp_info <- generate_combinations_vectorized(length(unique(block)), directed = network$gal$directed)
    tmp_mat[tmp_info[[2]]] <- tmp_info[[1]]
    if(clustering_with_features){
      paste_tmp <- paste(paste0(between_rhs,":nodemix('block', levels2 = tmp_mat)"), collapse = '+')
      between_formula <- as.formula(glue::glue("network ~ nodemix('block', levels2 = tmp_mat) + {paste_tmp}"))
      between_formula_mple <- as.formula(glue::glue("network ~ nodematch('block') + nodemix('block', levels2 = tmp_mat) + {paste_tmp}"))
    } else {
      if(length(between_rhs)>0){
        between_formula <- as.formula(glue::glue("network ~ nodemix('block', levels2 = tmp_mat) + {paste(between_rhs, collapse = '+')}"))
        between_formula_mple <- as.formula(glue::glue("network ~ nodematch('block') + nodemix('block', levels2 = tmp_mat) + {paste(between_rhs, collapse = '+')}"))  
      } else {
        between_formula <- as.formula(glue::glue("network ~ nodemix('block', levels2 = tmp_mat)"))
        between_formula_mple <- as.formula(glue::glue("network ~ nodematch('block') + nodemix('block', levels2 = tmp_mat)"))
      }
      
    }
  } else {
    between_formula <- as.formula(glue::glue("network ~ {paste(between_rhs, collapse = '+')}"))
    between_formula_mple <- as.formula(glue::glue("network ~ nodematch('block') + {paste(between_rhs, collapse = '+')}"))
  }
  
  # Estimate logit
  between_logit <- mple_estimator(formula = between_formula_mple,
                                  formula_full =  between_formula,
                                  value =  0, control = control.ergm(), n_dependent = 0)

  # Return the output
  return(between_logit)
}


mple_estimator <- function(formula,formula_full, value, control, n_dependent){
  tmp <- ergmMPLE(formula)
  include <- which(tmp$predictor[,1] == value)
  tmp$response <- tmp$response[include]
  tmp$predictor <- tmp$predictor[include,-1]
  tmp$weights <- tmp$weights[include]
  tmp_data <- data.frame(cbind(tmp$response, tmp$predictor, tmp$weights))
  if(is.vector(tmp$predictor)){
    names(tmp_data) <- c("response",  formula[[3]][[3]], "weight")
  } else {
    names(tmp_data) <- c("response",colnames(tmp$predictor), "weight")
  }
  # the -1 in the formula is to remove the intercept and the -weight is to remove the weight from the formula
  glm_fit <- glm(response~. -1-weight, family = binomial(), data = tmp_data, weights = tmp_data$weight)
  # This is to estimate the null model, which only includes the intercept
  # the -weight is to remove the weight from the formula
  glm_fit_null <- glm(response~ 1-weight, family = binomial(), data = tmp_data, weights = tmp_data$weight)
  glm_summary <- summary(glm_fit)
  res <- list()
  res$coefficients <- glm_fit$coefficients
  res$iterations <- glm_fit$iter
  res$MCMCtheta <- glm_fit$coefficients
  res$gradient <- rep(NA,length(glm_fit$coefficients))
  res$hessian <- -solve(glm_summary$cov.unscaled) 
  res$covar <- glm_summary$cov.unscaled
  res$failure <- !glm_fit$converged
  res$mple.lik <- logLik(glm_fit)
  res$mple.lik.null <- logLik(glm_fit_null)
  res$mle.lik <- logLik(glm_fit)
  res$null.lik <- logLik(glm_fit_null)
  res$estimate <- "MPLE"
  res$control <- control
  res$ergm_version <- as.package_version("4.6.0")
  res$formula <- formula_full
  res$info <- list(terms_dind = FALSE, space_dind = TRUE, n_info_dyads = sum(tmp$weight), obs = FALSE,valued = FALSE, MPLE_is_MLE = n_dependent>0)
  res$etamap$offsettheta <- rep(FALSE,length(res$coefficients))
  res$between_formula <- formula
  class(res) <- "ergm"
  return(res)
}

