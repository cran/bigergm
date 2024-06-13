# Get a list of sparse feature matrices from a formula and network. 
# @description
# Get a list of sparse feature adjacency matrix from a formula. 
# These matrices can be given to the  \code{\link{simulate_bigergm}} and \code{\link{gof}} functions as parameters. 
# The respective parameter is only needed for both functions if 
# \code{use_fast_between_simulation} is set to \code{TRUE}.
# Generally, this function should only be used if users are working with large networks.
#  
#' @importFrom foreach foreach %do%
# @param network a network object from which nodal covariates are extracted.
# @param formula a network model to be considered
# @return The list of sparse matrices of feature matrices that are used for the first step of the estimation.
# @examples
# data(toyNet)
# model_formula <- toyNet ~ edges + nodematch("x") + nodematch("y") 
# preprocessed_features <- 
#   get_features(toyNet, model_formula)
# @export
get_features <- function(network, formula) {
  term_names <- rownames(attr(terms(formula),"factors"))
  # Get variable names from formula (extract strings sandwiched by double quotes)
  if(sum(grepl("N\\(",term_names))>0){
    formula <- update(formula, new = paste("~. -",paste(term_names[grepl("N\\(",term_names)], collapse = "-")))   
  }
  # Exclude the offset terms
  formula <- no.offset(formula)
  
  ergm_model_info <- ergm::ergm_model(formula, nw = network)
  # Extract nodematch terms
  nodematch_terms <- unlist(lapply(ergm_model_info$terms, function(x)x$name == "nodematch"))
  i <- 0
  
  if(sum(nodematch_terms)>0){
    list_var_names <- unlist(lapply(ergm_model_info$terms[nodematch_terms], 
                                     function(x) gsub(pattern = "nodematch.", 
                                                      replacement = "", 
                                                      x$coef.names)))
    
    list_var <-
      foreach(i = 1:length(list_var_names)) %do% {
        feature <- network::get.vertex.attribute(x = network, list_var_names[i])
        return(feature)
      }
    
    tmp_data <- data.frame(list_var)
    node_data <- data.frame("vertex.names" = network %v% "vertex.names",
                            "block" = network %v% "block",
                            tmp_data)
    names(node_data)[-c(1,2)] <-  list_var_names
    # Create a list of sparse feature adjacency matrices
    list_sparse_feature_adjmat <-
      foreach(i = 1:length(list_var)) %do% {
        if (is.numeric(list_var[[i]])) {
          output <- get_sparse_feature_adjmat_fast(list_var[[i]])
        } else {
          output <- get_sparse_feature_adjmat_fast(list_var[[i]])
        }
      }

    # Attach variable names to each matrix
    names(list_sparse_feature_adjmat) <- names(list_var)
    
    # Return the output
    return(list(list_sparse_feature_adjmat = list_sparse_feature_adjmat,
                node_data = node_data))
  }else{
    return(list(node_data = data.frame("vertex.names" = network %v% "vertex.names",
                                           "block" = network %v% "block")))
  }
 
 
}

get_sparse_feature_adjmat_fast <- function(x){
  tmp <- table(x)
  outcomes <- names(tmp)[tmp>1]
  indicators <-lapply(outcomes, FUN = function(y, x) t(combn(which(x == y),2)), x = x)
  indicators <- do.call(rbind, indicators)
  return(sparseMatrix(indicators[,1], indicators[,2], symmetric = TRUE, dims = c(length(x),length(x))))
}

decimal_to_binary_vector_R <- function(decimal, vec_length) {
  n <- decimal
  output <- numeric(vec_length)
  
  for (i in 1:vec_length) {
    output[i] <- n %% 2
    n <- n %/% 2
  }
  return(as.vector(output))
}

get_matrix_for_denominator_R <- function(numOfVertices, list_feature_adjmat) {
  n_feature <- length(list_feature_adjmat)
  n_item <- 2 ^ n_feature
  res <- sparseMatrix(i = 1,j = 1, x = 0,dims = c(numOfVertices,numOfVertices))

  for (s in 1:(n_item - 1)) {
    index <- decimal_to_binary_vector_R(s, n_feature)
    k <- sum(index)
    counter <- 0
    for (t in 1:n_feature) {
      if (index[t] == 1) {
        counter <- counter + 1
        if (counter == 1) {
          X <- list_feature_adjmat[[t]]
        } else {
          X <- X*list_feature_adjmat[[t]]
        }
      }
      gc(full = TRUE)  
    }
    res <- res + (-1) ^ k * X  
    gc(full = TRUE)  
  }
  return(Matrix::drop0(res))
}

get_elementwise_multiplied_matrices_R <- function(adjmat, list_feature_adjmat,memory_efficient = FALSE) {
  n_node <- nrow(adjmat)
  n_feature <- length(list_feature_adjmat)
  list_mat <- c(adjmat, list_feature_adjmat)
  
  n_matrix <- length(list_mat)
  length_output <- 2 ^ n_matrix
  output <- vector("list", length = length_output)
  
  output[[1]] <- get_matrix_for_denominator_R(n_node, list_feature_adjmat)
  rm(adjmat, list_feature_adjmat)
  if(memory_efficient){
    gc(full = TRUE)  
  }
  for (s in 1:(length_output-1)) {
    index <- decimal_to_binary_vector_R(s, n_matrix)
    X <- sparseMatrix(i = 1,j = 1, x = 0,dims = c(n_node,n_node))
    counter <- 0
    for (t in 1:n_matrix) {
      if (index[t] == 1) {
        counter <- counter + 1
        if (counter == 1) {
          X <- list_mat[[t]]
        } else {
          X <- X * list_mat[[t]]
        }
      }
    }

    for (t in 1:n_matrix) {
      if (index[t] == 0) {
        X <- X - X*list_mat[[t]]
      }
    }
    output[[s+1]] <- as(Matrix::drop0(X),"nMatrix")
    if(memory_efficient){
      gc(full = TRUE) 
    }
  }
  return(output)
}


# Preprocess a list of sparse feature adjacency matrix from a formula. 
# @description
# These matrices can be given to the \code{\link{bigergm}} function as parameters. 
# Generally, this function should only be used if users are working with large networks and 
# are planning to estimate the model in multiple steps.
# Since the \code{\link{bigergm}} function is designed to continue the estimation from a previous point, 
# one can thereby save time by using this function and thus not repeat the preprocessing steps with each restart.
# @param net a network object from which nodal covariates are extracted.
# @param preprocessed_features a list of feature adjacency matrices generated by the \code{\link{get_features}} function.
# @export
# @return A list of sparse matrices of multiplied feature matrices that are needed for carrying our the first step of the estimation if the covariates should be used.
# @examples
# data(toyNet)
# \donttest{
# model_formula <- toyNet ~ edges + nodematch("x") + nodematch("y") 
# preprocessed_features <- get_features(toyNet, model_formula)
# multiplied_feature_matrices <- 
#   preprocess(net = toyNet,
#   preprocessed_features = preprocessed_features)
# }
preprocess <- function(net, preprocessed_features) {
  adj <- network::as.edgelist(net)
  N <- net$gal$n
  adj <- Matrix::sparseMatrix(i = adj[, 1], j = adj[, 2], x = 1, dims = c(N, N), symmetric = !net$gal$directed)
  get_elementwise_multiplied_matrices_R(adj, preprocessed_features)
}
