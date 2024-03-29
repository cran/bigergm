#' Get a list of sparse feature adjacency matrix from a formula
#' @importFrom foreach foreach %do%
#' @param network a network object from which nodal covariates are extracted.
#' @param formula a network model to be considered
#' @return The list of sparse matrices of feature matrices that are used for the first step of the estimation.
#' @examples
#' data(toyNet)
#' model_formula <- toyNet ~ edges + nodematch("x") + nodematch("y") 
#' list_feature_matrices <- 
#'   get_list_sparse_feature_adjmat(toyNet, model_formula)
#' @export
get_list_sparse_feature_adjmat <- function(network, formula) {
  # Get variable names from formula (extract strings sandwiched by double quotes)
  list_varname <- as.character(formula)[3]
  list_varname <- unlist(stringr::str_extract_all(string = list_varname, pattern = '"[^"]*"'))
  list_varname <- stringr::str_remove_all(string = list_varname, pattern = '\"')

  i <- 1
  # Extract variable from a network object
  list_var <-
    foreach(i = 1:length(list_varname)) %do% {
      feature <- network::get.vertex.attribute(x = network, list_varname[i])
      return(feature)
    }

  # Create a list of sparse feature adjacency matrices
  list_sparse_feature_adjmat <-
    foreach(i = 1:length(list_var)) %do% {
      if (is.numeric(list_var[[i]])) {
        output <- get_sparse_feature_adjmat(list_var[[i]])
      } else {
        output <- get_sparse_feature_adjmat_from_string(list_var[[i]])
      }
    }

  # Attach variable names to each matrix
  names(list_sparse_feature_adjmat) <- list_varname

  # Return the output
  return(list_sparse_feature_adjmat)
}

#' Get a list of sparse feature adjacency matrix from a formula. 
#' @description
#' These matrices can be given to the \code{\link{hergm}} function as parameters. Generally, this function should only be used if users are working with large networks and are planning to continually estimate the model.
#' @param net a network object from which nodal covariates are extracted.
#' @param list_feature_matrices a list of feature adjacency matrices generated by `get_list_sparse_feature_adjmat()`.
#' @export
#' @return A list of sparse matrices of multiplied feature matrices that are needed for carrying our the first step of the estimation if the covariates should be used.
#' @examples
#' data(toyNet)
#' \donttest{
#' model_formula <- toyNet ~ edges + nodematch("x") + nodematch("y") 
#' list_feature_matrices <- get_list_sparse_feature_adjmat(toyNet, model_formula)
#' multiplied_feature_matrices <- 
#'   compute_multiplied_feature_matrices(net = toyNet,
#'   list_feature_matrices = list_feature_matrices)
#' }
compute_multiplied_feature_matrices <- function(net, list_feature_matrices) {
  adj <- network::as.edgelist(net)
  N <- net$gal$n
  adj <- Matrix::sparseMatrix(i = adj[, 1], j = adj[, 2], x = 1, dims = c(N, N), symmetric = TRUE)
  get_elementwise_multiplied_matrices(adj, list_feature_matrices)
}
