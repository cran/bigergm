
#' @export
summary.bigergm <- function(object, ...) {
  cat("Call:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cat("Found", length(unique(object$block)), "clusters of relative sizes: \n", table(object$block)/length(object$block))
  
  if(object$estimate_parameters){
    cat("\n\nResults of within-cluster estimation: \n")
    tmp_est_within <- summary(object$est_within)
    print(tmp_est_within, print.message = FALSE, print.call = FALSE, print.deviance = FALSE)
    cat("\nResults of between-cluster estimation: \n")
    tmp_est_between <- summary(object$est_between)
    print(tmp_est_between, print.message = FALSE, print.call = FALSE, print.deviance = FALSE)
    cat("AIC:", tmp_est_within$aic + tmp_est_between$aic,"BIC: ",tmp_est_within$bic + tmp_est_between$bic)
  }
}


# Copy all vertex.attributes from net_tmp to net
# @param net a statnet network object
# @param net_tmp a statnet network object
copy_vertex_attributes <- function(net, net_tmp){
  attributes_tmp <- network::list.vertex.attributes(net_tmp)
  match_vertexnames <- match(net%v% "vertex.names",net_tmp%v% "vertex.names")
  for(i in attributes_tmp){
    set.vertex.attribute(net, i, network::get.vertex.attribute(net_tmp, i)[match_vertexnames])
  }
  return(net)
}

# Function to check if the first entry of a coef_between_block vector corresponds to the edge term and moves it forward if needed 
check_edge_term <- function(coef_between_block){
  # Check whether the first element of coef_between_block is "edges". If not, move it to the first position.
  if("edges" %in% names(coef_between_block)){
    position_edges <- which(names(coef_between_block) == "edges")
    if(position_edges != 1){
      coef_between_block <- c(coef_between_block[position_edges],
                              coef_between_block[-position_edges])
    }
  }
  return(coef_between_block)
}

# Function to turn a sparse matrix to an edgelist 
sparse_matrix_to_edgelist <- function(matrix, directed, vertex.names, N_node){
  # Convert the sparse matrix into an edgelist.
  matrix <-
    as.data.frame(Matrix::summary(matrix)) %>%
    dplyr::select("i", "j") %>%
    dplyr::rename(X1 = "i", X2 = "j")
  
  matrix <-  as.matrix(matrix)
  attr(matrix, "n") <- N_node
  attr(matrix, "vnames") <- vertex.names
  attr(matrix, "directed") <- directed
  attr(matrix, "bipartite") <- FALSE
  attr(matrix, "loops") <- FALSE
  attr(matrix, "class") <- c("edgelist", "matrix")
  # Return the output
  return(matrix)
}

# Function that deletes all offset terms from a formula
no.offset <- function(x) {
  k <- 0
  proc <- function(x) {
    if (length(x) == 1) return(x)
    if (x[[1]] == as.name("offset")) return(x[[1]])
    replace(x, -1, lapply(x[-1], proc))
  }
  update(proc(x), . ~ . - offset)
}

# Function used in est_within and draw_within_sample to generate the combinations of the block labels
generate_combinations_vectorized <- function(n, directed) {
  # Generate combinations of two numbers from 1 to n
  combinations <- t(combn(1:n, 2))
  if(directed){
    combinations <- rbind(combinations, cbind(combinations[,2], combinations[,1]))
  }
  # Convert the combinations to a character vector with an underscore separator
  combinationnames <- paste(combinations[,1], combinations[,2], sep = "_")
  # Return the vector of combinations
  return(list(combinationnames,combinations))
}

#' @export
print.bigergm <- function(x, ...) {
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cat("Found", length(unique(x$block)), "clusters of relative sizes: \n", 
      table(x$block)/length(x$block), "\n")
  
  if(x$estimate_parameters){
    cat("\nResults of within-cluster estimation: \n")
    print(x$est_within$coefficients)
    cat("\nResults of between-cluster estimation: \n")
    print(x$est_between$coefficients)
  }
}

# Check which nodes are in the same cluster
same_clusters <- function(z_memb) {
  n_nodes <- length(z_memb)
  block_sizes <- table(z_memb)
  memb_inds <- as.numeric(names(block_sizes))
  indicator = matrix('grey', n_nodes, n_nodes)
  for (i in memb_inds)
  {
    indicator[z_memb == i, z_memb == i] <- 'black'
  }
  indicator
}

#' Plot the network with the found clusters
#' @description
#' This function plots the network with the found clusters. The nodes are colored according to the found clusters.
#' Note that the function uses the `network` package for plotting the network and should therefore not be used for large networks with more than 1-2 K vertices
#' @param x The output of the bigergm function
#' @param ... Additional arguments, to be passed to lower-level functions
#' @importFrom network network
#' @importFrom stats var
#' @export
plot.bigergm <- function(x, ...) {
  plot(
    x$checkpoint$network,
    vertex.col = x$block,
    edge.col = same_clusters(x$block),
    main = ""
  )
}

# Function that projects the block labels to the canonical form, 
# which is the order of appearance of the blocks in the block vector.
relabel <- function(block, n_blocks) {
  levels(block) <-match(levels(block),block[which(!duplicated(block))])
  return(factor(block, levels = 1:n_blocks))
}


#' Install optional Python dependencies for bigergm
#' @description
#' Install Python dependencies needed for using the Python implementation of infomap. 
#' The code uses the `reticulate` package to install the Python packages `infomap` and `numpy`.
#' These packages are needed for the `bigergm` function when `use_infomap_python = TRUE` else the Python implementation is not needed.
#' 
#' @param envname The name, or full path, of the environment in which Python packages 
#' are to be installed. When NULL (the default), the active environment as set 
#' by the RETICULATE_PYTHON_ENV variable will be used; if that is unset, then the 
#' r-reticulate environment will be used.
#' @param method Installation method. By default, 
#' "auto" automatically finds a method that will
#'  work in the local environment. Change the 
#'  default to force a specific installation method. 
#'  Note that the "virtualenv" method is not available on Windows.
#' @param ... Additional arguments, to be passed to lower-level functions
#' @return No return value, called for installing the Python dependencies 'infomap' and 'numpy' 
#' @export
py_dep <- function(envname = "r-bigergm", method = "auto", ...) {
  reticulate::py_install("infomap", envname = envname, 
                         method = method, pip = TRUE, ...)
  reticulate::py_install("numpy", envname = envname, 
                         method = method, pip = TRUE, ...)
}

#' Compute the adjusted rand index (ARI) between two clusterings
#' @description
#' This function computes the adjusted rand index (ARI) of the true and estimated block membership (its definition can be found here \url{https://en.wikipedia.org/wiki/Rand_index}).
#' The adjusted rand index is used as a measure of association between two group membership vectors. 
#' The more similar the two partitions z_star and z are, the closer the ARI is to 1.
#' @param z_star The true block membership
#' @param z The estimated block membership
#' @return The adjusted rand index
#' @export
#' @examples
#' data(toyNet)
#' set.seed(123)
#' ari(z_star = toyNet%v% "block",
#' z = sample(c(1:4),size = 200,replace = TRUE))
ari <- function(z_star, z) {
  # Compute the contingency table
  cont_table <- table(factor(z_star, levels = unique(z_star)), 
                      factor(z, levels = unique(z)))
  
  # Sum of squares of sums of the contingency table
  sum_comb_c <- sum(choose(colSums(cont_table), 2))
  sum_comb_k <- sum(choose(rowSums(cont_table), 2))
  sum_comb <- sum(choose(cont_table, 2))
  
  n <- sum(cont_table)
  total_comb <- choose(n, 2)
  
  # Expected index and max index
  expected_index <- (sum_comb_k * sum_comb_c) / total_comb
  max_index <- (sum_comb_k + sum_comb_c) / 2
  
  # Adjusted Rand Index
  ari <- (sum_comb - expected_index) / (max_index - expected_index)
  return(ari)
}

#' Obtain the within-block networks defined by the block attribute.
#' 
#' @description
#' Function to return a list of networks, each network representing the within-block network of a block.
#' 
#' @param network a network object
#' @param block a vector of integers representing the block of each node
#' @param combined_networks a boolean indicating whether the between-block networks should be returned as a `combined_networks` object or not (default is TRUE)
#' @return a list of networks
#' @examples
#' # Load an embedded network object.
#' data(toyNet)
#' get_within_networks(toyNet, toyNet %v% "block")
#' @export
get_within_networks <- function(network, block, combined_networks = TRUE){
  block <- as.integer(block)
  is_directed <- network$gal$directed
  old_vertexnames <- network%v% "vertex.names"
  network%v% "vertex.names" <- 1:length(network%v% "vertex.names")
  # Get the edge list and constrain them to be edges where the involved actors are in the same block
  edge_list <- network::as.edgelist(network)
  edge_list <- cbind(edge_list,block[edge_list[,1]],block[edge_list[,2]])
  # Constraint the edge list to be within block
  edge_list <- edge_list[edge_list[,3] == edge_list[,4],]
  # Get all attributes in a data.format
  attributes <- list.vertex.attributes(network)
  attributes <- attributes[attributes != "vertex.names"]
  attributes <- attributes[attributes != "block"]
  data_attribute <- data.frame(vertex.names = network %v% "vertex.names", "block" = block)
  for(i in 1:length(attributes)){
    data_attribute <- cbind(data_attribute, get.vertex.attribute(network, attributes[i]))
    names(data_attribute)[i+2] <- attributes[i]
  }
  if(!combined_networks){
    edge_list <- as.data.frame(edge_list[,c(1,2)])
    tmp <- as.network(edge_list, vertices = data_attribute, directed = is_directed)
    tmp%v% "vertex.names" <- old_vertexnames[tmp%v% "vertex.names"]
    return(tmp)
  }  else{
    networks <- list()
    unique_block <- sort(unique(block))
    for(i in 1:length(unique_block)){
      # cat(i,"\n")
      tmp_edges <- as.data.frame(matrix(edge_list[edge_list[,3] == unique_block[i],1:2],ncol = 2))
      
      tmp_attributes <- data_attribute[data_attribute$block == unique_block[i],]
      # tmp_edges[,1] <- match(tmp_edges[,1],tmp_attributes$vertex.names)
      # tmp_edges[,2] <- match(tmp_edges[,2],tmp_attributes$vertex.names)
      if(nrow(tmp_edges)>0){
        tmp <- as.network(tmp_edges, vertices = tmp_attributes, directed = is_directed)
        tmp%v% "vertex.names" <- old_vertexnames[tmp%v% "vertex.names"]
        networks[[i]] <- tmp
      } 
      else {
        tmp <- as.network(data.frame(tmp_attributes$vertex.names[1],
                                               tmp_attributes$vertex.names[2]),
                                    vertices = tmp_attributes, directed = is_directed)
        delete.edges(tmp,eid = 1)
        tmp%v% "vertex.names" <- old_vertexnames[tmp%v% "vertex.names"]
        networks[[i]] <- tmp
      }
    }
    return(Networks(networks))
  }
}


#' Obtain the between-block networks defined by the block attribute.
#' 
#' @description
#' Function to return a list of networks, each network representing the within-block network of a block.
#' 
#' @param network a network object
#' @param block a vector of integers representing the block of each node
#' @return a list of networks
#' @examples
#' # Load an embedded network object.
#' data(toyNet)
#' get_within_networks(toyNet, toyNet %v% "block")
#' @export
get_between_networks <- function(network, block){
  block <- as.integer(block)
  is_directed <- network$gal$directed
  
  # Get the edge list and constrain them to be edges where the involved actors are in the same block
  edge_list <- network::as.edgelist(network)
  edge_list <- cbind(edge_list,block[edge_list[,1]],block[edge_list[,2]])
  # Constraint the edge list to be within block
  edge_list <- edge_list[edge_list[,3] != edge_list[,4],]
  # Get all attributes in a data.format
  attributes <- list.vertex.attributes(network)
  attributes <- attributes[attributes != "vertex.names"]
  attributes <- attributes[attributes != "block"]
  data_attribute <- data.frame(vertex.names = network %v% "vertex.names", "block" = block)
  for(i in 1:length(attributes)){
    data_attribute <- cbind(data_attribute, get.vertex.attribute(network, attributes[i]))
    names(data_attribute)[i+2] <- attributes[i]
  }
  return(as.network(edge_list, vertices = data_attribute, directed = is_directed))
}



utils::globalVariables(".") 
