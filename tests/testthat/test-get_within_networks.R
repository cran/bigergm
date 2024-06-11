test_that("getting within networks works", {
  # Which node belongs to which block
  df <-
    tibble::tibble(
      node_id = c("A", "B", "C", "D", "E", "F", "G", "H"), 
      block = c(1, 1, 1, 2, 2, 2, 2, 1)
    )

  # Edgelist
  edgelist <-
    data.frame(
      source_id = c("A", "A", "A", "B"),
      target_id = c("B", "C", "E", "F")
    )

  network <- as.network(edgelist, vertices = df, directed = FALSE)
  # network %v% "vertex.names" <- 1:length(network %v% "vertex.names")
  # When not all nodes are isolated in a block
  subgraphs <- get_within_networks(network =network,block = network%v% "block")
  adj_true <- matrix(c(0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), nrow = 4, ncol = 4)
  expect_equal(adj_true, network::as.matrix.network.adjacency(subgraphs)[1:4,1:4], check.attributes = FALSE)
  adj_true <- matrix(0, nrow = 4, ncol = 4)
  expect_equal(adj_true, network::as.matrix.network.adjacency(subgraphs)[1:4+4,1:4+4], check.attributes = FALSE)
  # check if the vertex.names are correct (the order must be changed since D is in block 2 and H in block 1) 
  expect_equal(subgraphs %v% "vertex.names", c("A", "B", "C", "H", "D", "E", "F", "G"))
  # Next test the function where no combined network is returned just the same network without the between-block connections
  within_networks <- get_within_networks(network =network,block = network%v% "block",combined_networks = FALSE)
  # The vertex.names should be unchanged
  expect_equal(within_networks %v% "vertex.names",df$node_id)
  # The edges as well (there are only two within-block edges -> from A-B and A-C, i.e., vertices 1-2 and 1-3)
  expect_equal(as.edgelist(within_networks)[1,],c(1,2))
  expect_equal(as.edgelist(within_networks)[2,],c(1,3))
  
})
