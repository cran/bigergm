test_that("computing the lower bound for directed networks works", {
  set.seed(1234)
  #################### 
  # Setup
  ####################
  # Number of nodes
  N <- 12
  # Number of clusters
  K <- 3

  # Create an adjacency matrix
  edgelist <-
    tibble::tibble(
      tail = 1:N,
      head = 1:N
    ) %>%
    tidyr::expand(tail, head) %>%
    dplyr::filter(tail != head) %>%
    dplyr::mutate(connect = as.integer(unlist(rbinom(n = nrow(.), prob = 0.5, size = 1)))) %>%
    dplyr::filter(connect == 1)

  net <- network::network(edgelist, matrix.type = "edgelist", directed = TRUE)
  adj <- network::as.matrix.network.adjacency(net)
  adj <- as(adj, "dgCMatrix")

  # Create feature matrices
  x <- as.integer(unlist(rbinom(n = N,size = 1,prob = 0.5)))
  S <- Matrix::sparseMatrix(i = {}, j = {}, dims = c(N, N))
  S <- as(S, "dMatrix")
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        s_ij <- ifelse(x[i] == x[j], 1, 0)
        S[i, j] <- s_ij
      }
    }
  }

  y <- as.integer(unlist(rbinom(size = 1,prob = 0.5,n = N)))
  V <- Matrix::sparseMatrix(i = {}, j = {}, dims = c(N, N))
  V <- as(V, "dgCMatrix")
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        v_ij <- ifelse(y[i] == y[j], 1, 0)
        V[i, j] <- v_ij
      }
    }
  }

  z <- as.integer(unlist(rbinom(size = 1,prob = 0.5,n = N)))
  W <- Matrix::sparseMatrix(i = {}, j = {}, dims = c(N, N))
  W <- as(W, "dgCMatrix")
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        w_ij <- ifelse(z[i] == z[j], 1, 0)
        W[i, j] <- w_ij
      }
    }
  }

  # Create a N x K matrix whose (i, k) element represents the probability that node i belongs to block k.
  tau <-
    matrix(c(
      0.2, 0.5, 0.3,
      0.4, 0.4, 0.2,
      0.1, 0.4, 0.5,
      0.4, 0.4, 0.2,
      0.1, 0.1, 0.8,
      0.05, 0.05, 0.9,
      0.8, 0.1, 0.1,
      0.3, 0.4, 0.3,
      0.1, 0.8, 0.1,
      0.5, 0.4, 0.1,
      0.3, 0.3, 0.4,
      0.8, 0.1, 0.1
    ),
    nrow = K, ncol = N
    )
  tau <- t(tau)

  # Compute gamma (parameter of multinomial distribution)
  alpha <- colMeans(tau)

  ###########################################################
  # Compute the lower bound in a naive way
  ###########################################################
  # Compute pi for D_ij = 1
  minPi <- 1e-4
  list_pi <- list()
  for (w in 0:1) {
    for (v in 0:1) {
      for (s in 0:1) {
        print(glue::glue("Compute pi for pi_s{s}v{v}w{w}"))
        denom <- matrix(0, nrow = K, ncol = K)
        num <- matrix(0, nrow = K, ncol = K)
        index <- s + 2 * v + 4 * w + 1
        for (k in 1:K) {
          for (l in 1:K) {
            for (i in 1:N) {
              for (j in 1:N) {
                if (i != j & S[i, j] == s & V[i, j] == v & W[i, j] == w) {
                  denom[k, l] <- denom[k, l] + tau[i, k] * tau[j, l]
                }
                if (i != j & adj[i, j] == 1 & S[i, j] == s & V[i, j] == v & W[i, j] == w) {
                  num[k, l] <- num[k, l] + tau[i, k] * tau[j, l]
                }
              }
            }
          }
        }
        pi <- num / denom
        # Remove extremely small elements in pi
        for (k in 1:K) {
          for (l in 1:K) {
            if (pi[k, l] < minPi) {
              pi[k, l] <- minPi
            }
          }
        }
        list_pi[[index]] <- pi
      }
    }
  }

  # Compute the true lower bound
  LB_true <- 0
  # First term
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        # For each ij, determine which pi must be used.
        index_ij <- S[i, j] + 2 * V[i, j] + 4 * W[i, j] + 1
        pi_ij <- list_pi[[index_ij]]
        # if D_ij = 0, replace pi with 1 - pi.
        if (adj[i, j] == 0) {
          pi_ij <- 1 - pi_ij
        }
        for (k in 1:K) {
          for (l in 1:K) {
            LB_true <- LB_true + tau[i, k] * tau[j, l] * log(pi_ij[k, l])
          }
        }
      }
    }
  }

  # Second term
  for (i in 1:N) {
    for (k in 1:K) {
      LB_true <- LB_true + tau[i, k] * (log(alpha[k]) - log(tau[i, k]))
    }
  }

  ###########################################################
  # Compute the lower bound using the cpp function
  ###########################################################
  list_feature_adjmat <- list(S, V, W)
  list_multiplied_feature_adjmat <- get_elementwise_multiplied_matrices_R(adj, list_feature_adjmat)
  denom <- get_matrix_for_denominator_R(N, list_feature_adjmat)
  list_multiplied_feature_adjmat[[1]] <- denom
  list_multiplied_feature_adjmat <- lapply(list_multiplied_feature_adjmat, function(x) {x*1})
  
  alpha_LB <- run_MM_with_features(N, K, alpha, 
                                   list_multiplied_feature_adjmat, tau, verbose = 3, directed = TRUE)

  
  # Check if it works
  expect_equal(alpha_LB[[2]], LB_true, tolerance = 1e-10)


  
  })


test_that("computing the lower bound for directed networks without features works", {
  ####################
  # Setup
  ####################
  # Number of nodes
  N <- 12
  # Number of clusters
  K <- 3
  set.seed(123)
  # Create an adjacency matrix
  edgelist <-
    tibble::tibble(
      tail = 1:N,
      head = 1:N
    ) %>%
    tidyr::expand(tail, head) %>%
    dplyr::filter(tail != head) %>%
    dplyr::mutate(connect = as.integer(unlist(rbinom(size = 1,prob = 0.5,n = nrow(.))))) %>%
    dplyr::filter(connect == 1)

  net <- network::network(edgelist, matrix.type = "edgelist", directed = TRUE)
  adj <- network::as.matrix.network.adjacency(net)
  adj <- as(adj, "dgCMatrix")

  # Create a N x K matrix whose (i, k) element represents the probability that node i belongs to block k.
  tau <-
    matrix(c(
      0.2, 0.5, 0.3,
      0.4, 0.4, 0.2,
      0.1, 0.4, 0.5,
      0.4, 0.4, 0.2,
      0.1, 0.1, 0.8,
      0.05, 0.05, 0.9,
      0.8, 0.1, 0.1,
      0.3, 0.4, 0.3,
      0.1, 0.8, 0.1,
      0.5, 0.4, 0.1,
      0.3, 0.3, 0.4,
      0.8, 0.1, 0.1
    ),
    nrow = K, ncol = N
    )
  tau <- t(tau)

  # Compute gamma (parameter of multinomial distribution)
  alpha <- colMeans(tau)

  ###########################################################
  # Compute the lower bound in a naive way
  ###########################################################
  # Compute pi for D_ij = 1
  minPi <- 1e-4
  denom <- matrix(0, nrow = K, ncol = K)
  num <- matrix(0, nrow = K, ncol = K)
  for (k in 1:K) {
    for (l in 1:K) {
      for (i in 1:N) {
        for (j in 1:N) {
          if (i != j) {
            denom[k, l] <- denom[k, l] + tau[i, k] * tau[j, l]
            if (i != j & adj[i, j] == 1) {
              num[k, l] <- num[k, l] + tau[i, k] * tau[j, l]
            }
          }
        }
      }
    }
  }
  pi <- num / denom

  # Remove extremely small elements in pi
  for (k in 1:K) {
    for (l in 1:K) {
      if (pi[k, l] < minPi) {
        pi[k, l] <- minPi
      }
    }
  }


  # Compute the true lower bound
  LB_true_part_1 <- 0
  # First term
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        # if D_ij = 0, replace pi with 1 - pi.
        if (adj[i, j] == 0) {
          pi_ij <- 1 - pi
        } else {
          pi_ij <- pi
        }
        for (k in 1:K) {
          for (l in 1:K) {
            LB_true_part_1 <- LB_true_part_1 + tau[i, k] * tau[j, l] * log(pi_ij[k, l])
          }
        }
      }
    }
  }

  LB_true_part_2 <- 0
  # Second term
  for (i in 1:N) {
    for (k in 1:K) {
      LB_true_part_2 <- LB_true_part_2 + tau[i, k] * (log(alpha[k]) - log(tau[i, k]))
    }
  }

  ###########################################################
  # Compute the lower bound using the cpp function
  ###########################################################
  alpha_LB <- run_MM_without_features(numOfVertices = N, numOfClasses = K, 
                                      alpha = alpha, tau = tau, network = adj, verbose = 0, directed = T)


  expect_equal(alpha_LB[[2]],   LB_true_part_1 +LB_true_part_2, tolerance = 1e-10)
  
  
  })
