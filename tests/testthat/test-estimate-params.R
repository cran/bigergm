set.seed(334)
# Prepare data
edgelist <-
  tibble::tribble(
    ~head, ~tail,
    1, 9,
    2, 6,
    2, 7,
    2, 9,
    3, 5,
    3, 9,
    4, 7,
    4, 11,
    4, 15,
    5, 11,
    5, 15,
    7, 8,
    7, 16,
    9, 13,
    9, 14,
    9, 16,
    10, 14,
    11, 15,
    13, 15,
    13, 16
  )
edgelist <-
  as.matrix(edgelist)
attr(edgelist, "n") <- 16
attr(edgelist, "vnames") <-
  c(
    "Acciaiuoli", "Albizzi", "Barbadori", "Bischeri", "Castellani", "Ginori",
    "Guadagni", "Lamberteschi", "Medici", "Pazzi", "Peruzzi", "Pucci", "Ridolfi",
    "Salviati", "Strozzi", "Tornabuoni"
  )
attr(edgelist, "directed") <- FALSE
attr(edgelist, "bipartite") <- FALSE
attr(edgelist, "loops") <- FALSE
attr(edgelist, "class") <- c("edgelist", "matrix")
g <- network::network(edgelist, matrix.type = "edgelist", directed = FALSE)
g%v%"vertex.names" <- 1:length(g%v%"vertex.names")
x1 <- as.integer(unlist(rbinom(size = 1,prob = 0.5,n = g$gal$n)))
network::set.vertex.attribute(x = g, attrname = "x1", value = x1)

# Cluster
z_memb <- rep(1:4, each = 4)
network::set.vertex.attribute(x = g, attrname = "block", value = z_memb)

# Create dataset for test
g_link <- intergraph::asDF(g)$edges
g_attr <- intergraph::asDF(g)$vertexes



df_g <-
  tibble::tibble(
    head = 1:g$gal$n,
    tail = 1:g$gal$n
  ) %>%
  tidyr::expand(.data$tail, .data$head) %>%
  dplyr::filter(.data$tail < .data$head) %>%
  dplyr::left_join(., g_attr, by = c("tail" = "intergraph_id")) %>%
  dplyr::left_join(., g_attr, by = c("head" = "intergraph_id")) %>%
  dplyr::mutate(
    nodematch.x1 = ifelse(x1.x == x1.y, 1, 0),
    same_block = ifelse(block.x == block.y, 1, 0)
  ) %>%
  dplyr::select("tail", "head", "nodematch.x1":"same_block") %>%
  dplyr::left_join(., g_link, by = c("tail" = "V1", "head" = "V2")) %>%
  dplyr::mutate(connected = ifelse(is.na(na), 0, 1)) %>%
  dplyr::select(-"na")



# Estimate the model
formula <- g ~ edges + nodematch("x1") + triangle + kstar(2)
g %v% "block" <- z_memb

est_between <- est_between(
  formula = formula,
  network = g,
  add_intercepts = FALSE
)

est_within <- est_within(
  formula = g ~ edges + nodematch("x1"),
  network = g, seed = 1,
  method_within = "MPLE", 
  add_intercepts = FALSE
)


test_that("estimating between-block parameters by logit works", {
  # Check if between-block connections are all zero.
  g_logit <- g
  edgelist <- intergraph::asDF(g_logit)$edges

  true_edgelist <-
    df_g %>%
    dplyr::filter(same_block == 0 & connected == 1) %>%
    dplyr::select("tail", "head") %>%
    dplyr::arrange(tail, head)
  
  # Check if estimates for between-block parameters are the same.
  param_est <- stats::coef(est_between)
  
  logit_true_between <- glm(
    formula = connected ~ nodematch.x1,
    data = df_g[df_g$same_block ==0,],
    family = "binomial"
  )
  logit_true_within <- glm(
    formula = connected ~ nodematch.x1,
    data = df_g[df_g$same_block ==1,],
    family = "binomial"
  )
  
  # Does it work?
  expect_equal(est_within$coefficients, 
               logit_true_within$coefficients,
               check.attributes = FALSE, tolerance = 1e-7)
  expect_equal(est_between$coefficients, 
               logit_true_between$coefficients, 
               check.attributes = FALSE, tolerance = 1e-7)
  
})
