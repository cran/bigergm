## ----include = FALSE----------------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  out.width = "100%"
)

## ----eval=FALSE---------------------------------------------------------------
# install.packages("bigergm")

## ----setup, results=FALSE, message=FALSE, warning=FALSE, echo=FALSE-----------
library(bigergm)
library(ergm)
library(dplyr)


## ----message=FALSE------------------------------------------------------------
# Load the network object.
data(toyNet)
# Plot the network.
plot(toyNet, vertex.col = rep(c("tomato", "steelblue", "darkgreen", "black"),
                        each = toyNet$gal$n/4))

## -----------------------------------------------------------------------------
model_formula <- toyNet ~ edges + nodematch("x") + nodematch("y") + triangle

## ----message=FALSE------------------------------------------------------------
res <-bigergm(
    # The model you would like to estimate
    object = model_formula,
    # The number of blocks
    n_blocks =  4, 
    # The maximum number of MM algorithm steps
    n_MM_step_max = 100,
    # The tolarence for the MM algorithm
    tol_MM_step = 1e-6,
    # Perform parameter estimation after the block recovery step
    estimate_parameters = TRUE,
    # Indicate that clustering must take into account nodematch on characteristics
    clustering_with_features = TRUE,
    # Keep track of block memberships at each EM iteration
    check_block_membership = TRUE, 
    # Name the heuristic algorithm used for initializing the block memberships
    initialization = "walktrap"
)


## -----------------------------------------------------------------------------
plot(1:length(res$MM_lower_bound),
     res$MM_lower_bound, type = "l", xlab = "Iterations", ylab = "Lower Bound")

## -----------------------------------------------------------------------------
plot(res)

## -----------------------------------------------------------------------------
# For the between networks
summary(res$est_between)

## -----------------------------------------------------------------------------
# For the within networks
summary(res$est_within)

## ----message=FALSE, echo=TRUE-------------------------------------------------
simulate(res, seed = 1)

## ----message=FALSE, echo=TRUE-------------------------------------------------
sim_net <- bigergm::simulate_bigergm(
  formula = model_formula,
  # The coefficients for the between connections
  coef_between = res$est_between$coefficients,
   # The coefficients for the within connections
  coef_within = res$est_within$coefficients,
  # Number of simulations to return
  nsim = 1,
  # If `stats` a list with network statistics 
  # for the between and within connections is returned
  output = "network"
)

## -----------------------------------------------------------------------------
plot(sim_net)

## ----message=FALSE------------------------------------------------------------
gof_res <- gof(
  # The object returned by bigergm::bigergm()
  object = res,
  # The number of simulations to use
  nsim = 100, 
  # Compute the geodesic distance for the observed and each simulated network
  compute_geodesic_distance = TRUE,
  # Set a seed for reproducibility
  seed = 1234,
  # Start at the observed network
  start_from_observed = TRUE, type = "within",
  # The control parameters for the simulation
  control_within = ergm::control.simulate.formula(MCMC.burnin = 1000, MCMC.interval = 1000)
)

## ----message=FALSE, warning=FALSE---------------------------------------------
degree_gof <- 
  gof_res$simulated$degree_dist %>%
  dplyr::group_by(degree) %>%
  dplyr::summarise(log_mean_share = mean(log(share)),
                   log_sd_share = sd(log(share))) %>%
  dplyr::ungroup()
plot(degree_gof$degree, degree_gof$log_mean_share,
     xlab = "Degree", ylab = "Log Prop. of Nodes",
     ylim = c(-5.5,-1.8), xlim = c(6,17), type = "l", lty = 2)
lines(degree_gof$degree, degree_gof$log_mean_share+ 1.96 * degree_gof$log_sd_share, type = "l", lty = 2)
lines(degree_gof$degree, degree_gof$log_mean_share- 1.96 * degree_gof$log_sd_share, type = "l", lty = 2)
tmp_info <- gof_res$original$degree_dist %>% 
  dplyr::filter(share > 0 & degree < 22)
lines(tmp_info$degree, log(tmp_info$share), lty = 1)

## -----------------------------------------------------------------------------
plot(gof_res)

## ----message=FALSE------------------------------------------------------------
res_second <-
  bigergm::bigergm(object = res)

