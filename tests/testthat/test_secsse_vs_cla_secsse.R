test_that("secsse gives the same result as cla_secsse", {
  Sys.unsetenv("R_TESTS")
  
  utils::data("example_phy_GeoSSE", package = "secsse")
  traits <- as.numeric(example_phy_GeoSSE$tip.state)
  
  lambdas <- list()
  lambdas[[1]] <- matrix(0, ncol = 3, nrow = 3, byrow = TRUE)
  lambdas[[2]] <- lambdas[[1]]
  lambdas[[3]] <- lambdas[[1]]
  lambdas[[1]][1, 1] <- 1.5
  lambdas[[2]][2, 2] <- 0.5
  lambdas[[3]][3, 3] <- 1
  
  mus <- c(0.7, 0.7, 0.7)
  
  q <- matrix(0, ncol = 3, nrow = 3, byrow = TRUE)
  q[2, 1] <- 1.4
  q[3, 1] <- 1.3
  q[1, 2] <- 0.7
  q[1, 3] <- 0.7
  
  parameter <- list()
  parameter[[1]] <- lambdas
  parameter[[2]] <- mus
  parameter[[3]] <- q
  
  num_concealed_states <- 3
  
  num_modeled_traits <- ncol(q) / floor(num_concealed_states)
  
  setting_calculation <- build_initStates_time(phy = example_phy_GeoSSE,
                                               traits = traits,
                                               num_concealed_states = num_concealed_states,
                                               sampling_fraction = c(1, 1, 1),
                                               is_complete_tree = FALSE,
                                               mus = mus,
                                               num_unique_traits = num_modeled_traits,
                                               first_time = TRUE)
  states <- setting_calculation$states
  d <- ncol(states) / 2
  new_states <- states[, c(1, 2, 3, 10, 11, 12)]
  states <- new_states
  setting_calculation$states <- states

  cla_secsse_LL <- cla_secsse_loglik(parameter = parameter,
                                     phy = example_phy_GeoSSE,
                                     traits = traits,
                                     num_concealed_states = 3,
                                     cond = "proper_cond",
                                     root_state_weight = "proper_weights",
                                     sampling_fraction = c(1, 1, 1),
                                     setting_calculation = setting_calculation,
                                     see_ancestral_states = FALSE,
                                     loglik_penalty = 0,
                                     is_complete_tree = FALSE,
                                     num_threads = 1,
                                     method = "odeint::bulirsch_stoer",
                                     atol = 1e-8,
                                     rtol = 1e-7)
  
  pars <- parameter
  pars[[1]] <- c(lambdas[[1]][1, 1], lambdas[[2]][2, 2], lambdas[[3]][3, 3])
  
  secsse_LL <- secsse_loglik(parameter = pars,
                             phy = example_phy_GeoSSE,
                             traits = traits,
                             num_concealed_states = num_concealed_states,
                             cond = "proper_cond",
                             root_state_weight = "proper_weights",
                             sampling_fraction = c(1,1,1),
                             setting_calculation = setting_calculation,
                             see_ancestral_states = FALSE,
                             loglik_penalty = 0,
                             is_complete_tree = FALSE,
                             num_threads = 1,
                             atol = 1e-8,
                             rtol = 1e-7,
                             method = "odeint::bulirsch_stoer")
  
  testthat::expect_equal(secsse_LL, cla_secsse_LL,tolerance = 1e-5)
})
