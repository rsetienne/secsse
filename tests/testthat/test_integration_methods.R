test_that("loglik for different integrators", {
  set.seed(42)
  out <- DDD::dd_sim(pars = c(0.4, 0.1, 40), age = 15)
  phy <- out$tes
  traits <- sample(c(0, 1), ape::Ntip(phy), replace = TRUE)
  b <- c(0.04, 0.04)  # lambda
  d <- rep(0, 2)
  userTransRate <- 0.2 # transition rate among trait states
  num_concealed_states <- 2
  sampling_fraction <- c(1, 1)
  toCheck <- secsse::id_paramPos(traits, num_concealed_states)
  toCheck[[1]][] <- b
  toCheck[[2]][] <- d
  toCheck[[3]][, ] <- userTransRate
  diag(toCheck[[3]]) <- NA
  root_state_weight <- "maddison_weights"
  cond <- "noCondit"
  
  loglik1 <- testthat::expect_warning(as.numeric(secsse::secsse_loglik(parameter = toCheck,
                                      phy = phy,
                                      traits = traits,
                                      num_concealed_states =
                                        num_concealed_states,
                                      cond = cond,
                                      root_state_weight = root_state_weight,
                                      sampling_fraction = sampling_fraction)
  ))
  
  for (integ_method in c("odeint::runge_kutta_cash_karp54", 
                         "odeint::runge_kutta_fehlberg78", 
                         "odeint::runge_kutta_dopri5", 
                         "odeint::bulirsch_stoer",
                         "odeint::runge_kutta4")) {
    loglik2 <- testthat::expect_warning(as.numeric(secsse::secsse_loglik(parameter = toCheck,
                                        phy = phy,
                                        traits = traits,
                                        num_concealed_states =
                                          num_concealed_states,
                                        cond = cond,
                                        root_state_weight = root_state_weight,
                                        sampling_fraction = sampling_fraction,
                                        method = integ_method)))
    testthat::expect_equal(loglik1, loglik2, tolerance = 0.01)
  }
})  