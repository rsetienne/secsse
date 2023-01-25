context("test_cla_timezones")

test_that("secsse gives the same result as GeoSSE", {
  pars <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
  names(pars) <- c("sA", "sB", "sAB", "xA", "xB", "dA", "dB")
  #set.seed(5)
  #phy <- diversitree::tree.geosse(pars, max.t=4, x0=0)
  phy <- NULL
  rm(phy)
  utils::data('example_phy_GeoSSE', package = 'secsse')
  traits <- as.numeric(phy$tip.state)
  testit::assert(!is.null(phy))
  lik.g <- diversitree::make.geosse(phy, phy$tip.state)
  pars.g <- c(1.5, 0.5, 1.0, 0.7, 0.7, 1.4, 1.3)
  names(pars.g) <- diversitree::argnames(lik.g)
  lik.c <- diversitree::make.classe(phy, phy$tip.state + 1, 3)
  pars.c <- 0 * diversitree::starting.point.classe(phy, 3)
  pars.c['lambda222'] <- pars.c['lambda112'] <- pars.g['sA']
  pars.c['lambda333'] <- pars.c['lambda113'] <- pars.g['sB']
  pars.c['lambda123'] <- pars.g['sAB']
  pars.c['mu2'] <- pars.c['q13'] <- pars.g['xA']
  pars.c['mu3'] <- pars.c['q12'] <- pars.g['xB']
  pars.c['q21'] <- pars.g['dA']
  pars.c['q31'] <- pars.g['dB']
  lik.g(pars.g) # -175.7685
  classe_diversitree_LL <- lik.c(pars.c) # -175.7685
  
  ## Secsse part 
  lambdas <- list()
  lambdas[[1]] <- matrix(0,ncol = 3,nrow = 3,byrow = TRUE)
  #lambdas[[1]][1,1] <- 1.5
  lambdas[[1]][2,1] <- 1.5
  lambdas[[1]][3,1] <- 0.5
  lambdas[[1]][3,2] <- 1
  lambdas[[2]] <- matrix(0,ncol = 3,nrow = 3,byrow = TRUE)
  lambdas[[2]][2,2] <- 1.5
  #lambdas[[2]][2,2] <- 1.1
  lambdas[[3]] <- matrix(0,
                         ncol = 3,
                         nrow = 3,
                         byrow = TRUE)
  lambdas[[3]][3,3] <- 0.5
  #lambdas[[3]][3,3] <- 1
  
  mus <- c(0, 0.7, 0.7)
  
  q <- matrix(0, ncol = 3, nrow = 3, byrow = TRUE)
  q[2, 1] <- 1.4
  q[3, 1] <- 1.3
  q[1, 2] <- 0.7
  q[1, 3] <- 0.7
  
  parameter <- list()
  parameter[[1]] <- lambdas
  parameter[[2]] <- mus
  parameter[[3]] <- q
  
  num_concealed_states <- 3.1
  
  secsse_cla_LL <- secsse::cla_secsse_loglik(parameter,
                                             phy,
                                             traits,
                                             num_concealed_states,
                                             cond = "maddison_cond",
                                             root_state_weight = "maddison_weights",
                                             sampling_fraction = c(1, 1, 1),
                                             setting_calculation = NULL,
                                             see_ancestral_states = FALSE,
                                             loglik_penalty = 0)
  
  # now, we add tests for compound / arbitrary length:
  
  parameter <- list()
  parameter[[1]] <- lambdas
  parameter[[2]] <- mus
  parameter[[3]] <- q
  
  param_compound <- list()
  param_compound[[1]] <- parameter
  param_compound[[2]] <- parameter
  
  secsse_cla_timezones_LL <- 
    cla_secsse_loglik_timezones(param_compound,
                                        phy,
                                        traits,
                                        critical_t =  c(0.4),
                                        num_concealed_states,
                                        cond = "maddison_cond",
                                        root_state_weight = "maddison_weights",
                                        sampling_fraction = c(1, 1, 1),
                                        setting_calculation = NULL,
                                        see_ancestral_states = FALSE,
                                        loglik_penalty = 0)
  
  testthat::expect_equal(secsse_cla_LL, secsse_cla_timezones_LL)
  
  # we change the second parameter set, likelihood should change
  param_compound[[2]][[2]] <- c(0.0, 0.0, 0.0)
  
  secsse_cla_timezones_LL <- 
    secsse::cla_secsse_loglik_timezones(param_compound,
                                        phy,
                                        traits,
                                        critical_t =  c(1.5),
                                        num_concealed_states,
                                        cond = "maddison_cond",
                                        root_state_weight = "maddison_weights",
                                        sampling_fraction = c(1, 1, 1),
                                        setting_calculation = NULL,
                                        see_ancestral_states = FALSE,
                                        loglik_penalty = 0)
  
  testthat::expect_false(secsse_cla_LL == secsse_cla_timezones_LL)
})
