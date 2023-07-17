context("test_secsse_sim")

test_that("test secsse_sim", {
  testthat::skip_on_cran()
  
  parenthesis <- "(((6:0.2547423371,(1:0.0496153503,4:0.0496153503):0.2051269868):0.1306304758,(9:0.2124135406,5:0.2124135406):0.1729592723):1.151205247,(((7:0.009347664296,3:0.009347664296):0.2101416075,10:0.2194892718):0.1035186448,(2:0.2575886319,8:0.2575886319):0.06541928469):1.213570144);" # nolint
  phylotree <- ape::read.tree(file = "", parenthesis)
  traits <- c(2, 0, 1, 0, 2, 0, 1, 2, 2, 0)
  num_concealed_states <- 3
  idparslist <- cla_id_paramPos(traits, num_concealed_states)
  idparslist$lambdas[2, ] <- rep(1, 9)
  idparslist[[2]][] <- 4
  masterBlock <- matrix(5, ncol = 3, nrow = 3, byrow = TRUE)
  diag(masterBlock) <- NA
  diff.conceal <- FALSE
  idparslist[[3]] <- q_doubletrans(traits, masterBlock, diff.conceal)
  testthat::expect_output(
    startingpoint <- DDD::bd_ML(brts = ape::branching.times(phylotree))
  )
  intGuessLamba <- startingpoint$lambda0
  intGuessMu <- startingpoint$mu0
  idparsopt <- c(1)
  initparsopt <- c(rep(intGuessLamba, 1))
  idparsfix <- c(0, 4, 5)
  parsfix <- c(0, 0, 0.01)
  tol <- c(1e-04, 1e-05, 1e-07)
  maxiter <- 1000 * round((1.25)^length(idparsopt))
  optimmethod <- "subplex"
  cond <- "proper_cond"
  root_state_weight <- "proper_weights"
  sampling_fraction <- c(1, 1, 1)
  
  testthat::expect_warning(
    model_R <- secsse::cla_secsse_ml(
      phylotree,
      traits,
      num_concealed_states,
      idparslist,
      idparsopt,
      initparsopt,
      idparsfix,
      parsfix,
      cond,
      root_state_weight,
      sampling_fraction,
      tol,
      maxiter,
      optimmethod,
      num_cycles = 1,
      verbose = FALSE)
  )
  
  qs <- model_R$MLpars[[3]]
  diag(qs) <- 0
  
  lambdas <- model_R$MLpars[[1]]
  mus <- model_R$MLpars[[2]]
  maxSpec <- 10000
  num_repl <- 100
  
  max_time <- 1
  
  tree1 <- secsse::secsse_sim(lambdas = lambdas,
                              mus = mus,
                              qs = qs,
                              num_concealed_states = num_concealed_states,
                              crown_age = max_time,
                              maxSpec = maxSpec,
                              conditioning = "obs_states")
  
  all_obs_present <- c(0, 1, 2) %in% tree1$obs_traits
  testthat::expect_equal(sum(all_obs_present), 3)
  
  tree2 <- secsse::secsse_sim(lambdas = lambdas,
                              mus = mus,
                              qs = qs,
                              num_concealed_states = num_concealed_states,
                              crown_age = max_time,
                              maxSpec = maxSpec,
                              conditioning = "true_states")
  
  all_obs_present <- names(mus) %in% tree2$true_traits
  testthat::expect_equal(sum(all_obs_present), 9)
  
  if (requireNamespace("ape")) {
    testthat::expect_equal(max(ape::branching.times(tree1$phy)), 1)
  }
})

test_that("test secsse_sim() with extinct species", {
  
  spec_matrix <- c(0, 0, 0, 1)
  spec_matrix <- rbind(spec_matrix, c(1, 1, 1, 1))
  lambda_list <- secsse::create_lambda_list(state_names = c(0, 1),
                                            num_concealed_states = 2,
                                            transition_matrix = spec_matrix,
                                            model = "CR")
  
  mu_vector <- secsse::create_mu_vector(state_names = c(0, 1),
                                        num_concealed_states = 2,
                                        model = "CR",
                                        lambda_list = lambda_list)
  
  shift_matrix <- c(0, 1, 3)
  shift_matrix <- rbind(shift_matrix, c(1, 0, 4))
  
  q_matrix <- secsse::create_q_matrix(state_names = c(0, 1),
                                      num_concealed_states = 2,
                                      shift_matrix = shift_matrix,
                                      diff.conceal = FALSE)
  
  
  speciation_rate <- 0.5
  extinction_rate <- 0.05
  q_ab <- 0.1
  q_ba <- 0.1
  used_params <- c(speciation_rate, extinction_rate, q_ab, q_ba)
  
  sim_lambda_list <- secsse::fill_in(lambda_list, used_params)
  sim_mu_vector   <- secsse::fill_in(mu_vector, used_params)
  sim_q_matrix    <- secsse::fill_in(q_matrix, used_params)
  
  sim_tree <- testthat::expect_silent(
    secsse::secsse_sim(lambdas = sim_lambda_list,
                       mus = sim_mu_vector,
                       qs = sim_q_matrix,
                       crown_age = 5,
                       num_concealed_states = 2,
                       seed = 5, 
                       drop_extinct = FALSE) # Keep extinct species
  )
})
