test_that("test secsse_sim", {
  
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
  max_spec <- 10000
  num_repl <- 100
  
  max_time <- 1
  
  tree1 <- secsse::secsse_sim(lambdas = lambdas,
                              mus = mus,
                              qs = qs,
                              num_concealed_states = num_concealed_states,
                              crown_age = max_time,
                              max_spec = max_spec,
                              conditioning = "obs_states",
                              seed = 42)
  
  all_obs_present <- c(0, 1, 2) %in% tree1$obs_traits
  testthat::expect_equal(sum(all_obs_present), 3)
  
  tree2 <- secsse::secsse_sim(lambdas = lambdas,
                              mus = mus,
                              qs = qs,
                              num_concealed_states = num_concealed_states,
                              crown_age = max_time,
                              max_spec = max_spec,
                              conditioning = "true_states",
                              seed = 43)
  
  all_obs_present <- names(mus) %in% tree2$true_traits
  testthat::expect_equal(sum(all_obs_present), 9)
  
  if (requireNamespace("ape")) {
    testthat::expect_equal(max(ape::branching.times(tree1$phy)), 1)
  }
  
  # custom conditioning
  tree3 <- secsse::secsse_sim(lambdas = lambdas,
                              mus = mus,
                              qs = qs,
                              num_concealed_states = num_concealed_states,
                              crown_age = max_time,
                              max_spec = max_spec,
                              conditioning = c(0, 1),
                              seed = 444)
  traits_present <- c(0, 1) %in% tree3$obs_traits
  testthat::expect_equal(sum(traits_present), 2)
})

test_that("test secsse_sim 2", {
  lambda_shift <- secsse::create_default_lambda_transition_matrix()
  lambda_list <- secsse::create_lambda_list(transition_matrix = lambda_shift)
  mus <- secsse::create_mu_vector(state_names = c(0, 1),
                                  num_concealed_states = 2,
                                  lambda_list = lambda_list)
  q_mat <- secsse::create_default_shift_matrix(mu_vector = mus)
  q_mat <- secsse::create_q_matrix(state_names = c(0, 1),
                                   num_concealed_states = 2,
                                   shift_matrix = q_mat)  
  
  pars <- c(0.5, 0.3, 0.7, 0.1, 0.1)
  lambda_p <- secsse::fill_in(lambda_list, pars)
  mu_p <- secsse::fill_in(mus, pars)
  q_mat_p <- secsse::fill_in(q_mat, pars)
  
  focal_tree <- secsse::secsse_sim(lambdas = lambda_p,
                                   mus = mu_p,
                                   qs = q_mat_p,
                                   crown_age = 10,
                                   num_concealed_states = 2,
                                   max_spec = 100,
                                   seed = 21,
                                   drop_extinct = FALSE)
  if (requireNamespace("geiger")) {
    vx <- geiger::is.extinct(focal_tree$phy)
    testthat::expect_true(length(vx) > 0)
  }
})


test_that("test secsse_sim pool_init_states and complete tree", {
  lambda_shift <- secsse::create_default_lambda_transition_matrix()
  lambda_list <- secsse::create_lambda_list(transition_matrix = lambda_shift)
  mus <- secsse::create_mu_vector(state_names = c(0, 1),
                                  num_concealed_states = 2,
                                  lambda_list = lambda_list)
  q_mat <- secsse::create_default_shift_matrix(mu_vector = mus)
  q_mat <- secsse::create_q_matrix(state_names = c(0, 1),
                                   num_concealed_states = 2,
                                   shift_matrix = q_mat)  
  
  pars <- c(0.5, 0.3, 0.7, 0.1, 0.1)
  lambda_p <- secsse::fill_in(lambda_list, pars)
  mu_p <- secsse::fill_in(mus, pars)
  q_mat_p <- secsse::fill_in(q_mat, pars)
  
  focal_tree <- secsse::secsse_sim(lambdas = lambda_p,
                                   mus = mu_p,
                                   qs = q_mat_p,
                                   crown_age = 10,
                                   num_concealed_states = 2,
                                   pool_init_states = c("0A"),
                                   max_spec = 100,
                                   seed = 21,
                                   drop_extinct = FALSE)
  testthat::expect_true(focal_tree$initialState == "0A")
  
  focal_tree <- secsse::secsse_sim(lambdas = lambda_p,
                                   mus = mu_p,
                                   qs = q_mat_p,
                                   crown_age = 10,
                                   num_concealed_states = 2,
                                   pool_init_states = c("0A", "1B"),
                                   max_spec = 100,
                                   seed = 21,
                                   drop_extinct = FALSE)
  testthat::expect_true(focal_tree$initialState %in% c("0A", "1B"))
  
  
  pars <- c(0.5, 0.3, 0.2, 0.1, 0.1)
  # lambda_p <- secsse::fill_in(lambda_list, pars)
  mu_p <- secsse::fill_in(mus, pars)
  q_mat_p <- secsse::fill_in(q_mat, pars)
  focal_tree <- secsse::secsse_sim(lambdas = lambda_p,
                                   mus = mu_p,
                                   qs = q_mat_p,
                                   crown_age = 10,
                                   num_concealed_states = 2,
                                   max_spec = 100,
                                   min_spec = 99,
                                   max_species_extant = FALSE,
                                   seed = 21,
                                   drop_extinct = FALSE)
  testthat::expect_equal(length(focal_tree$phy$tip.label), 100)
  if (requireNamespace("geiger")) {
    vx <- geiger::is.extinct(focal_tree$phy)
    testthat::expect_true(length(vx) > 0)
  }
})
