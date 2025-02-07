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
  
  focal_tree <- secsse::secsse_sim(lambdas = lambda_p,
                                   mus = mu_p,
                                   qs = q_mat_p,
                                   crown_age = 10,
                                   num_concealed_states = 2,
                                   pool_init_states = c("0"),
                                   max_spec = 100,
                                   seed = 21,
                                   drop_extinct = FALSE)
  testthat::expect_true(focal_tree$initialState %in% c("0A", "0B"))
  
  pars <- c(0.5, 0.3, 0.2, 0.1, 0.1)
  
  mu_p <- secsse::fill_in(mus, pars)
  q_mat_p <- secsse::fill_in(q_mat, pars)

  focal_tree <- secsse::secsse_sim(lambdas = lambda_p,
                                   mus = mu_p,
                                   qs = q_mat_p,
                                   crown_age = 10,
                                   num_concealed_states = 2,
                                   max_spec = 100,
                                   min_spec = 100,
                                   max_species_extant = FALSE,
                                   seed = 21,
                                   drop_extinct = FALSE,
                                   tree_size_hist = TRUE,
                                   verbose = FALSE)
  
  testthat::expect_equal(length(focal_tree$phy$tip.label), 100)
  if (requireNamespace("geiger")) {
    vx <- geiger::is.extinct(focal_tree$phy)
    testthat::expect_true(length(vx) > 0)
  }
  testthat::expect_gt(length(focal_tree$size_hist), 0)
})


test_that("test trait shift", {
  states <- c("S", "G")
  
  spec_matrix <- c("S", "S", "S", 1)
  spec_matrix <- rbind(spec_matrix, c("G", "G", "G", 2))
  lambda_list <- secsse::create_lambda_list(state_names = states,
                                            num_concealed_states = 2,
                                            transition_matrix = spec_matrix,
                                            model = "ETD")
  
  mu_vector <- secsse::create_mu_vector(state_names = states,
                                        num_concealed_states = 2,
                                        model = "ETD",
                                        lambda_list = lambda_list)
  
  shift_matrix <- c("S", "G", 5)
  shift_matrix <- rbind(shift_matrix, c("G", "S", 6))
  
  q_matrix <- secsse::create_q_matrix(state_names = states,
                                      num_concealed_states = 2,
                                      shift_matrix = shift_matrix,
                                      diff.conceal = TRUE)
  idparslist <- list()
  idparslist[[1]] <- lambda_list
  idparslist[[2]] <- mu_vector
  idparslist[[3]] <- q_matrix
  
  spec_S <- 0
  spec_G <- 0
  ext_S <- ext_G <- 0.0
  q_SG <- 10.0
  q_GS <- 0
  used_params <- c(spec_S, spec_G, ext_S, ext_G, q_SG, q_GS, 0, 0)
  
  crown_age_used <- 5
  
  sim_lambda_list_etd <- secsse::fill_in(idparslist[[1]], used_params)
  sim_mu_vector_etd   <- secsse::fill_in(idparslist[[2]], used_params)
  sim_q_matrix_etd    <- secsse::fill_in(idparslist[[3]], used_params)
  
  
  testthat::expect_error(
  sim_tree <- secsse::secsse_sim(lambdas = sim_lambda_list_etd,
                                  mus = sim_mu_vector_etd,
                                  qs = sim_q_matrix_etd,
                                  crown_age = crown_age_used,
                                  num_concealed_states = 2,
                                  conditioning = "none",
                                  pool_init_states = c("S"))
  )

  spec_S <- 1e-10
  used_params <- c(spec_S, spec_G, ext_S, ext_G, q_SG, q_GS, 0, 0)
  sim_lambda_list_etd <- secsse::fill_in(idparslist[[1]], used_params)
  
  
  sim_tree <- secsse::secsse_sim(lambdas = sim_lambda_list_etd,
                                 mus = sim_mu_vector_etd,
                                 qs = sim_q_matrix_etd,
                                 crown_age = 100,
                                 num_concealed_states = 2,
                                 conditioning = "none",
                                 pool_init_states = c("S"))
  
  testthat::expect_true(length(which(sim_tree$obs_traits == "S")) == 0)
})

test_that("test mutate away shift", {
  states <- c("S", "G", "D")
  
  spec_matrix <- c("S", "S", "S", 1)
  spec_matrix <- rbind(spec_matrix, c("G", "G", "G", 2))
  spec_matrix <- rbind(spec_matrix, c("D", "S", "G", 3))
  lambda_list <- secsse::create_lambda_list(state_names = states,
                                            num_concealed_states = 3,
                                            transition_matrix = spec_matrix,
                                            model = "ETD")
  
  mu_vector <- secsse::create_mu_vector(state_names = states,
                                        num_concealed_states = 3,
                                        model = "ETD",
                                        lambda_list = lambda_list)
  
  shift_matrix <- c("S", "G", 7)
  shift_matrix <- rbind(shift_matrix, c("G", "S", 8))
  
  q_matrix <- secsse::create_q_matrix(state_names = states,
                                      num_concealed_states = 3,
                                      shift_matrix = shift_matrix,
                                      diff.conceal = FALSE)
  idparslist <- list()
  idparslist[[1]] <- lambda_list
  idparslist[[2]] <- mu_vector
  idparslist[[3]] <- q_matrix
  
  used_params <- rep(0, 8)
  
  used_params[3] <- 1 # spec_D
  
  crown_age_used <- 5
  
  sim_lambda_list_etd <- secsse::fill_in(idparslist[[1]], used_params)
  sim_mu_vector_etd   <- secsse::fill_in(idparslist[[2]], used_params)
  sim_q_matrix_etd    <- secsse::fill_in(idparslist[[3]], used_params)
  
  sim_tree <- secsse::secsse_sim(lambdas = sim_lambda_list_etd,
                                 mus = sim_mu_vector_etd,
                                 qs = sim_q_matrix_etd,
                                 crown_age = crown_age_used,
                                 num_concealed_states = 2,
                                 conditioning = "none",
                                 pool_init_states = c("D"))
  
  testthat::expect_true(length(sim_tree$phy$tip.label) == 2)
  
  testthat::expect_true(length(which(sim_tree$obs_traits == "D")) == 0)
  testthat::expect_true(length(which(sim_tree$obs_traits == "S")) == 1)
  testthat::expect_true(length(which(sim_tree$obs_traits == "G")) == 1)
})

test_that("test comparison musse", {
  testthat::skip_on_ci()
  testthat::skip_on_cran()
  
  if (!requireNamespace("diversitree")) testthat::skip()
  
  states <- c("S", "G", "D")
  
  spec_matrix <- c("S", "S", "S", 1)
  spec_matrix <- rbind(spec_matrix, c("G", "G", "G", 2))
  spec_matrix <- rbind(spec_matrix, c("D", "D", "D", 3))
  lambda_list <- secsse::create_lambda_list(state_names = states,
                                            num_concealed_states = 3,
                                            transition_matrix = spec_matrix,
                                            model = "ETD")
  
  mu_vector <- secsse::create_mu_vector(state_names = states,
                                        num_concealed_states = 3,
                                        model = "ETD",
                                        lambda_list = lambda_list)
  
  shift_matrix <- c("S", "G", 7)
  shift_matrix <- rbind(shift_matrix, c( "S", "D", 8))
  
  shift_matrix <- rbind(shift_matrix, c("G", "S", 9))
  shift_matrix <- rbind(shift_matrix, c("G", "D", 10))
  
  shift_matrix <- rbind(shift_matrix, c( "D", "S", 11))
  shift_matrix <- rbind(shift_matrix, c( "D", "G", 12))
  
  q_matrix <- secsse::create_q_matrix(state_names = states,
                                      num_concealed_states = 3,
                                      shift_matrix = shift_matrix,
                                      diff.conceal = TRUE)
  idparslist <- list()
  idparslist[[1]] <- lambda_list
  idparslist[[2]] <- mu_vector
  idparslist[[3]] <- q_matrix
  
  idparslist_ETD <- idparslist
  
  spec_S <- 0.1
  spec_G <- 0.15
  spec_D <- 0.2
  
  
  ext_S <- ext_G <- ext_D <- 0.0
  
  q_SG <- 0.05
  q_SD <- 0.0
  
  q_GS <- 0.05
  q_GD <- 0.05
  
  q_DS <- 0.0
  q_DG <- 0.05
  
  used_params <- c(spec_S, spec_G, spec_D,
                   ext_S, ext_G, ext_D,
                   q_SG, q_SD,
                   q_GS, q_GD,
                   q_DS, q_DG,
                   0, 0, 0, 0, 0, 0)
  
  sim_lambda_list_etd <- secsse::fill_in(idparslist_ETD[[1]], used_params)
  sim_mu_vector_etd   <- secsse::fill_in(idparslist_ETD[[2]], used_params)
  sim_q_matrix_etd    <- secsse::fill_in(idparslist_ETD[[3]], used_params)
  
  crown_age_used = 20
  
  num_repl <- 500 # increase if you want higher precision.
  
  found1 <- c()
  found2 <- c()
  
  for (r in 1:num_repl) {
    sim_tree <- secsse::secsse_sim(lambdas = sim_lambda_list_etd,
                                   mus = sim_mu_vector_etd,
                                   qs = sim_q_matrix_etd,
                                   crown_age = crown_age_used,
                                   num_concealed_states = 3,
                                   conditioning = "none",
                                   min_spec = 10,
                                   pool_init_states = c("S"),
                                   start_at_crown = FALSE)
    
    local_stats <- c(length(sim_tree$phy$tip.label))
    freq_1 <- sum(sim_tree$obs_traits == "S") / length(sim_tree$obs_traits)
    freq_2 <- sum(sim_tree$obs_traits == "G") / length(sim_tree$obs_traits)
    freq_3 <- sum(sim_tree$obs_traits == "D") / length(sim_tree$obs_traits)

    to_add <- c(unlist(local_stats), freq_1, freq_2, freq_3)
    found1 <- rbind(found1, to_add)
  }
  
  
  for (r in 1:num_repl) {
    sim_tree <- diversitree::tree.musse(pars = used_params[c(1:12)],
                                        max.t = crown_age_used,
                                        x0 = 1)
    
    while(length(sim_tree$tip.label) < 10) {
      sim_tree <- diversitree::tree.musse(pars = used_params[c(1:12)],
                                          max.t = crown_age_used,
                                          x0 = 1)
    }
    
    if (class(sim_tree) == "phylo") {
      local_stats <- c(length(sim_tree$tip.label))
      freq_1 <- sum(sim_tree$tip.state == 1) / length(sim_tree$tip.state)
      freq_2 <- sum(sim_tree$tip.state == 2) / length(sim_tree$tip.state)
      freq_3 <- sum(sim_tree$tip.state == 3) / length(sim_tree$tip.state)
      
      to_add <- c(unlist(local_stats), freq_1, freq_2, freq_3)
      found2 <- rbind(found2, to_add)
    }
  }
  
  for (i in 1:4) {
    a1 <- found1[, i]
    a2 <- found2[, i]
    b <- t.test(a1, a2)
    testthat::expect_true(b$p.value > 0.05)
  #  cat(i, mean(a1), mean(a2), b$p.value, "\n")
  }
})

test_that("test comparison bisse", {
  testthat::skip_on_ci()
  testthat::skip_on_cran()
  
  if (!requireNamespace("diversitree")) testthat::skip()
  
  states <- c("S", "G")
  
  spec_matrix <- c("S", "S", "S", 1)
  spec_matrix <- rbind(spec_matrix, c("G", "G", "G", 2))
  lambda_list <- secsse::create_lambda_list(state_names = states,
                                            num_concealed_states = 2,
                                            transition_matrix = spec_matrix,
                                            model = "ETD")
  
  mu_vector <- secsse::create_mu_vector(state_names = states,
                                        num_concealed_states = 2,
                                        model = "ETD",
                                        lambda_list = lambda_list)
  
  shift_matrix <- c("S", "G", 5)
  shift_matrix <- rbind(shift_matrix, c("G", "S", 6))
  
  q_matrix <- secsse::create_q_matrix(state_names = states,
                                      num_concealed_states = 2,
                                      shift_matrix = shift_matrix,
                                      diff.conceal = TRUE)
  idparslist <- list()
  idparslist[[1]] <- lambda_list
  idparslist[[2]] <- mu_vector
  idparslist[[3]] <- q_matrix
  
  idparslist_ETD <- idparslist
  
  spec_S <- 0.5
  spec_G <- 1.0
  ext_S <- ext_G <- 0.0
  q_SG <- 0.05
  q_GS <- 0.05
  used_params <- c(spec_S, spec_G, ext_S, ext_G, q_SG, q_GS, 0, 0)
  
  sim_lambda_list_etd <- secsse::fill_in(idparslist_ETD[[1]], used_params)
  sim_mu_vector_etd   <- secsse::fill_in(idparslist_ETD[[2]], used_params)
  sim_q_matrix_etd    <- secsse::fill_in(idparslist_ETD[[3]], used_params)
  
  crown_age_used = 5
  
  num_repl <- 100 # increase if you want higher precision.
  
  found1 <- c()
  found2 <- c()
  
  for (r in 1:num_repl) {
    sim_tree <- secsse::secsse_sim(lambdas = sim_lambda_list_etd,
                                   mus = sim_mu_vector_etd,
                                   qs = sim_q_matrix_etd,
                                   crown_age = crown_age_used,
                                   num_concealed_states = 2,
                                   conditioning = "none",
                                   min_spec = 10,
                                   pool_init_states = c("S"),
                                   start_at_crown = FALSE,
                                   verbose = TRUE)
    
    local_stats <- c(length(sim_tree$phy$tip.label))
    freq_1 <- sum(sim_tree$obs_traits == "S") / length(sim_tree$obs_traits)
    freq_2 <- sum(sim_tree$obs_traits == "G") / length(sim_tree$obs_traits)
    
    to_add <- c(unlist(local_stats), freq_1, freq_2)
    found1 <- rbind(found1, to_add)
  }
  
  
  for (r in 1:num_repl) {
    sim_tree <- diversitree::tree.musse(pars = used_params[1:6],
                                        max.t = crown_age_used,
                                        x0 = 1)
    
    while (length(sim_tree$tip.label) < 10) {
      sim_tree <- diversitree::tree.musse(pars = used_params[1:6],
                                          max.t = crown_age_used,
                                          x0 = 1)
    }
    
    if (class(sim_tree) == "phylo") {
      local_stats <- c(length(sim_tree$tip.label))
      freq_1 <- sum(sim_tree$tip.state == 1) / length(sim_tree$tip.state)
      freq_2 <- sum(sim_tree$tip.state == 2) / length(sim_tree$tip.state)
      
      to_add <- c(unlist(local_stats), freq_1, freq_2)
      found2 <- rbind(found2, to_add)
    }
  }
  
  for (i in 1:3) {
    a1 <- found1[, i]
    a2 <- found2[, i]
    b <- t.test(a1, a2)
    testthat::expect_true(b$p.value > 0.05)
 #   cat(i, mean(a1), mean(a2), b$p.value, "\n")
  }
})

test_that("test comparison classe", {
  
  states <- c("S", "G", "D")
  
  spec_matrix <- c("S", "S", "S", 1)
  spec_matrix <- rbind(spec_matrix, c("G", "G", "G", 2))
  spec_matrix <- rbind(spec_matrix, c("D", "D", "D", 3))
  
  spec_matrix <- rbind(spec_matrix, c("G", "S", "D", 4))
  
  lambda_list <- secsse::create_lambda_list(state_names = states,
                                            num_concealed_states = 3,
                                            transition_matrix = spec_matrix,
                                            model = "ETD")
  
  mu_vector <- secsse::create_mu_vector(state_names = states,
                                        num_concealed_states = 3,
                                        model = "ETD",
                                        lambda_list = lambda_list)
  
  shift_matrix <- c("S", "G", 8)
  shift_matrix <- rbind(shift_matrix, c( "S", "D", 9))
  
  shift_matrix <- rbind(shift_matrix, c("G", "S", 10))
  shift_matrix <- rbind(shift_matrix, c("G", "D", 11))
  
  shift_matrix <- rbind(shift_matrix, c( "D", "S", 12))
  shift_matrix <- rbind(shift_matrix, c( "D", "G", 13))
  
  q_matrix <- secsse::create_q_matrix(state_names = states,
                                      num_concealed_states = 3,
                                      shift_matrix = shift_matrix,
                                      diff.conceal = TRUE)
  idparslist <- list()
  idparslist[[1]] <- lambda_list
  idparslist[[2]] <- mu_vector
  idparslist[[3]] <- q_matrix
  
  idparslist_ETD <- idparslist
  
  spec_S <- 0.1
  spec_G <- 0.0
  spec_D <- 0.1
  
  spec_G_SD <- 0.2
  
  ext_S <- ext_G <- ext_D <- 0.0
  
  q_SG <- 0.0
  q_SD <- 0.0
  
  q_GS <- 0.0
  q_GD <- 0.0
  
  q_DS <- 0.0
  q_DG <- 0.0
  
  num_repl <- 300
  crown_age_used <- 30
  
  used_params <- c(spec_S, spec_G, spec_D, spec_G_SD,
                   ext_S, ext_G, ext_D,
                   q_SG, q_SD,
                   q_GS, q_GD,
                   q_DS, q_DG,
                   0, 0, 0, 0, 0, 0)
  
  sim_lambda_list_etd <- secsse::fill_in(idparslist_ETD[[1]], used_params)
  sim_mu_vector_etd   <- secsse::fill_in(idparslist_ETD[[2]], used_params)
  sim_q_matrix_etd    <- secsse::fill_in(idparslist_ETD[[3]], used_params)
  
  
  lambda_S <- c(spec_S, 0, 0, 0, 0, 0)         # lambda111 lambda112 lambda113 lambda122 lambda123 lambda133
  lambda_G <- c(0, 0, spec_G_SD, spec_G, 0, 0) # lambda211 lambda212 lambda213 lambda222 lambda223 lambda233 
  lambda_D <- c(0, 0, 0, 0, 0, spec_D)         # lambda311 lambda312 lambda313 lambda322 lambda323 lambda333
  
  
  params <- c(lambda_S, 
              lambda_G, 
              lambda_D,
              ext_S, ext_G, ext_D,
              q_SG, q_SD,
              q_GS, q_GD,
              q_DS, q_DG)
  
  found1 <- c()
  found2 <- c()
  
  for (r in 1:num_repl) {
    sim_tree <- secsse::secsse_sim(lambdas = sim_lambda_list_etd,
                                   mus = sim_mu_vector_etd,
                                   qs = sim_q_matrix_etd,
                                   crown_age = crown_age_used,
                                   num_concealed_states = 3,
                                   conditioning = "none",
                                   min_spec = 10,
                                   pool_init_states = c("G"),
                                   start_at_crown = FALSE)
    
    freq_1 <- sum(sim_tree$obs_traits == "S") / length(sim_tree$obs_traits)
    freq_2 <- sum(sim_tree$obs_traits == "G") / length(sim_tree$obs_traits)
    freq_3 <- sum(sim_tree$obs_traits == "D") / length(sim_tree$obs_traits)
    
    to_add <- c(freq_1, freq_2, freq_3)
    found1 <- rbind(found1, to_add)
  }
  
  for (r in 1:num_repl) {
    sim_tree <- diversitree::tree.classe(pars = params,
                                         max.t = crown_age_used,
                                         x0 = 2)
    
    while (length(sim_tree$tip.label) < 10) {
      sim_tree <- diversitree::tree.classe(pars = params,
                                           max.t = crown_age_used,
                                           x0 = 2)
    }
    
    if (class(sim_tree) == "phylo") {
      freq_1 <- sum(sim_tree$tip.state == 1) / length(sim_tree$tip.state)
      freq_2 <- sum(sim_tree$tip.state == 2) / length(sim_tree$tip.state)
      freq_3 <- sum(sim_tree$tip.state == 3) / length(sim_tree$tip.state)
    
      to_add <- c(freq_1, freq_2, freq_3)
      found2 <- rbind(found2, to_add)
    }
  }
  
  testthat::expect_true(max(found1[, 2]) == max(found2[, 2]))
  testthat::expect_true(max(found1[, 2]) == 0)
})
