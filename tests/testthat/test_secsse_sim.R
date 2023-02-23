context("test_secsse_sim")

test_that("compare R and Cpp simulation", {
  testthat::skip_on_cran()

  parenthesis <- "(((6:0.2547423371,(1:0.0496153503,4:0.0496153503):0.2051269868):0.1306304758,(9:0.2124135406,5:0.2124135406):0.1729592723):1.151205247,(((7:0.009347664296,3:0.009347664296):0.2101416075,10:0.2194892718):0.1035186448,(2:0.2575886319,8:0.2575886319):0.06541928469):1.213570144);"
  phylotree <- ape::read.tree(file = "",parenthesis)
  traits <- c(2, 0, 1, 0, 2, 0, 1, 2, 2, 0)
  num_concealed_states <- 3
  idparslist <- cla_id_paramPos(traits,num_concealed_states)
  idparslist$lambdas[2,] <- rep(1,9)
  idparslist[[2]][] <- 4
  masterBlock <- matrix(5,ncol = 3,nrow = 3,byrow = TRUE) 
  diag(masterBlock) <- NA
  diff.conceal <- FALSE
  idparslist[[3]] <- q_doubletrans(traits,masterBlock,diff.conceal)
  testthat::expect_output(
    startingpoint <- DDD::bd_ML(brts = ape::branching.times(phylotree))
  )
  intGuessLamba <- startingpoint$lambda0
  intGuessMu <- startingpoint$mu0
  idparsopt <- c(1)
  initparsopt <- c(rep(intGuessLamba,1))
  idparsfix <- c(0,4,5)
  parsfix <- c(0,0,0.01)
  tol <- c(1e-04, 1e-05, 1e-07)
  maxiter <- 1000 * round((1.25)^length(idparsopt))
  optimmethod <- "simplex"
  cond <- "proper_cond"
  root_state_weight <- "proper_weights"
  sampling_fraction <- c(1,1,1)
  
  testthat::expect_output(
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
  
  lambdas = model_R$MLpars[[1]]
  mus <- model_R$MLpars[[2]]
  maxSpec = 1000
  num_repl <- 100
  
  # now we simulate:
  for (max_time in c(1, 2)) {
    found_R <- c()
    for (r in 1:num_repl) {
      # and now we simulate:
      tree <- secsse::secsse_sim_cond(states = names(mus),
                                      lambdas = lambdas,
                                      mus = mus,
                                      timeSimul = max_time,
                                      qs = qs,
                                      pool_init_states = names(mus),
                                      maxSpec = maxSpec)
      if (length(tree$phy) > 1) {
        found_R[r] <- length(tree$phy$tip.label)
      }
    }
    
    found_cpp <- c()
    for (r in 1:num_repl) {
      # and now we simulate:
      tree <- secsse::secsse_sim_cond_cpp(states = names(mus),
                                          lambdas = lambdas,
                                          mus = mus,
                                          timeSimul = max_time,
                                          qs = qs,
                                          maxSpec = maxSpec)
      if (length(tree$phy) > 1) {
        found_cpp[r] <- length(tree$phy$tip.label)
      }
    }
  
    mean_r <- mean(found_R, na.rm = TRUE)
    mean_cpp <- mean(found_cpp, na.rm = TRUE)
    
    testthat::expect_equal(mean_r, mean_cpp, tolerance = 0.1)
  }
  
  # now with extinction:
  mus <- model_R$MLpars[[2]]
  old_mus <- mus
  mus2 <- rep(0.2, length(mus))
  names(mus2) <- names(mus)
  mus <- mus2
  # now we simulate:
  for (max_time in c(1, 2)) {
    found_R <- c()
    for (r in 1:num_repl) {
      # and now we simulate:
      tree <- secsse::secsse_sim_cond(states = names(mus),
                                      lambdas = lambdas,
                                      mus = mus,
                                      timeSimul = max_time,
                                      qs = qs,
                                      pool_init_states = names(mus),
                                      maxSpec = maxSpec)
      if (length(tree$phy) > 1) {
        found_R[r] <- length(tree$phy$tip.label)
      }
    }
    
    found_cpp <- c()
    for (r in 1:num_repl) {
      # and now we simulate:
      tree <- secsse::secsse_sim_cond_cpp(states = names(mus),
                                          lambdas = lambdas,
                                          mus = mus,
                                          timeSimul = max_time,
                                          qs = qs,
                                          maxSpec = maxSpec)
      if (length(tree$phy) > 1) {
        found_cpp[r] <- length(tree$phy$tip.label)
      }
    }
    
    mean_r <- mean(found_R, na.rm = TRUE)
    mean_cpp <- mean(found_cpp, na.rm = TRUE)
    
    testthat::expect_equal(mean_r, mean_cpp, tolerance = 0.1)
  }
})