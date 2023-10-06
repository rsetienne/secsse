test_that("trying a short ML search: cla_secsse", {
  Sys.unsetenv("R_TESTS")
  parenthesis <- "(((6:0.2547423371,(1:0.0496153503,4:0.0496153503):0.2051269868):0.1306304758,(9:0.2124135406,5:0.2124135406):0.1729592723):1.151205247,(((7:0.009347664296,3:0.009347664296):0.2101416075,10:0.2194892718):0.1035186448,(2:0.2575886319,8:0.2575886319):0.06541928469):1.213570144);" #nolint
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
  maxiter <- 1000 * round((1.25) ^ length(idparsopt))
  optimmethod <- "subplex"
  cond <- "proper_cond"
  root_state_weight <- "proper_weights"
  sampling_fraction <- c(1, 1, 1)

  # Expect warning because some transitions are set to be impossible
  testthat::expect_warning(
    model_R <- cla_secsse_ml(
      phy = phylotree,
      traits = traits,
      num_concealed_states = num_concealed_states,
      idparslist = idparslist,
      idparsopt = idparsopt,
      initparsopt = initparsopt,
      idparsfix = idparsfix,
      parsfix = parsfix,
      cond = cond,
      root_state_weight = root_state_weight,
      sampling_fraction = sampling_fraction,
      tol = tol,
      maxiter = maxiter,
      optimmethod = optimmethod,
      num_cycles = 1,
      verbose = FALSE)
  )

  testthat::expect_equal(model_R$ML, -16.1342246206186)
  
  # we have to translate to lambda matrices to test the following:
  param_posit <- idparslist
  param_posit[[1]] <- secsse::prepare_full_lambdas(traits,
                                                   num_concealed_states,
                                                   idparslist[[1]])
  
  found_pars <- secsse::extract_par_vals(param_posit,
                                         model_R$MLpars)
  testthat::expect_equal(length(found_pars),
                         max(param_posit[[3]], na.rm = TRUE))
})
