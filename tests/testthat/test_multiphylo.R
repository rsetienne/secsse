test_that("multi phylo", {
  set.seed(42)

  focal_tree <- ape::rphylo(n = 3, birth = 0.3 ,death = 0)

  focal_tree$root.edge <- NULL
  traits <- c(1, 1, 1)
  
  num_concealed_states <- 2
  idparslist <- cla_id_paramPos(c(1, 2), num_concealed_states)
  idparslist$lambdas[1, ] <- rep(1, 2)
  idparslist[[2]][] <- 2
  masterBlock <- matrix(3, ncol = 2, nrow = 2, byrow = TRUE)
  diag(masterBlock) <- NA
  diff.conceal <- FALSE
  idparslist[[3]] <- q_doubletrans(c(1, 2), masterBlock, diff.conceal)
  idparslist[[1]] <- secsse::prepare_full_lambdas(c(1, 2),
                                                  num_concealed_states,
                                                  idparslist[[1]])
  
  # Expect warning because some transitions are set to be impossible
  
  params <- c(0.3, 0.0, 0.0) # extinction and shifts to zero, to allow direct
  # comparison
  
  lambdas <- secsse::fill_in(idparslist[[1]], params)
  mus <- secsse::fill_in(idparslist[[2]], params)
  q_mat <- secsse::fill_in(idparslist[[3]], params)
  
  parslist <- list()
  parslist[[1]] <- lambdas
  parslist[[2]] <- mus
  parslist[[3]] <- q_mat
  
  sf <- c(1, 1)
  
  focal_tree$root.edge <- NULL
  testthat::expect_warning(
  res1 <- secsse::cla_secsse_loglik(parameter = parslist,
                                    phy = focal_tree,
                                    traits = traits,
                                    num_concealed_states = num_concealed_states,
                                    sampling_fraction = sf,
                                    cond = "no_cond")
  )
  trees <- list()
  trees[[1]]<- focal_tree
  trees[[2]]<- focal_tree
  
  class(trees) <- "multiPhylo"
  
  trait_list <- list()
  trait_list[[1]] <- traits
  trait_list[[2]] <- traits
  
  testthat::expect_warning(
  res2 <- secsse::cla_secsse_loglik(parameter = parslist,
                                    phy = trees,
                                    traits = trait_list,
                                    num_concealed_states = num_concealed_states,
                                    sampling_fraction = sf,
                                    cond = "no_cond")
  )
  testthat::expect_equal(2*res1, res2)
  
  sf_list <- list()
  sf_list[[1]] <- sf
  sf_list[[2]] <- sf
  
  testthat::expect_warning(
    res3 <- secsse::cla_secsse_loglik(parameter = parslist,
                                      phy = trees,
                                      traits = trait_list,
                                      num_concealed_states = num_concealed_states,
                                      sampling_fraction = sf_list,
                                      cond = "no_cond")
  )
  testthat::expect_equal(2*res1, res3)
})

test_that("multi phylo ML", {
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
  
  # and now let's do multi phylo stuff!
  
  phylo_list <- list()
  trait_list <- list()
  sf_list    <- list()
  
  for (r in 1:3) {
    phylo_list[[r]] <- phylotree
    trait_list[[r]] <- traits
    sf_list[[r]] <- sampling_fraction
  }
  
  class(phylo_list) <- "multiPhylo"
  testthat::expect_warning(
  multi_R <- cla_secsse_ml(
    phy = phylo_list,
    traits = trait_list,
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
  
  testthat::expect_equal(3 * model_R$ML, multi_R$ML)
  testthat::expect_true(all.equal(model_R$MLpars, multi_R$MLpars))
  
  # repeat with list:
  testthat::expect_warning(
    multi_R <- cla_secsse_ml(
      phy = phylo_list,
      traits = trait_list,
      num_concealed_states = num_concealed_states,
      idparslist = idparslist,
      idparsopt = idparsopt,
      initparsopt = initparsopt,
      idparsfix = idparsfix,
      parsfix = parsfix,
      cond = cond,
      root_state_weight = root_state_weight,
      sampling_fraction = sf_list,
      tol = tol,
      maxiter = maxiter,
      optimmethod = optimmethod,
      num_cycles = 1,
      verbose = FALSE)
  )
  
  testthat::expect_equal(3 * model_R$ML, multi_R$ML)
  testthat::expect_true(all.equal(model_R$MLpars, multi_R$MLpars))
})