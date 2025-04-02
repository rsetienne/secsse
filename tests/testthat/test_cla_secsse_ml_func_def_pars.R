test_that("multiplication works", {
  
  set.seed(16)
  phylotree <- ape::rbdtree(0.07,0.001,Tmax=50)
  startingpoint <- testthat::expect_output(
    DDD::bd_ML(brts = ape::branching.times(phylotree))
  )
  intGuessLamba <- startingpoint$lambda0
  intGuessMu <- startingpoint$mu0
  traits <-  sample(c(0,1,2),
                    ape::Ntip(phylotree), replace = TRUE) # get some traits
  num_concealed_states <- 3
  idparslist <- cla_id_paramPos(traits, num_concealed_states)
  idparslist$lambdas[1,] <- c(1,2,3,1,2,3,1,2,3)
  idparslist[[2]][] <- 4
  masterBlock <- matrix(c(5,6,5,6,5,6,5,6,5),ncol = 3, nrow=3, byrow = TRUE)
  diag(masterBlock) <- NA
  diff.conceal <- FALSE
  idparslist[[3]] <- q_doubletrans(traits,masterBlock,diff.conceal)
  idparsfuncdefpar <- c(3,5,6)
  idparsopt <- c(1,2)
  idparsfix <- c(0,4)
  initparsopt <- c(rep(intGuessLamba,2))
  parsfix <- c(0,0)
  idfactorsopt <- 1
  initfactors <- 4
  
  functions_defining_params <- list()
  functions_defining_params[[1]] <- function() {
    par_3 <- par_1 + par_2
  }
  functions_defining_params[[2]] <- function() {
    par_5 <- par_1 * factor_1
  }
  functions_defining_params[[3]] <- function() {
    par_6 <- par_3 * factor_1
  }
  
  tol = c(1e-02, 1e-03, 1e-04)
  maxiter = 100 * round((1.25)^length(idparsopt))
  optimmethod = 'subplex'
  cond <- 'proper_cond'
  root_state_weight <- 'proper_weights'
  sampling_fraction <- c(1, 1, 1)
  model <- testthat::expect_warning(
    cla_secsse_ml_func_def_pars(
    phylotree,
    traits,
    num_concealed_states,
    idparslist,
    idparsopt,
    initparsopt,
    idfactorsopt,
    initfactors,
    idparsfix,
    parsfix,
    idparsfuncdefpar,
    functions_defining_params,
    cond,
    root_state_weight,
    sampling_fraction,
    tol,
    maxiter,
    optimmethod,
    num_cycles = 1,
    verbose = 0
  ))
  
  testthat::expect_equal(model$ML, -136.4534, tol = 1e-4)
  testthat::expect_length(model, 3)
  testthat::expect_length(model$MLpars, 3)
  testthat::expect_equal(model$MLpars[[2]],
               c("0A" = 0,
                 "1A" = 0,
                 "2A" = 0,
                 "0B" = 0,
                 "1B" = 0,
                 "2B" = 0,
                 "0C" = 0,
                 "1C" = 0,
                 "2C" = 0
               )
  )
})