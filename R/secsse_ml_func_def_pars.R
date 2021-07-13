#' Maximum likehood estimation under Several examined and concealed States-dependent Speciation and Extinction (SecSSE) where some paramaters are functions of other parameters and/or factors.
#' @title Maximum likehood estimation for (SecSSE) with parameter as complex functions.
#' @param phy phylogenetic tree of class phylo, ultrametric, rooted and with branch lengths.
#' @param traits a vector with trait states for each tip in the phylogeny.
#' @param num_concealed_states number of concealed states, generally equivalent to the number of examined states in the dataset.
#' @param idparslist overview of parameters and their values.
#' @param idparsopt id of parameters to be estimated.
#' @param initparsopt initial guess of the parameters to be estimated.
#' @param idfactorsopt id of the factors that will be optimized. There are not fixed factors, so use a constant within 'functions_defining_params'.
#' @param initfactors the initial guess for a factor (it should be set to NULL when no factors).
#' @param idparsfix id of the fixed parameters (it should be set to NULL when there are no factors).
#' @param parsfix value of the fixed parameters.
#' @param idparsfuncdefpar id of the parameters which will be a function of optimized and/or fixed parameters. The order of id should match functions_defining_params
#' @param functions_defining_params a list of functions. Each element will be a function which defines a parameter e.g. id_3 <- (id_1+id_2)/2. See example and vigenette
#' @param cond condition on the existence of a node root: "maddison_cond","proper_cond"(default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:"maddison_weights","proper_weights"(default) or "equal_weights". It can also be specified the root state:the vector c(1,0,0) indicates state 1 was the root state.
#' @param sampling_fraction vector that states the sampling proportion per trait state. It must have as many elements as there are trait states.
#' @param tol maximum tolerance. Default is "c(1e-04, 1e-05, 1e-05)".
#' @param maxiter max number of iterations. Default is "1000 *round((1.25)^length(idparsopt))".
#' @param use_fortran Should the Fortran code for numerical integration be called? Default is TRUE.
#' @param methode method used for integration calculation. Default is "ode45".
#' @param optimmethod method used for optimization. Default is "simplex".
#' @param num_cycles number of cycles of the optimization (default is 1).
#' @param run_parallel should the routine to run in parallel be called? Read note below
#' @param loglik_penalty the size of the penalty for all parameters; default is 0 (no penalty)
#' @param is_complete_tree whether or not a tree with all its extinct species is provided
#' @param func function to be used in solving the ODE system. Currently only for testing purposes.
#' @note To run in parallel it is needed to load the following libraries when windows: apTreeshape, doparallel and foreach. When unix, it requires: apTreeshape, doparallel, foreach and doMC
#' @return Parameter estimated and maximum likelihood
#' @examples
#'# Example of how to set the arguments for a ML search.
#'rm(list=ls(all=TRUE))
#'library(secsse)
#'library(DDD)
#'set.seed(16)
#'phylotree <- ape::rbdtree(0.07,0.001,Tmax=50)
#'startingpoint<-bd_ML(brts = ape::branching.times(phylotree))
#'intGuessLamba <- startingpoint$lambda0
#'intGuessMu <- startingpoint$mu0
#'traits <- sample(c(0,1,2), ape::Ntip(phylotree),replace=TRUE) #get some traits
#'num_concealed_states<-3
#'idparslist<-id_paramPos(traits, num_concealed_states)
#'idparslist[[1]][c(1,4,7)] <- 1
#'idparslist[[1]][c(2,5,8)] <- 2
#'idparslist[[1]][c(3,6,9)] <- 3
#'idparslist[[2]][] <- 4
#'masterBlock <- matrix(c(5,6,5,6,5,6,5,6,5),ncol = 3,nrow = 3,byrow = TRUE)
#'diag(masterBlock) <- NA
#'diff.conceal <- FALSE
#'idparslist[[3]] <- q_doubletrans(traits,masterBlock,diff.conceal)
#'idparsfuncdefpar <- c(3,5,6)
#'idparsopt <- c(1,2)
#'idparsfix <- c(0,4)
#'initparsopt <- c(rep(intGuessLamba,2))
#'parsfix <- c(0,0)
#'idfactorsopt <- 1
#'initfactors <- 4
#'# functions_defining_params is a list of functions. Each function has no arguments and to refer
#'# to parameters ids should be indicated as "par_" i.e. par_3 refers to parameter 3. When a
#'# function is defined, be sure that all the parameters involved are either estimated, fixed or
#'# defined by previous functions (i.e, a function that defines parameter in
#'# 'functions_defining_params'). The user is responsible for this. In this example, par_3
#'# (i.e., parameter 3) is needed to calculate par_6. This is correct because par_3 is defined in
#'# the first function of 'functions_defining_params'. Notice that factor_1 indicates a value
#'# that will be estimated to satisfy the equation. The same factor can be shared to define
#'# several parameters.
#'functions_defining_params <- list()
#'functions_defining_params[[1]] <- function(){
#'  par_3 <- par_1 + par_2
#'}
#'functions_defining_params[[2]] <- function(){
#'  par_5 <- par_1 * factor_1
#'}
#'functions_defining_params[[3]] <- function(){
#'  par_6 <- par_3 * factor_1
#'}

#'
#'tol = c(1e-04, 1e-05, 1e-07)
#'maxiter = 1000 * round((1.25)^length(idparsopt))
#'use_fortran = TRUE
#'methode = "ode45"
#'optimmethod = "simplex"
#'run_parallel = FALSE
#'cond<-"proper_cond"
#'root_state_weight <- "proper_weights"
#'sampling_fraction <- c(1,1,1)
#'#model <- secsse_ml_func_def_pars(phylotree, traits, num_concealed_states, idparslist, idparsopt,
#'#                               initparsopt, idfactorsopt, initfactors, idparsfix, parsfix,
#'#                               idparsfuncdefpar, functions_defining_params, cond,
#'#                               root_state_weight, sampling_fraction, tol, maxiter, use_fortran,
#'#                               methode, optimmethod, num_cycles = 1, run_parallel)
#'
#'# ML -136.5796
#' @export

secsse_ml_func_def_pars <- function(phy,
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
                                    cond = "proper_cond",
                                    root_state_weight = "proper_weights",
                                    sampling_fraction,
                                    tol = c(1E-4, 1E-5, 1E-7),
                                    maxiter = 1000 * round((1.25) ^ length(idparsopt)),
                                    use_fortran = TRUE,
                                    methode = "ode45",
                                    optimmethod = 'simplex',
                                    num_cycles = 1,
                                    run_parallel = FALSE,
                                    loglik_penalty = 0,
                                    is_complete_tree = FALSE,
                                    func = ifelse(is_complete_tree,
                                                  "secsse_runmod_ct",
                                                  ifelse(use_fortran == FALSE,
                                                         secsse_loglik_rhs,
                                                         "secsse_runmod2")
                                    )
) {
  
  structure_func <- list()
  structure_func[[1]] <- idparsfuncdefpar
  structure_func[[2]] <- functions_defining_params
  if(is.null(idfactorsopt)){
    structure_func[[3]] <- "noFactor"
  } else {
    structure_func[[3]] <- idfactorsopt
  }
  
  see_ancestral_states <- FALSE
  if (is.null(idfactorsopt) == FALSE) {
    if (length(initfactors) != length(idfactorsopt)) {
      stop("idfactorsopt should have the same length as initfactors.")
    }
  }
  
  if (is.list(functions_defining_params) == FALSE) {
    stop(
      "The argument functions_defining_params should be a list of functions. See example and vignette"
    )
  }
  
  if (length(functions_defining_params) != length(idparsfuncdefpar)) {
    stop(
      "The argument functions_defining_params should have the same length than idparsfuncdefpar"
    )
  }
  
  if (is.matrix(traits)) {
    cat("You are setting a model where some species had more than one trait state \n")
  }
  
  if (length(initparsopt) != length(idparsopt)) {
    stop(
      "initparsopt must be the same length as idparsopt. Number of parameters to optimize does not match the number of initial values for the search"
    )
  }
  
  if (length(idparsfix) != length(parsfix)) {
    stop(
      "idparsfix and parsfix must be the same length. Number of fixed elements does not match the fixed figures"
    )
  }
  
  if (anyDuplicated(c(idparsopt, idparsfix, idparsfuncdefpar)) != 0) {
    stop("At least one element was asked to be fixed, estimated or a function at the same time")
  }
  
  if (identical(as.numeric(sort(
    c(idparsopt, idparsfix, idparsfuncdefpar)
  )), as.numeric(sort(unique(
    unlist(idparslist)
  )))) == FALSE) {
    stop(
      "All elements in idparslist must be included in either idparsopt or idparsfix or idparsfuncdefpar "
    )
  }
  
  if (anyDuplicated(c(unique(sort(
    as.vector(idparslist[[3]])
  )), idparsfix[which(parsfix == 0)])) != 0) {
    cat("You set some transitions as impossible to happen", "\n")
  }
  
  
  #options(warn = -1)
  cat("Calculating the likelihood for the initial parameters.", "\n")
  utils::flush.console()
  
  initparsopt2 <- c(initparsopt, initfactors)
  
  trparsopt <- initparsopt2 / (1 + initparsopt2)
  trparsopt[which(initparsopt2 == Inf)] <- 1
  trparsfix <- parsfix / (1 + parsfix)
  trparsfix[which(parsfix == Inf)] <- 1
  
  mus <- calc_mus(is_complete_tree, idparslist, idparsfix, parsfix, idparsopt, initparsopt)
  
  optimpars <- c(tol, maxiter)

  if(run_parallel == TRUE){
    setting_calculation <-
      build_initStates_time_bigtree(phy, traits, num_concealed_states, sampling_fraction, is_complete_tree, mus)
    setting_parallel <- 1
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(2)
      doParallel::registerDoParallel(cl)
      # pass libPath to workers
      # see https://stackoverflow.com/questions/6412459/how-to-specify-the-location-of-r-packages-in-foreach-packages-pkg-do
      # https://gitlab.com/CarlBrunius/MUVR/-/issues/11
      parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
      on.exit(parallel::stopCluster(cl))
    }
    if (.Platform$OS.type == "unix") {
      doMC::registerDoMC(2)
    }
  } else {
    setting_calculation <-
      build_initStates_time(phy, traits, num_concealed_states, sampling_fraction, is_complete_tree, mus)
    setting_parallel <- NULL
  }
  
  if(optimmethod == 'subplex') {verbose <- TRUE} else {verbose <- FALSE}
  initloglik <-
    secsse_loglik_choosepar(
      trparsopt = trparsopt,
      trparsfix = trparsfix,
      idparsopt = idparsopt,
      idparsfix = idparsfix,
      idparslist = idparslist,
      structure_func = structure_func,
      phy = phy,
      traits = traits,
      num_concealed_states =
        num_concealed_states,
      use_fortran = use_fortran,
      methode = methode,
      cond = cond,
      root_state_weight = root_state_weight,
      sampling_fraction = sampling_fraction,
      setting_calculation =
        setting_calculation,
      run_parallel = run_parallel,
      setting_parallel = setting_parallel,
      see_ancestral_states = see_ancestral_states,
      loglik_penalty = loglik_penalty,
      is_complete_tree = is_complete_tree,
      func = func,
      verbose = verbose
    )
  cat("The loglikelihood for the initial parameter values is",
      initloglik,
      "\n")
  if (initloglik == -Inf)
  {
    stop(
      "The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values."
    )
  } else {
    cat("Optimizing the likelihood - this may take a while.", "\n")
    utils::flush.console()
    cat(setting_parallel, "\n")
    out <-
      DDD::optimizer(
        optimmethod = optimmethod,
        optimpars = optimpars,
        fun = secsse_loglik_choosepar,
        trparsopt = trparsopt,
        idparsopt = idparsopt,
        trparsfix = trparsfix,
        idparsfix = idparsfix,
        idparslist = idparslist,
        structure_func = structure_func,
        phy = phy,
        traits = traits,
        num_concealed_states = num_concealed_states,
        use_fortran = use_fortran,
        methode = methode,
        cond = cond,
        root_state_weight = root_state_weight,
        sampling_fraction = sampling_fraction,
        setting_calculation = setting_calculation,
        run_parallel = run_parallel,
        setting_parallel = setting_parallel,
        see_ancestral_states = see_ancestral_states,
        num_cycles = num_cycles,
        loglik_penalty = loglik_penalty,
        is_complete_tree = is_complete_tree,
        func = func,
        verbose = verbose
      )
    if (out$conv != 0)
    {
      stop("Optimization has not converged. Try again with different initial values.\n")
    } else {
      MLpars1 <-
        secsse_transform_parameters(
          as.numeric(unlist(out$par)),
          trparsfix,
          idparsopt,
          idparsfix,
          idparslist,
          structure_func
        )
      out2 <-
        list(MLpars = MLpars1, ML = as.numeric(unlist(out$fvalues)),conv = out$conv)
    }
  }
  return(out2)
}
