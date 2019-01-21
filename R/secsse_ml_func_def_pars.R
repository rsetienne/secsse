secsse_transform_parameters_funcdefpar <-
  function(trparsopt,
           trparsfix,
           idparsopt,
           idparsfix,
           idparslist,
           idparsfuncdefpar,
           functions_defining_params,
           idfactosopt) {
    trparfuncdefpar <- NULL
    ids_all <- c(idparsfix, idparsopt)
    
    values_all <- c(trparsfix / (1 - trparsfix), trparsopt / (1 - trparsopt))
    
    for (jj in 1:length(functions_defining_params)) {
      myfunc <- functions_defining_params[[jj]]
      
      calculate_param_function <-
        function(myfunc,
                 ids_all,
                 values_all,
                 idfactosopt) {
          x <- as.list(values_all) ## To declare all the ids as variables
          
          if (is.null(idfactosopt)) {
            names(x) <- paste0("par_", ids_all)
          } else {
            names(x) <-
              c(paste0("par_", ids_all),
                paste0("factor_", idfactosopt))
          }
          list2env(x, envir = .GlobalEnv)
          myfunc()
        }
      value_func_defining_parm <-
        calculate_param_function(myfunc, ids_all, values_all, idfactosopt)
      
      ## Now, declare the variable that is just calculated, so it is
      ## available for the next calculation if needed
      y <- as.list(value_func_defining_parm)
      names(y) <- paste0("par_", idparsfuncdefpar[jj])
      list2env(y, envir = .GlobalEnv)
      
      if (is.numeric(value_func_defining_parm) == FALSE) {
        stop("something went with the calculation of parameters in 'functions_param_struct'")
      }
      trparfuncdefpar <- c(trparfuncdefpar, value_func_defining_parm)
    }
    trparfuncdefpar <- trparfuncdefpar / (1 + trparfuncdefpar)
    
    trpars1 <- idparslist
    for (j in 1:3) {
      trpars1[[j]][] = NA
    }
    if (length(idparsfix) != 0) {
      for (i in 1:length(idparsfix)) {
        for (j in 1:3) {
          id <- which(idparslist[[j]] == idparsfix[i])
          trpars1[[j]][id] <- trparsfix[i]
          
        }
      }
    }
    for (i in 1:length(idparsopt)) {
      for (j in 1:3)
      {
        id <- which(idparslist[[j]] == idparsopt[i])
        trpars1[[j]][id] <- trparsopt[i]
        
      }
    }
    for (i in 1:length(idparsfuncdefpar)) {
      ## structure part
      for (j in 1:3)
      {
        id <- which(idparslist[[j]] == idparsfuncdefpar[i])
        trpars1[[j]][id] <- trparfuncdefpar[i]
        
      }
    }
    pars1 <- list()
    for (j in 1:3) {
      pars1[[j]] <- trpars1[[j]] / (1 - trpars1[[j]])
    }
    return(pars1)
  }

secsse_loglik_choosepar_funcdefpar <-
  function(trparsopt,
           trparsfix,
           idparsopt,
           idparsfix,
           idparslist,
           idparsfuncdefpar = idparsfuncdefpar,
           functions_defining_params = functions_defining_params,
           idfactosopt = idfactosopt,
           phy = phy,
           traits = traits,
           num_concealed_states = num_concealed_states,
           use_fortran = use_fortran,
           methode,
           cond = cond,
           root_state_weight = root_state_weight,
           sampling_fraction = sampling_fraction,
           setting_calculation = setting_calculation,
           run_parallel = run_parallel,
           setting_parallel = setting_parallel) {
    alltrpars <- c(trparsopt, trparsfix)
    if (max(alltrpars) > 1 | min(alltrpars) < 0) {
      loglik = -Inf
    } else {
      pars1 <-
        secsse_transform_parameters_funcdefpar(
          trparsopt,
          trparsfix,
          idparsopt,
          idparsfix,
          idparslist,
          idparsfuncdefpar,
          functions_defining_params,
          idfactosopt
        )
      
      loglik <-
        secsse_loglik(
          parameter = pars1,
          phy = phy,
          traits = traits,
          num_concealed_states = num_concealed_states,
          use_fortran = use_fortran,
          methode = methode,
          cond = cond,
          root_state_weight = root_state_weight,
          sampling_fraction = sampling_fraction,
          run_parallel = run_parallel,
          setting_calculation = setting_calculation,
          setting_parallel = setting_parallel
        )
      
      if (is.nan(loglik) || is.na(loglik)) {
        print(trparsopt) ## new thing
        cat("There are parameter values used which cause numerical problems.\n")
        loglik <- -Inf
      }
    }
    return(loglik)
  }

#' Maximum likehood estimation under Several examined and concealed States-dependent Speciation and Extinction (SecSSE) where some paramaters are functions of other parameters and/or factors.
#' @title Maximum likehood estimation for (SecSSE) with parameter as complex functions.
#' @param phy phylogenetic tree of class phylo, ultrametric, rooted and with branch lengths.
#' @param traits a vector with trait states for each tip in the phylogeny.
#' @param num_concealed_states number of concealed states, generally equivalent to the number of examined states in the dataset.
#' @param idparslist overview of parameters and their values.
#' @param idparsopt id of parameters to be estimated.
#' @param initparsopt initial guess of the parameters to be estimated.
#' @param idfactosopt id of the factors that will be optimized. There are not fixed factors, so use a constant within 'functions_defining_params'.
#' @param initfactos the initial guess for a factor (it should be set to NULL when no factors).
#' @param idparsfix id of the fixed parameters (it should be set to NULL when no factors).
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
#' @param run_parallel should the routine to run in parallel be called? Read note below
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
#'traits <-  sample(c(0,1,2), ape::Ntip(phylotree),replace=TRUE) #get some traits
#'num_concealed_states<-3
#'idparslist<-id_paramPos(traits, num_concealed_states)
#'idparslist[[1]][c(1,4,7)]<-1
#'idparslist[[1]][c(2,5,8)]<-2
#'idparslist[[1]][c(3,6,9)]<-3
#'idparslist[[2]][]<-4
#'masterBlock<-matrix(c(5,6),ncol=3,nrow=3,byrow=TRUE)
#'diag(masterBlock)<-NA
#'diff.conceal <- FALSE
#'idparslist[[3]]<-q_doubletrans(traits,masterBlock,diff.conceal)
#'idparsfuncdefpar<-c(3,5,6)
#'idparsopt<-c(1,2)
#'idparsfix<-c(0,4)
#'initparsopt<-c(rep(intGuessLamba,2))
#'parsfix<-c(0,0)
#'idfactosopt<-1
#'initfactos<- 4
#'# functions_defining_params is a list of functions. Each function has no arguments
#'# and referring to parameters ids should be indicated as "par_", i.e. par_3 refers to
#'# parameter 3. When a function is defined, all parameters involved should be either
#'# estimated, fixed or defined by previous functions (i.e, a function that
#'# defines parameter in 'functions_defining_params').
#'# In this example, par_3 (i.e., parameter 3) is needed to calculate par_6. This is correct
#'# because par_3 is defined in the first function of 'functions_defining_params'.
#'# Note that factor_1 indicates a value that will be estimated to satisfy the equation.
#'# The same factor can be shared to define several parameters.
#'functions_defining_params<-list()
#'functions_defining_params[[1]]<-function(){
#'  par_3 <- par_1 + par_2
#'}
#'functions_defining_params[[2]]<-function(){
#'  par_5 <- par_1 * factor_1
#'}
#'functions_defining_params[[3]]<-function(){
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
#'sampling_fraction<-c(1,1,1)
#'#model<-secsse_ml_func_def_pars(
#'#phylotree, traits,
#'#num_concealed_states,
#'#idparslist,
#'#idparsopt,
#'#initparsopt,
#'#idfactosopt,
#'#initfactos,
#'#idparsfix,
#'#parsfix,
#'#idparsfuncdefpar,
#'#functions_defining_params,
#'#cond,
#'#root_state_weight,
#'#sampling_fraction,
#'#tol,
#'#maxiter,
#'#use_fortran,
#'#methode,
#'#optimmethod,
#'#run_parallel)
#'
#'# ML -136.5796
#' @export

secsse_ml_func_def_pars <- function(phy,
                                    traits,
                                    num_concealed_states,
                                    idparslist,
                                    idparsopt,
                                    initparsopt,
                                    idfactosopt,
                                    initfactos,
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
                                    run_parallel = FALSE) {
  if (is.null(idfactosopt) == FALSE) {
    if (length(initfactos) != length(idfactosopt)) {
      stop("idfactosopt should have the same length than initfactos.")
    }
  }
  
  if (is.list(functions_defining_params) == FALSE) {
    stop(
      "the argument functions_defining_params should be a list of functions. See example and vignette"
    )
  }
  
  if (length(functions_defining_params) != length(idparsfuncdefpar)) {
    stop(
      "the argument functions_defining_params should have the same length than idparsfuncdefpar"
    )
  }
  
  if (is.matrix(traits)) {
    cat("you are setting a model where some species had more than one trait state \n")
  }
  
  if (length(initparsopt) != length(idparsopt)) {
    stop(
      "initparsopt must be the same length as idparsopt. Number of parameters to optimize does not match the number of initial values for the search"
    )
  }
  
  if (length(idparsfix) != length(parsfix)) {
    stop(
      "idparsfix and parsfix must be the same length.Number of fixed elements does not match the fixed figures"
    )
  }
  
  if (anyDuplicated(c(idparsopt, idparsfix, idparsfuncdefpar)) != 0) {
    stop("at least one element was asked to be fixed, estimated or a function at the same time")
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
    cat("You set some transition states as non possible to happen", "\n")
  }
  
  
  options(warn = -1)
  cat("Calculating the likelihood for the initial parameters.", "\n")
  utils::flush.console()
  
  initparsopt2 <- c(initparsopt, initfactos)
  
  trparsopt <- initparsopt2 / (1 + initparsopt2)
  trparsopt[which(initparsopt2 == Inf)] <- 1
  trparsfix <- parsfix / (1 + parsfix)
  trparsfix[which(parsfix == Inf)] <- 1
  
  
  optimpars <- c(tol, maxiter)
  
  
  
  if (.Platform$OS.type == "windows" && run_parallel == TRUE) {
    cl <- parallel::makeCluster(2)
    doParallel::registerDoParallel(cl)
    setting_calculation <-
      build_initStates_time_bigtree(phy, traits, num_concealed_states, sampling_fraction)
    setting_parallel <- 1
    on.exit(parallel::stopCluster(cl))
  }
  
  if (.Platform$OS.type == "unix" && run_parallel == TRUE) {
    doMC::registerDoMC(2)
    setting_calculation <-
      build_initStates_time_bigtree(phy, traits, num_concealed_states, sampling_fraction)
    setting_parallel <- 1
  }
  
  if (run_parallel == FALSE) {
    setting_calculation <-
      build_initStates_time(phy, traits, num_concealed_states, sampling_fraction)
    setting_parallel <- NULL
  }
  
  
  initloglik <-
    secsse_loglik_choosepar_funcdefpar(
      trparsopt = trparsopt,
      trparsfix = trparsfix,
      idparsopt = idparsopt,
      idparsfix = idparsfix,
      idparslist = idparslist,
      idparsfuncdefpar = idparsfuncdefpar,
      functions_defining_params = functions_defining_params,
      idfactosopt = idfactosopt,
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
      setting_parallel = setting_parallel
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
        fun = secsse_loglik_choosepar_funcdefpar,
        trparsopt = trparsopt,
        idparsopt = idparsopt,
        trparsfix = trparsfix,
        idparsfix = idparsfix,
        idparslist = idparslist,
        idparsfuncdefpar = idparsfuncdefpar ,
        functions_defining_params = functions_defining_params,
        idfactosopt = idfactosopt,
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
        setting_parallel = setting_parallel
      )
    if (out$conv != 0)
    {
      stop("Optimization has not converged. Try again with different initial values.\n")
    } else {
      MLpars1 <-
        secsse_transform_parameters_funcdefpar(
          as.numeric(unlist(out$par)),
          trparsfix,
          idparsopt,
          idparsfix,
          idparslist,
          idparsfuncdefpar,
          functions_defining_params,
          idfactosopt
        )
      out2 <-
        list(MLpars = MLpars1, ML = as.numeric(unlist(out$fvalues)))
    }
  }
  return(out2)
}
