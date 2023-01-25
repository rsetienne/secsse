#' Maximum likehood estimation under Several examined and concealed 
#' States-dependent Speciation and Extinction (SecSSE) with cladogenetic option
#' @title Maximum likehood estimation for (SecSSE)
#' @param phy phylogenetic tree of class phylo, ultrametric, rooted and with 
#' branch lengths.
#' @param traits a vector with trait states for each tip in the phylogeny.
#' @param num_concealed_states number of concealed states, generally equivalent 
#' to the number of examined states in the dataset.
#' @param idparslist overview of parameters and their values.
#' @param idparsopt id of parameters to be estimated.
#' @param initparsopt initial guess of the parameters to be estimated.
#' @param idparsfix id of the fixed parameters.
#' @param parsfix value of the fixed parameters.
#' @param critical_t vector of time points at which a new parameter regime is
#' used, if applicable.
#' @param cond condition on the existence of a node root: 'maddison_cond',
#' 'proper_cond'(default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:'maddison_weights','proper_weights'(default) or 'equal_weights'. It can also be specified the 
#' root state:the vector c(1,0,0) indicates state 1 was the root state.
#' @param sampling_fraction vector that states the sampling proportion per 
#' trait state. It must have as many elements as there are trait states.
#' @param tol maximum tolerance. Default is 'c(1e-04, 1e-05, 1e-05)'.
#' @param maxiter max number of iterations. Default is 
#' '1000*round((1.25)^length(idparsopt))'.
#' @param optimmethod method used for optimization. Default is 'simplex'.
#' @param num_cycles number of cycles of the optimization (default is 1).
#' @param loglik_penalty the size of the penalty for all parameters; default is
#'  0 (no penalty)
#' @param is_complete_tree whether or not a tree with all its extinct species 
#' is provided
#' @param verbose sets verbose output; default is verbose when optimmethod is 
#' 'subplex'
#' @param num_threads number of threads. Set to -1 to use all available threads. 
#' Default is one thread.
#' @param atol absolute tolerance of integration
#' @param rtol relative tolerance of integration
#' @param method integration method used, available are: 
#' "odeint::runge_kutta_cash_karp54", "odeint::runge_kutta_fehlberg78", 
#' "odeint::runge_kutta_dopri5", "odeint::bulirsch_stoer" and 
#' "odeint::runge_kutta4". Default method is:"odeint::bulirsch_stoer".
#' @return Parameter estimated and maximum likelihood
#' @examples
#'# Example of how to set the arguments for a ML search.
#'library(secsse)
#'library(DDD)
#'set.seed(13)
#'# Check the vignette for a better working exercise. 
#'# lambdas for 0A and 1A and 2A are the same but need to be estimated 
#'# (CTD model, see Syst Biol paper)
#'# mus are fixed to zero, 
#'# the transition rates are constrained to be equal and fixed 0.01
#'phylotree <- ape::rcoal(31, tip.label = 1:31)
#'traits <-  sample(c(0,1,2), ape::Ntip(phylotree),replace=TRUE) #get some traits
#'num_concealed_states <- 3
#'idparslist <- cla_id_paramPos(traits,num_concealed_states)
#'idparslist$lambdas[1,] <- c(1,1,1,2,2,2,3,3,3)
#'idparslist[[2]][] <- 4
#'masterBlock <- matrix(5,ncol = 3,nrow = 3,byrow = TRUE) 
#'diag(masterBlock) <- NA
#'diff.conceal <- FALSE
#'idparslist[[3]] <- q_doubletrans(traits,masterBlock,diff.conceal)
#'startingpoint <- bd_ML(brts = ape::branching.times(phylotree))
#'intGuessLamba <- startingpoint$lambda0
#'intGuessMu <- startingpoint$mu0
#'idparsopt <- c(1,2,3)
#'initparsopt <- c(rep(intGuessLamba,3))
#'idparsfix <- c(0,4,5)
#'parsfix <- c(0,0,0.01)
#'tol <- c(1e-04, 1e-05, 1e-07)
#'maxiter <- 1000 * round((1.25) ^ length(idparsopt))
#'optimmethod <- 'simplex'
#'cond <- 'proper_cond'
#'root_state_weight <- 'proper_weights'
#'sampling_fraction <- c(1,1,1)
#' \dontrun{ 
#'model <- cla_secsse_ml(
#'  phylotree,
#'  traits,
#'  num_concealed_states,
#'  idparslist,
#'  idparsopt,
#'  initparsopt,
#'  idparsfix,
#'  parsfix,
#'  cond,
#'  root_state_weight,
#'  sampling_fraction,
#'  tol,
#'  maxiter,
#'  optimmethod,
#'  num_cycles = 1) }
#' # [1] -90.97626
#' @export
cla_secsse_timezones_ml <- function(phy,
                                    traits,
                                    num_concealed_states,
                                    idparslist,
                                    idparsopt,
                                    initparsopt,
                                    idparsfix,
                                    parsfix,
                                    critical_t = NULL,
                                    cond = "proper_cond", 
                                    root_state_weight = "proper_weights",
                                    sampling_fraction,
                                    tol = c(1e-04, 1e-05, 1e-07),
                                    maxiter = 1000 * round((1.25)^length(idparsopt)),
                                    optimmethod = "simplex",
                                    num_cycles = 1, 
                                    loglik_penalty = 0, 
                                    is_complete_tree = FALSE, 
                                    verbose = (optimmethod == "subplex"),
                                    num_threads = 1,
                                    atol = 1e-12,
                                    rtol = 1e-12,
                                    method = "odeint::bulirsch_stoer") {
  
  structure_func <- NULL
  if (is.matrix(traits)) {
    cat("you are setting a model where some species have more than one trait state \n")
  }
  
  if (length(initparsopt) != length(idparsopt)) {
    stop("initparsopt must be the same length as idparsopt. Number of parameters to optimize does not match the number of initial values for the search")
  }
  
  if (length(idparsfix) != length(parsfix)) {
    stop("idparsfix and parsfix must be the same length.Number of fixed elements does not match the fixed figures")
  }
  
  if (anyDuplicated(c(idparsopt, idparsfix)) != 0) {
    stop("at least one element was asked to be both fixed and estimated ")
  }
  
  if (identical(as.numeric(sort(c(idparsopt, idparsfix))), as.numeric(sort(unique(unlist(idparslist))))) == FALSE) {
    stop("All elements in idparslist must be included in either idparsopt or idparsfix ")
  }
  
  if (is.null(critical_t)) {
    if (anyDuplicated(c(unique(sort(as.vector(idparslist[[3]]))), idparsfix[which(parsfix == 0)])) != 0) {
      cat("Note: you set some transitions as impossible to happen.", "\n")
    }
  } else {
    for (j in 1:length(idparslist)) {
      if (anyDuplicated(c(unique(sort(as.vector(idparslist[[j]][[3]]))), idparsfix[which(parsfix == 0)])) != 0) {
        cat("Note: you set some transitions as impossible to happen.", "\n")
      } 
    }
  }
  
  for (t in 1:length(idparslist)) {
    if (is.matrix(idparslist[[t]][[1]])) {
      ## it is a tailor case otherwise
      idparslist[[t]][[1]] <- prepare_full_lambdas(traits, num_concealed_states, idparslist[[t]][[1]])
    }
  }
  
  if (!is.null(critical_t)) {
    len_t <- length(critical_t)
    len_par <- length(idparslist)
    if (1 + len_t != len_par) {
      stop("You have set transitions, but have not provided the right number of parameter sets")
    }
  } else {
    stop("you have not set time points at which new parameter regimes start")
  }
  
  see_ancestral_states <- FALSE
  
  # options(warn=-1)
  cat("Calculating the likelihood for the initial parameters ...", "\n")
  utils::flush.console()
  trparsopt <- initparsopt/(1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] <- 1
  trparsfix <- parsfix/(1 + parsfix)
  trparsfix[which(parsfix == Inf)] <- 1
  
  mus <- calc_mus(is_complete_tree,
                  idparslist[[1]],
                  idparsfix,
                  parsfix,
                  idparsopt,
                  initparsopt)
  
  
  optimpars <- c(tol, maxiter)
  
  setting_calculation <- build_initStates_time(phy, 
                                               traits,
                                               num_concealed_states,
                                               sampling_fraction,
                                               is_complete_tree, 
                                               mus)
  
  initloglik <- secsse_loglik_choosepar_timezones(trparsopt = trparsopt,
                                                  trparsfix = trparsfix,
                                                  idparsopt = idparsopt,
                                                  idparsfix = idparsfix, 
                                                  idparslist = idparslist,
                                                  critical_t = critical_t,
                                                  structure_func = structure_func,
                                                  phy = phy,
                                                  traits = traits,
                                                  num_concealed_states = num_concealed_states,
                                                  cond = cond,
                                                  root_state_weight = root_state_weight,
                                                  sampling_fraction = sampling_fraction, 
                                                  setting_calculation = setting_calculation,
                                                  see_ancestral_states = see_ancestral_states, 
                                                  loglik_penalty = loglik_penalty,
                                                  is_complete_tree = is_complete_tree, 
                                                  verbose = verbose,
                                                  num_threads = num_threads,
                                                  atol = atol,
                                                  rtol = rtol,
                                                  method = method)
  cat("The loglikelihood for the initial parameter values is", initloglik, "\n")
  if (initloglik == -Inf) {
    stop("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.")
  } else {
    cat("Optimizing the likelihood - this may take a while.", "\n")
    utils::flush.console()
    out <- DDD::optimizer(optimmethod = optimmethod,
                          optimpars = optimpars,
                          fun = secsse_loglik_choosepar_timezones,
                          trparsopt = trparsopt, 
                          idparsopt = idparsopt,
                          trparsfix = trparsfix,
                          idparsfix = idparsfix,
                          idparslist = idparslist,
                          critical_t = critical_t,
                          structure_func = structure_func,
                          phy = phy,
                          traits = traits,
                          num_concealed_states = num_concealed_states, 
                          cond = cond,
                          root_state_weight = root_state_weight,
                          sampling_fraction = sampling_fraction,
                          setting_calculation = setting_calculation,
                          see_ancestral_states = see_ancestral_states,
                          num_cycles = num_cycles,
                          loglik_penalty = loglik_penalty,
                          is_complete_tree = is_complete_tree, 
                          verbose = verbose,
                          num_threads = num_threads,
                          atol = atol,
                          rtol = rtol,
                          method = method)
    if (out$conv != 0) {
      stop("Optimization has not converged. Try again with different initial values.\n")
    } else {
      MLpars1 <- secsse_transform_parameters_timezone(as.numeric(unlist(out$par)),
                                                     trparsfix,
                                                     idparsopt,
                                                     idparsfix,
                                                     idparslist, 
                                                     structure_func)
      out2 <- list(MLpars = MLpars1,
                   ML = as.numeric(unlist(out$fvalues)),
                   conv = out$conv)
    }
  }
  return(out2)
}

secsse_loglik_choosepar_timezones <- function(trparsopt, 
                                              trparsfix, 
                                              idparsopt, 
                                              idparsfix, 
                                              idparslist, 
                                              critical_t,
                                              structure_func = structure_func,
                                              phy = phy, 
                                              traits = traits, 
                                              num_concealed_states = num_concealed_states, 
                                              cond = cond,
                                              root_state_weight = root_state_weight, 
                                              sampling_fraction = sampling_fraction, 
                                              setting_calculation = setting_calculation,
                                              see_ancestral_states = see_ancestral_states, 
                                              loglik_penalty = loglik_penalty,
                                              is_complete_tree = is_complete_tree, 
                                              verbose = verbose,
                                              num_threads = num_threads,
                                              atol = atol,
                                              rtol = rtol,
                                              method = method) {
  alltrpars <- c(trparsopt, trparsfix)
  if (max(alltrpars) > 1 | min(alltrpars) < 0) {
    loglik <- -Inf
  } else {
    pars1 <- secsse_transform_parameters_timezone(trparsopt, 
                                                  trparsfix, 
                                                  idparsopt, 
                                                  idparsfix, 
                                                  idparslist, 
                                                  structure_func)
    
    if (is.list(pars1[[1]][[1]])) {
      # is the cla_ used?
      loglik <- cla_secsse_loglik_timezones(parameter = pars1,
                                            phy = phy,
                                            traits = traits,
                                            num_concealed_states = num_concealed_states,
                                            critical_t = critical_t,
                                            cond = cond,
                                            root_state_weight = root_state_weight,
                                            sampling_fraction = sampling_fraction,
                                            setting_calculation = setting_calculation, 
                                            see_ancestral_states = see_ancestral_states,
                                            loglik_penalty = loglik_penalty,
                                            is_complete_tree = is_complete_tree,
                                            num_threads = num_threads,
                                            atol = atol,
                                            rtol = rtol,
                                            method = method)
    } else {
      loglik <- secsse_loglik_timezones(parameter = pars1,
                                        phy = phy,
                                        traits = traits,
                                        num_concealed_states = num_concealed_states,
                                        cond = cond,
                                        critical_t = critical_t,
                                        root_state_weight = root_state_weight,
                                        sampling_fraction = sampling_fraction,
                                        setting_calculation = setting_calculation,
                                        see_ancestral_states = see_ancestral_states,
                                        loglik_penalty = loglik_penalty,
                                        is_complete_tree = is_complete_tree,
                                        num_threads = num_threads,
                                        atol = atol,
                                        rtol = rtol,
                                        method = method)
    }
    if (is.nan(loglik) || is.na(loglik)) {
      cat("There are parameter values used which cause numerical problems.\n")
      loglik <- -Inf
    }
  }
  if (verbose) {
    cat(c(trparsopt/(1 - trparsopt), loglik), "\n")
  }
  return(loglik)
}

secsse_transform_parameters_timezone <- function(trparsopt,
                                                 trparsfix,
                                                 idparsopt,
                                                 idparsfix,
                                                 idparslist,
                                                 structure_func) {
  
  pars1_l <- list()
  
  if (is.null(structure_func) == FALSE) {
    idparsfuncdefpar <- structure_func[[1]]
    functions_defining_params <- structure_func[[2]]
    # idfactorsopt <- structure_func[[3]] <- idfactorsopt
    
    if (length(structure_func[[3]]) > 1) {
      idfactorsopt <- structure_func[[3]]
    } else {
      if (structure_func[[3]] == "noFactor") {
        idfactorsopt <- NULL
      } else {
        idfactorsopt <- structure_func[[3]]
      }
    }
    
    trparfuncdefpar <- transf_funcdefpar_timezones(idparsfuncdefpar = idparsfuncdefpar,
                                         functions_defining_params = functions_defining_params,
                                         idfactorsopt = idfactorsopt,
                                         trparsfix = trparsfix,
                                         trparsopt = trparsopt,
                                         idparsfix = idparsfix,
                                         idparsopt = idparsopt)
  }
  
  if (is.list(idparslist[[1]][[1]])) {
    # when the ml function is called from cla_secsse
    for (t in seq_along(idparslist)) {
      trpars1 <- idparslist[[t]]
      
      for (j in 1:nrow(trpars1[[3]])) {
        trpars1[[1]][[j]][, ] <- NA
      }
      
      for (j in 2:3) {
        trpars1[[j]][] <- NA
      }
      
      
      if (length(idparsfix) != 0) {
        
        for (i in seq_along(idparsfix)) {
          
          for (j in 1:nrow(trpars1[[3]])) {
            id <- which(idparslist[[t]][[1]][[j]] == idparsfix[i])
            trpars1[[1]][[j]][id] <- trparsfix[i]
          }
          for (j in 2:3) {
            id <- which(idparslist[[t]][[j]] == idparsfix[i])
            trpars1[[j]][id] <- trparsfix[i]
          }
        }
      }
      
      for (i in 1:length(idparsopt)) {
        for (j in 1:nrow(trpars1[[3]])) {
          id <- which(idparslist[[t]][[1]][[j]] == idparsopt[i])
          trpars1[[1]][[j]][id] <- trparsopt[i]
        }
        
        
        for (j in 2:3) {
          id <- which(idparslist[[t]][[j]] == idparsopt[i])
          trpars1[[j]][id] <- trparsopt[i]
        }
      }
      
      ## structure_func part
      
      if (is.null(structure_func) == FALSE) {
        for (i in 1:length(idparsfuncdefpar)) {
          for (j in 1:nrow(trpars1[[3]])) {
            id <- which(idparslist[[t]][[1]][[j]] == idparsfuncdefpar[i])
            trpars1[[1]][[j]][id] <- trparfuncdefpar[i]
          }
          
          for (j in 2:3) {
            id <- which(idparslist[[t]][[1]][[j]] == idparsfuncdefpar[i])
            trpars1[[j]][id] <- trparfuncdefpar[i]
          }
        }
      }
      
      pre_pars1 <- list()
      
      for (j in 1:nrow(trpars1[[3]])) {
        pre_pars1[[j]] <- trpars1[[1]][[j]][, ]/(1 - trpars1[[1]][[j]][, ])
      }
      
      pars1 <- list()
      pars1[[1]] <- pre_pars1
      for (j in 2:3) {
        pars1[[j]] <- trpars1[[j]]/(1 - trpars1[[j]])
      }
      
      pars1_l[[t]] <- pars1
    }
  } else {
    #### when non-cla option is called
    for (t in 1:length(idparslist)) {
    
      trpars1 <- idparslist[[t]]
      for (j in 1:3) {
        trpars1[[j]][] <- NA
      }
      if (length(idparsfix) != 0) {
        for (i in seq_along(idparsfix)) {
          for (j in 1:3) {
            id <- which(idparslist[[t]][[j]] == idparsfix[i])
            trpars1[[j]][id] <- trparsfix[i]
          }
        }
      }
      for (i in seq_along(idparsopt)) {
        for (j in 1:3) {
          id <- which(idparslist[[t]][[j]] == idparsopt[i])
          trpars1[[j]][id] <- trparsopt[i]
        }
      }
      
      ## if structure_func part
      
      if (is.null(structure_func) == FALSE) {
        for (i in seq_along(idparsfuncdefpar)) {
          for (j in 1:3) {
            id <- which(idparslist[[t]][[j]] == idparsfuncdefpar[i])
            trpars1[[j]][id] <- trparfuncdefpar[i]
          }
        }
      }
      
      pars1_l[[t]] <- list()
      for (j in 1:3) {
        pars1_l[[t]][[j]] <- trpars1[[j]]/(1 - trpars1[[j]])
      }
    }
  }
  return(pars1_l)
}

transf_funcdefpar_timezones <- function(idparsfuncdefpar,
                                        functions_defining_params,
                                        idfactorsopt,
                                        trparsfix,
                                        trparsopt,
                                        idparsfix,
                                        idparsopt) {
  trparfuncdefpar <- NULL
  ids_all <- c(idparsfix, idparsopt)
  
  values_all <- c(trparsfix/(1 - trparsfix), trparsopt/(1 - trparsopt))
  a_new_envir <- new.env()
  x <- as.list(values_all)  ## To declare all the ids as variables
  
  if (is.null(idfactorsopt)) {
    names(x) <- paste0("par_", ids_all)
  } else {
    names(x) <- c(paste0("par_", ids_all), paste0("factor_", idfactorsopt))
  }
  list2env(x, envir = a_new_envir)
  
  for (jj in seq_along(functions_defining_params)) {
    myfunc <- functions_defining_params[[jj]]
    environment(myfunc) <- a_new_envir
    value_func_defining_parm <- local(myfunc(), envir = a_new_envir)
    
    ## Now, declare the variable that is just calculated, so it is available 
    ## for the next calculation if needed
    y <- as.list(value_func_defining_parm)
    names(y) <- paste0("par_", idparsfuncdefpar[jj])
    list2env(y, envir = a_new_envir)
    
    if (is.numeric(value_func_defining_parm) == FALSE) {
      stop("Something went wrong with the calculation of parameters in 'functions_param_struct'")
    }
    trparfuncdefpar <- c(trparfuncdefpar, value_func_defining_parm)
  }
  trparfuncdefpar <- trparfuncdefpar/(1 + trparfuncdefpar)
  rm(a_new_envir)
  return(trparfuncdefpar)
}
