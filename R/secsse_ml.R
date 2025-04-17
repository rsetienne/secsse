#' @keywords internal
master_ml <- function(phy,
                      traits,
                      num_concealed_states,
                      idparslist,
                      idparsopt,
                      initparsopt,
                      idparsfix,
                      parsfix,
                      idfactorsopt = NULL,
                      initfactors = NULL,
                      idparsfuncdefpar = NULL,
                      functions_defining_params = NULL,
                      cond = "proper_cond",
                      root_state_weight = "proper_weights",
                      sampling_fraction,
                      tol = c(1e-04, 1e-05, 1e-07),
                      maxiter = 1000 * round((1.25) ^ length(idparsopt)),
                      optimmethod = "simplex",
                      num_cycles = 1,
                      loglik_penalty = 0,
                      is_complete_tree = FALSE,
                      take_into_account_root_edge = take_into_account_root_edge,
                      verbose = FALSE,
                      num_threads = 1,
                      atol = 1e-8,
                      rtol = 1e-7,
                      method = "odeint::bulirsch_stoer",
                      use_normalization = TRUE) {
  
  structure_func <- NULL
  if (!is.null(functions_defining_params)) {
    structure_func <- set_and_check_structure_func(idparsfuncdefpar,
                                                   functions_defining_params,
                                                   idparslist,
                                                   idparsopt,
                                                   idfactorsopt,
                                                   idparsfix,
                                                   initfactors)
  } else {
    if (identical(as.numeric(sort(c(idparsopt, idparsfix))),
                  as.numeric(sort(unique(unlist(idparslist))))) == FALSE) {
      stop("All elements in idparslist must be included in either
             idparsopt or idparsfix ")
    }
  }
  
  check_ml_conditions(traits,
                      idparslist,
                      initparsopt,
                      idparsopt,
                      idparsfix,
                      parsfix)
  
  if (is.matrix(idparslist[[1]])) {
    ## it is a tailor case otherwise
    idparslist[[1]] <- prepare_full_lambdas(traits,
                                            num_concealed_states,
                                            idparslist[[1]])
  }
  
  see_ancestral_states <- FALSE
  if (!is.null(structure_func)) {
    initparsopt <- c(initparsopt, initfactors)
  }
  
  trparsopt <- initparsopt / (1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] <- 1
  trparsfix <- parsfix / (1 + parsfix)
  trparsfix[which(parsfix == Inf)] <- 1
  
  mus <- calc_mus(is_complete_tree,
                  idparslist,
                  idparsfix,
                  parsfix,
                  idparsopt,
                  initparsopt)
  
  optimpars <- c(tol, maxiter, verbose)
  
  num_modeled_traits <- length(idparslist[[1]]) / num_concealed_states
  
  if (!is.list(traits)) {
    
    setting_calculation <- build_initStates_time(phy,
                                                 traits,
                                                 num_concealed_states,
                                                 sampling_fraction,
                                                 is_complete_tree,
                                                 mus,
                                                 num_modeled_traits,
                                                 traitStates = 
                                                   get_trait_states(idparslist,
                                                                    num_concealed_states, FALSE))
  } else {
    setting_calculation <- list()
    for (i in 1:length(phy)) {
      
      input_phy <- phy[[i]]
      input_traits <- traits[[i]]
      
      if (is.list(sampling_fraction)) { # weirdly, ifelse(is.list) doesn't work
        input_sampling_fraction <- sampling_fraction[[i]]
      } else {
        input_sampling_fraction <- sampling_fraction
      }
      
      if (length(input_phy$tip.label) == 1) {
        fake_phy <- ape::rphylo(n = 2, birth = 1, death = 0)
        fake_phy$edge.length[1:2] <- input_phy$edge.length[1]
        input_phy <- fake_phy
        input_traits <- c(input_traits, input_traits)
      }
      
      setting_calculation[[i]] <- build_initStates_time(phy = input_phy,
                                                   traits = input_traits,
                                                   num_concealed_states =
                                                     num_concealed_states,
                                                   sampling_fraction =
                                                     input_sampling_fraction,
                                                   is_complete_tree =
                                                     is_complete_tree,
                                                   mus = mus,
                                                   num_unique_traits = 
                                                     num_modeled_traits,
                                                   first_time = FALSE,
                                                   traitStates = 
                                                     get_trait_states(idparslist,
                                                                      num_concealed_states, FALSE))
    }
  }
  
  ll_verbose <- ifelse(optimmethod == "subplex",
                       verbose,
                       FALSE)
  initloglik <- secsse_loglik_choosepar(trparsopt = trparsopt,
                                        trparsfix = trparsfix,
                                        idparsopt = idparsopt,
                                        idparsfix = idparsfix,
                                        idparslist = idparslist,
                                        structure_func = structure_func,
                                        phy = phy,
                                        traits = traits,
                                        num_concealed_states =
                                          num_concealed_states,
                                        cond = cond,
                                        root_state_weight = root_state_weight,
                                        sampling_fraction = sampling_fraction,
                                        setting_calculation =
                                          setting_calculation,
                                        see_ancestral_states =
                                          see_ancestral_states,
                                        loglik_penalty = loglik_penalty,
                                        is_complete_tree = is_complete_tree,
                                        take_into_account_root_edge =
                                          take_into_account_root_edge,
                                        num_threads = num_threads,
                                        atol = atol,
                                        rtol = rtol,
                                        method = method,
                                        display_warning = FALSE,
                                        verbose = ll_verbose,
                                        use_normalization = use_normalization)
  # Function here
  if (verbose) print_init_ll(initloglik = initloglik)

  if (initloglik == -Inf) {
    stop("The initial parameter values have a likelihood that is 
             equal to 0 or below machine precision. 
             Try again with different initial values.")
  } else {
    out <- DDD::optimizer(optimmethod = optimmethod,
                          optimpars = optimpars,
                          fun = secsse_loglik_choosepar,
                          trparsopt = trparsopt,
                          num_cycles = num_cycles,
                          idparsopt = idparsopt,
                          trparsfix = trparsfix,
                          idparsfix = idparsfix,
                          idparslist = idparslist,
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
                          take_into_account_root_edge = 
                            take_into_account_root_edge,
                          num_threads = num_threads,
                          atol = atol,
                          rtol = rtol,
                          method = method,
                          display_warning = FALSE,
                          verbose = ll_verbose,
                          use_normalization = use_normalization)
    if (out$conv != 0) {
      stop("Optimization has not converged. 
                 Try again with different initial values.")
    } else {
      ml_pars1 <- secsse_transform_parameters(as.numeric(unlist(out$par)),
                                              trparsfix,
                                              idparsopt,
                                              idparsfix,
                                              idparslist,
                                              structure_func)
      out2 <- list(MLpars = ml_pars1,
                   ML = as.numeric(unlist(out$fvalues)),
                   conv = out$conv)
    }
  }
  return(out2)
}

#' Maximum likehood estimation for (SecSSE)
#' 
#' Maximum likehood estimation under Several examined and concealed
#' States-dependent Speciation and Extinction (SecSSE)
#' @inheritParams default_params_doc
#' 
#' @return Parameter estimated and maximum likelihood
#' @examples
#'# Example of how to set the arguments for a ML search.
#'library(secsse)
#'library(DDD)
#'set.seed(13)
#'# lambdas for 0A and 1A and 2A are the same but need to be estimated
#'# mus are fixed to
#'# the transition rates are constrained to be equal and fixed 0.01
#'phylotree <- ape::rcoal(31, tip.label = 1:31)
#'traits <-  sample(c(0,1,2), ape::Ntip(phylotree),replace=TRUE)#get some traits
#'num_concealed_states<-3
#'idparslist <- id_paramPos(traits, num_concealed_states)
#'idparslist[[1]][c(1,4,7)] <- 1
#'idparslist[[1]][c(2,5,8)] <- 2
#'idparslist[[1]][c(3,6,9)] <- 3
#'idparslist[[2]][]<-4
#'masterBlock <- matrix(5,ncol = 3,nrow = 3,byrow = TRUE)
#'diag(masterBlock) <- NA
#'diff.conceal <- FALSE
#'idparslist[[3]] <- q_doubletrans(traits,masterBlock,diff.conceal)
#'startingpoint <- DDD::bd_ML(brts = ape::branching.times(phylotree))
#'intGuessLamba <- startingpoint$lambda0
#'intGuessMu <- startingpoint$mu0
#'idparsopt <- c(1,2,3,5)
#'initparsopt <- c(rep(intGuessLamba,3),rep((intGuessLamba/5),1))
#'idparsfix <- c(0,4)
#'parsfix <- c(0,0)
#'tol <- c(1e-02, 1e-03, 1e-04)
#'maxiter <- 1000 * round((1.25)^length(idparsopt))
#'optimmethod <- 'simplex'
#'cond <- 'proper_cond'
#'root_state_weight <- 'proper_weights'
#'sampling_fraction <- c(1,1,1)
#'model<-secsse_ml(
#'phylotree,
#'traits,
#'num_concealed_states,
#'idparslist,
#'idparsopt,
#'initparsopt,
#'idparsfix,
#'parsfix,
#'cond,
#'root_state_weight,
#'sampling_fraction,
#'tol,
#'maxiter,
#'optimmethod,
#'num_cycles = 1,
#'verbose = FALSE)
#'# model$ML
#'# [1] -16.04127
#' @export
secsse_ml <- function(phy,
                      traits,
                      num_concealed_states,
                      idparslist,
                      idparsopt,
                      initparsopt,
                      idparsfix,
                      parsfix,
                      cond = "proper_cond",
                      root_state_weight = "proper_weights",
                      sampling_fraction,
                      tol = c(1e-04, 1e-05, 1e-07),
                      maxiter = 1000 * round((1.25)^length(idparsopt)),
                      optimmethod = "simplex",
                      num_cycles = 1,
                      loglik_penalty = 0,
                      is_complete_tree = FALSE,
                      take_into_account_root_edge = FALSE,
                      verbose = FALSE,
                      num_threads = 1,
                      atol = 1e-8,
                      rtol = 1e-7,
                      method = "odeint::bulirsch_stoer",
                      use_normalization = TRUE) {
  master_ml(phy = phy,
            traits = traits,
            num_concealed_states = num_concealed_states,
            idparslist = idparslist,
            idparsopt = idparsopt,
            initparsopt = initparsopt,
            idparsfix = idparsfix,
            parsfix = parsfix,
            initfactors = NULL,
            idparsfuncdefpar = NULL,
            cond = cond,
            root_state_weight = root_state_weight,
            sampling_fraction = sampling_fraction,
            tol = tol,
            maxiter = maxiter,
            optimmethod = optimmethod,
            num_cycles = num_cycles,
            loglik_penalty = loglik_penalty,
            is_complete_tree = is_complete_tree,
            take_into_account_root_edge = take_into_account_root_edge,
            verbose = verbose,
            num_threads = num_threads,
            atol = atol,
            rtol = rtol,
            method = method,
            use_normalization = use_normalization)
}

#' @keywords internal
secsse_loglik_choosepar <- function(trparsopt,
                                    trparsfix,
                                    idparsopt,
                                    idparsfix,
                                    idparslist,
                                    structure_func,
                                    phy,
                                    traits,
                                    num_concealed_states,
                                    cond = cond,
                                    root_state_weight,
                                    sampling_fraction,
                                    setting_calculation,
                                    see_ancestral_states,
                                    loglik_penalty,
                                    is_complete_tree,
                                    take_into_account_root_edge,
                                    num_threads,
                                    atol,
                                    rtol,
                                    method,
                                    #structure_func = structure_func,
                                    #phy = phy,
                                    #traits = traits,
                                    #num_concealed_states = num_concealed_states,
                                    #cond = cond,
                                    #root_state_weight = root_state_weight,
                                    #sampling_fraction = sampling_fraction,
                                    #setting_calculation = setting_calculation,
                                    #see_ancestral_states = see_ancestral_states,
                                    #loglik_penalty = loglik_penalty,
                                    #is_complete_tree = is_complete_tree,
                                    #take_into_account_root_edge = take_into_account_root_edge,
                                    #num_threads = num_threads,
                                    #atol = atol,
                                    #rtol = rtol,
                                    #method = method,
                                    display_warning,
                                    verbose,
                                    use_normalization) {
  alltrpars <- c(trparsopt, trparsfix)
  if (max(alltrpars) > 1 || min(alltrpars) < 0) {
    loglik <- -Inf
  } else {
    pars1 <- secsse_transform_parameters(trparsopt, trparsfix,
                                         idparsopt, idparsfix,
                                         idparslist, structure_func)
    
    loglik <- master_loglik(parameter = pars1,
                            phy = phy,
                            traits = traits,
                            num_concealed_states =
                              num_concealed_states,
                            cond = cond,
                            root_state_weight =
                              root_state_weight,
                            sampling_fraction =
                              sampling_fraction,
                            setting_calculation =
                              setting_calculation,
                            see_ancestral_states =
                              see_ancestral_states,
                            loglik_penalty = loglik_penalty,
                            is_complete_tree =
                              is_complete_tree,
                            take_into_account_root_edge =
                              take_into_account_root_edge,
                            num_threads = num_threads,
                            method = method,
                            atol = atol,
                            rtol = rtol,
                            display_warning = display_warning,
                            use_normalization = use_normalization)
    
    if (is.nan(loglik) || is.na(loglik)) {
      warning("There are parameter values used which cause
                numerical problems.")
      loglik <- -Inf
    }
  }
  if (verbose) {
    out_print <- c(trparsopt / (1 - trparsopt), loglik)
    message(paste(out_print, collapse = " "))
  }
  return(loglik)
}

#' Maximum likehood estimation for (SecSSE)
#' 
#' Maximum likehood estimation under Several examined and concealed
#' States-dependent Speciation and Extinction (SecSSE) with cladogenetic option
#'
#' @inheritParams default_params_doc
#' 
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
#'#get some traits
#'traits <-  sample(c(0,1,2), ape::Ntip(phylotree), replace = TRUE)
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
#'  num_cycles = 1,
#'  num_threads = 1,
#'  verbose = FALSE)
#' # [1] -90.97626
#' @export
cla_secsse_ml <- function(phy,
                          traits,
                          num_concealed_states,
                          idparslist,
                          idparsopt,
                          initparsopt,
                          idparsfix,
                          parsfix,
                          cond = "proper_cond",
                          root_state_weight = "proper_weights",
                          sampling_fraction,
                          tol = c(1e-04, 1e-05, 1e-07),
                          maxiter = 1000 * round((1.25)^length(idparsopt)),
                          optimmethod = "simplex",
                          num_cycles = 1,
                          loglik_penalty = 0,
                          is_complete_tree = FALSE,
                          take_into_account_root_edge = FALSE,
                          verbose = FALSE,
                          num_threads = 1,
                          atol = 1e-8,
                          rtol = 1e-7,
                          method = "odeint::bulirsch_stoer",
                          use_normalization = TRUE) {
  master_ml(phy = phy,
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
            num_cycles = num_cycles,
            loglik_penalty = loglik_penalty,
            is_complete_tree = is_complete_tree,
            take_into_account_root_edge = take_into_account_root_edge,
            verbose = verbose,
            num_threads = num_threads,
            atol = atol,
            rtol = rtol,
            method = method,
            use_normalization = use_normalization)
}
