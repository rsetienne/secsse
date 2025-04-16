#' @keywords internal
master_loglik <- function(parameter,
                          phy,
                          traits,
                          num_concealed_states,
                          cond = "proper_cond",
                          root_state_weight = "proper_weights",
                          sampling_fraction,
                          setting_calculation = NULL,
                          see_ancestral_states = FALSE,
                          loglik_penalty = 0,
                          is_complete_tree = FALSE,
                          num_threads = 1,
                          atol = 1e-8,
                          rtol = 1e-7,
                          method = "odeint::bulirsch_stoer",
                          take_into_account_root_edge = FALSE,
                          display_warning = TRUE,
                          use_normalization = TRUE) {
  
  if (is.list(phy)) {
    if (!inherits(phy, "phylo")) {
      if (!inherits(phy, "multiPhylo")) {
        stop("when providing multiple phylogenies, make sure to use the multiPhylo class")
      }
    }
  }
  
  
  if (inherits(phy, "multiPhylo")) {
    if (!is.list(traits)) {
      stop("traits needs to be supplied as a list now that there are multiple phylogenies")
    }
    return(multi_loglik(parameter = parameter,
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
                        take_into_account_root_edge = take_into_account_root_edge,
                        num_threads = num_threads,
                        atol = atol,
                        rtol = rtol,
                        method = method,
                        display_warning = display_warning,
                        use_normalization = use_normalization))
  }
  
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  q_matrix <- parameter[[3]]
  
  using_cla <- is.list(lambdas)
  
  num_modeled_traits <- ncol(q_matrix) / floor(num_concealed_states)
  
  traitStates = get_trait_states(parameter,
                                 num_concealed_states, display_warning)
  
  if (is.null(setting_calculation)) {
    check_input(traits,
                phy,
                sampling_fraction,
                root_state_weight,
                is_complete_tree)
    setting_calculation <- build_initStates_time(phy,
                                                 traits,
                                                 num_concealed_states,
                                                 sampling_fraction,
                                                 is_complete_tree,
                                                 mus,
                                                 num_modeled_traits,
                                                 traitStates = traitStates)
  } 
  
  states <- setting_calculation$states
  forTime <- setting_calculation$forTime
  ances <- setting_calculation$ances
  
  d <- ncol(states) / 2
  
  # with a complete tree, we need to re-calculate the states every time we
  # run, because they are dependent on mu.
  if (is_complete_tree) {
    states <- build_states(phy = phy,
                           traits = traits,
                           num_concealed_states = num_concealed_states,
                           sampling_fraction = sampling_fraction,
                           is_complete_tree = is_complete_tree,
                           mus = mus,
                           num_unique_traits = num_modeled_traits,
                           traitStates = traitStates)
  }
  
  RcppParallel::setThreadOptions(numThreads = num_threads)
  calcul <- calc_ll_cpp(rhs = if (using_cla) "ode_cla" else "ode_standard",
                        ances = ances,
                        states = states,
                        forTime = forTime,
                        lambdas = lambdas,
                        mus = mus,
                        Q = q_matrix,
                        method = method,
                        atol = atol,
                        rtol = rtol,
                        is_complete_tree = is_complete_tree,
                        see_states = see_ancestral_states,
                        use_normalization = use_normalization)
  loglik <- calcul$loglik
  nodeM <- calcul$node_M
  mergeBranch <- calcul$merge_branch
  
  if (length(nodeM) > 2 * d) nodeM <- nodeM[1:(2 * d)]
  
  if (!is.null(phy$root.edge) && take_into_account_root_edge == TRUE ) {
    if (phy$root.edge > 0) {
      calcul2 <- calc_ll_single_branch_cpp(rhs = 
                                             if (using_cla) "ode_cla" else "ode_standard",
                                           states = c(nodeM[1:d], mergeBranch),
                                           forTime = c(0, phy$root.edge),
                                           lambdas = lambdas,
                                           mus = mus,
                                           Q = q_matrix,
                                           method = method,
                                           atol = atol,
                                           rtol = rtol,
                                           see_states = see_ancestral_states,
                                           use_normalization = use_normalization)
      loglik <- loglik + calcul2$loglik
      nodeM <- calcul2$states
      
      mergeBranch <- calcul2$merge_branch
    }
  }
  
  ## At the root
  weight_states <- get_weight_states(root_state_weight,
                                     num_concealed_states,
                                     mergeBranch,
                                     lambdas,
                                     nodeM,
                                     d,
                                     is_cla = using_cla)
  
  if (is_complete_tree) nodeM <- update_complete_tree(phy,
                                                      lambdas,
                                                      mus,
                                                      q_matrix,
                                                      method,
                                                      atol,
                                                      rtol,
                                                      length(mergeBranch),
                                                      use_normalization)
  
  mergeBranch2 <- condition(cond,
                            mergeBranch,
                            weight_states,
                            lambdas,
                            nodeM,
                            is_root_edge = take_into_account_root_edge)
  
  wholeLike <- sum((mergeBranch2) * (weight_states))
  
  LL <- log(wholeLike) +
    loglik -
    penalty(pars = parameter, loglik_penalty = loglik_penalty)
  
  if (see_ancestral_states == TRUE) {
    states <- calcul$states
    num_tips <- ape::Ntip(phy)
    ancestral_states <- states[(num_tips + 1):(nrow(states)), ]
    ancestral_states <-
      ancestral_states[, -1 * (1:(ncol(ancestral_states) / 2))]
    rownames(ancestral_states) <- ances
    return(list(ancestral_states = ancestral_states, LL = LL, states = states))
  } else {
    return(LL)
  }
}

#' @title Likelihood for SecSSE model
#' Loglikelihood calculation for the SecSSE model given a set of parameters and
#' data
#' 
#' @inheritParams default_params_doc
#' @return The loglikelihood of the data given the parameter.
#' @examples
#' rm(list = ls(all = TRUE))
#' library(secsse)
#' set.seed(13)
#' phylotree <- ape::rcoal(31, tip.label = 1:31)
#' traits <- sample(c(0,1,2),ape::Ntip(phylotree),replace = TRUE)
#' num_concealed_states <- 2
#' cond <- "proper_cond"
#' root_state_weight <- "proper_weights"
#' sampling_fraction <- c(1,1,1)
#' drill <- id_paramPos(traits,num_concealed_states)
#' drill[[1]][] <- c(0.12,0.01,0.2,0.21,0.31,0.23)
#' drill[[2]][] <- 0
#' drill[[3]][,] <- 0.1
#' diag(drill[[3]]) <- NA
#' secsse_loglik(parameter = drill,
#' phylotree,
#' traits,
#' num_concealed_states,
#' cond,
#' root_state_weight,
#' sampling_fraction,
#' see_ancestral_states = FALSE)
#'
#' #[1] -113.1018
#' @export
secsse_loglik <- function(parameter,
                          phy,
                          traits,
                          num_concealed_states,
                          cond = "proper_cond",
                          root_state_weight = "proper_weights",
                          sampling_fraction,
                          setting_calculation = NULL,
                          see_ancestral_states = FALSE,
                          loglik_penalty = 0,
                          is_complete_tree = FALSE,
                          take_into_account_root_edge = FALSE,
                          num_threads = 1,
                          atol = 1e-8,
                          rtol = 1e-7,
                          method = "odeint::bulirsch_stoer",
                          display_warning = TRUE,
                          use_normalization = TRUE) {
  master_loglik(parameter = parameter,
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
                take_into_account_root_edge = take_into_account_root_edge,
                num_threads = num_threads,
                atol = atol,
                rtol = rtol,
                method = method,
                display_warning = display_warning,
                use_normalization = use_normalization)
}


#' @keywords internal
multi_loglik <- function(parameter,
                         phy,
                         traits,
                         num_concealed_states,
                         cond = "proper_cond",
                         root_state_weight = "proper_weights",
                         sampling_fraction,
                         setting_calculation = NULL,
                         see_ancestral_states = FALSE,
                         loglik_penalty = 0,
                         is_complete_tree = FALSE,
                         take_into_account_root_edge = FALSE,
                         num_threads = 1,
                         atol = 1e-8,
                         rtol = 1e-7,
                         method = "odeint::bulirsch_stoer",
                         display_warning = FALSE,
                         use_normalization = TRUE) {
  
  res <- list()
  for (i in 1:length(phy)) {
    if (is.list(sampling_fraction)) {
      focal_sampling_fraction <- sampling_fraction[[i]]
    } else {
      focal_sampling_fraction <- sampling_fraction
    }
    
    if (is.list(root_state_weight)) {
      if (sum(is.na(root_state_weight[[i]])) || 
          length(root_state_weight[[i]]) < 1) {
        focal_root_state_weight <- "proper_weights"
      } else {
        focal_root_state_weight <- root_state_weight[[i]]
      }
    } else {
      focal_root_state_weight <- root_state_weight
    }
    
    focal_setting_calculation <- NULL
    if (is.list(setting_calculation)) {
      focal_setting_calculation <- setting_calculation[[i]]
    }
    
    if (length(phy[[i]]$tip.label) == 1) {
      res[[i]] <- secsse::secsse_single_branch_loglik(parameter = parameter,
                                                      phy = phy[[i]],
                                                      traits = traits[[i]],
                                                      num_concealed_states =
                                                        num_concealed_states,
                                                      cond = cond,
                                                      root_state_weight = 
                                                        focal_root_state_weight,
                                                      sampling_fraction = 
                                                        focal_sampling_fraction,
                                                      setting_calculation = 
                                                        focal_setting_calculation,
                                                      see_ancestral_states = FALSE,
                                                      loglik_penalty = loglik_penalty,
                                                      is_complete_tree = 
                                                        is_complete_tree,
                                                      take_into_account_root_edge = 
                                                        take_into_account_root_edge,
                                                      num_threads = num_threads,
                                                      atol = atol,
                                                      rtol = rtol,
                                                      method = method,
                                                      display_warning = display_warning,
                                                      use_normalization = use_normalization)$loglik
    } else {
      res[[i]] <- secsse_loglik(parameter = parameter,
                                phy = phy[[i]],
                                traits = traits[[i]],
                                num_concealed_states = num_concealed_states,
                                cond = cond,
                                root_state_weight = focal_root_state_weight,
                                sampling_fraction = focal_sampling_fraction,
                                setting_calculation = focal_setting_calculation,
                                see_ancestral_states = FALSE,
                                loglik_penalty = loglik_penalty,
                                is_complete_tree = is_complete_tree,
                                take_into_account_root_edge = 
                                  take_into_account_root_edge,
                                num_threads = num_threads,
                                atol = atol,
                                rtol = rtol,
                                method = method,
                                display_warning = display_warning,
                                use_normalization = use_normalization) 
    }
  }
  
  ll <- do.call(sum, res)
  
  return(ll) 
}

#' @title Likelihood for SecSSE model, using Rcpp
#' Loglikelihood calculation for the cla_SecSSE model given a set of parameters
#' and data using Rcpp
#' 
#' @inheritParams default_params_doc
#' 
#' @return The loglikelihood of the data given the parameters
#' @examples
#'rm(list=ls(all=TRUE))
#'library(secsse)
#'set.seed(13)
#'phylotree <- ape::rcoal(12, tip.label = 1:12)
#'traits <- sample(c(0,1,2),ape::Ntip(phylotree),replace=TRUE)
#'num_concealed_states <- 3
#'sampling_fraction <- c(1,1,1)
#'phy <- phylotree
#'# the idparlist for a ETD model (dual state inheritance model of evolution)
#'# would be set like this:
#'idparlist <- cla_id_paramPos(traits,num_concealed_states)
#'lambd_and_modeSpe <- idparlist$lambdas
#'lambd_and_modeSpe[1,] <- c(1,1,1,2,2,2,3,3,3)
#'idparlist[[1]] <- lambd_and_modeSpe
#'idparlist[[2]][] <- 0
#'masterBlock <- matrix(4,ncol=3,nrow=3,byrow=TRUE)
#'diag(masterBlock) <- NA
#'idparlist [[3]] <- q_doubletrans(traits,masterBlock,diff.conceal = FALSE)
#'# Now, internally, clasecsse sorts the lambda matrices, so they look like:
#'prepare_full_lambdas(traits,num_concealed_states,idparlist[[1]])
#'# which is a list with 9 matrices, corresponding to the 9 states
#'# (0A,1A,2A,0B,etc)
#'# if we want to calculate a single likelihood:
#'parameter <- idparlist
#'lambda_and_modeSpe <- parameter$lambdas
#'lambda_and_modeSpe[1,] <- c(0.2,0.2,0.2,0.4,0.4,0.4,0.01,0.01,0.01)
#'parameter[[1]] <- prepare_full_lambdas(traits,num_concealed_states,
#'lambda_and_modeSpe)
#'parameter[[2]] <- rep(0,9)
#'masterBlock <- matrix(0.07, ncol=3, nrow=3, byrow=TRUE)
#'diag(masterBlock) <- NA
#'parameter [[3]] <- q_doubletrans(traits,masterBlock,diff.conceal = FALSE)
#'cla_secsse_loglik(parameter, phy, traits, num_concealed_states,
#'                  cond = 'maddison_cond',
#'                  root_state_weight = 'maddison_weights', sampling_fraction,
#'                  setting_calculation = NULL,
#'                  see_ancestral_states = FALSE,
#'                  loglik_penalty = 0)
#'# LL = -42.18407
#' @export
cla_secsse_loglik <- function(parameter,
                              phy,
                              traits,
                              num_concealed_states,
                              cond = "proper_cond",
                              root_state_weight = "proper_weights",
                              sampling_fraction,
                              setting_calculation = NULL,
                              see_ancestral_states = FALSE,
                              loglik_penalty = 0,
                              is_complete_tree = FALSE,
                              take_into_account_root_edge = FALSE,
                              num_threads = 1,
                              method = "odeint::bulirsch_stoer",
                              atol = 1e-8,
                              rtol = 1e-7,
                              display_warning = TRUE,
                              use_normalization = TRUE) {
  master_loglik(parameter = parameter,
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
                take_into_account_root_edge = take_into_account_root_edge,
                num_threads = num_threads,
                atol = atol,
                rtol = rtol,
                method = method,
                display_warning = display_warning,
                use_normalization = use_normalization)
}
