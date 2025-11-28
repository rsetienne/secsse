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
                          method = "odeint::runge_kutta_cash_karp54",
                          take_into_account_root_edge = FALSE,
                          display_warning = TRUE,
                          use_normalization = TRUE,
                          return_root_state = FALSE) {
  
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
                        use_normalization = use_normalization,
                        return_root_state = return_root_state))
  }
  
  if (length(phy$tip.label) == 1) {
    return(secsse::secsse_single_branch_loglik(
                       parameter = parameter,
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
                       see_ancestral_states = see_ancestral_states,
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
                       use_normalization = use_normalization,
                       return_root_state = return_root_state))
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
  
  d <- ncol(states) / 3
  
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
  return_states = see_ancestral_states
  if (return_root_state) return_states = TRUE
  
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
                        see_states = return_states,
                        use_normalization = use_normalization)
  loglik <- calcul$loglik
  nodeM <- calcul$node_M
  mergeBranch <- calcul$merge_branch
  
  E <- nodeM[1:d]
  S <- nodeM[(2 * d + 1):(3 * d)]
  
  if (using_cla && !is_complete_tree) {
    
    if (any(is.na(S))) {
      S <- 1 - E
    } 
  } else {
    # currently, S is not implemented in complete_tree LL
    S <- 1 - E
  }
  
  if (!is.null(phy$root.edge) && take_into_account_root_edge == TRUE ) {
    if (phy$root.edge > 0) {
      calcul2 <- calc_ll_single_branch_cpp(rhs = 
                                             if (using_cla) "ode_cla" else "ode_standard",
                                           states = c(E, mergeBranch, S),
                                           forTime = c(0, phy$root.edge),
                                           lambdas = lambdas,
                                           mus = mus,
                                           Q = q_matrix,
                                           method = method,
                                           atol = atol,
                                           rtol = rtol,
                                           see_states = return_states,
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
  
  if (is_complete_tree) {
    nodeM <- update_complete_tree(phy,
                                  lambdas,
                                  mus,
                                  q_matrix,
                                  method,
                                  atol,
                                  rtol,
                                  length(mergeBranch),
                                  use_normalization)
    # TODO: fix this cheating way of implementing survival for CT
    E <- nodeM[1:d]
    S <- 1 - E
  }
  
  
  mergeBranch2 <- condition(cond,
                            mergeBranch,
                            weight_states,
                            lambdas,
                            is_root_edge = take_into_account_root_edge,
                            S)
  
  wholeLike <- sum( (mergeBranch2) * (weight_states) )
  
  LL <- log(wholeLike) +
    loglik -
    penalty(pars = parameter, loglik_penalty = loglik_penalty)
  
  # reset number of threads:
  RcppParallel::setThreadOptions(numThreads = 1)

  result <- list()
  result$LL <- LL
  if (see_ancestral_states == TRUE) {
    states <- calcul$states
    num_tips <- ape::Ntip(phy)
    ancestral_states <- states[(num_tips + 1):(nrow(states)), ]
    ancestral_states <-
      ancestral_states[, (1/3 * ncol(ancestral_states) + 1):(2/3 * ncol(ancestral_states))]
    
    rownames(ancestral_states) <- ances
    colnames(ancestral_states) <- names(mus)
    colnames(states) <- c(paste("E",names(mus), sep = '_'),
                          paste("D",names(mus), sep = '_'),
                          paste("S",names(mus), sep = '_'))
    result$ancestral_states <- ancestral_states
    result$states = states
  } 
  
  if (return_root_state) {
    root_state <- get_root_state(calcul$states,
                                 phy,
                                 mus,
                                 d)
    result$root_state <- root_state
  }
  return(result)
}

#' @title Likelihood for SecSSE model
#' Loglikelihood calculation for the SecSSE model given a set of parameters and
#' data
#' 
#' @inheritParams default_params_doc
#' @return A list with the following elements:
#' $LL the loglikelihood of the data (phylogeny + tip states) given the
#' parameters (speciation, extinction, transition rates).
#' If see_ancestral_states = TRUE, then there will be two additional elements:
#' $ancestral_states: a matrix with the probabilities of each state at the
#' internal nodes
#' $states: a matrix with the probabilities E, D (normalized) and S that are used
#' in the calculations. The ancestral_states matrix is a submatrix of this matrix.
#' This matrix is mostly used for package developers.
#' If return_root_state = TRUE, then there will be one additional element:
#' $root_state: vector with probabilities of each state at the root. This vector
#' is the same as the top row of $ancestral_states
#' We have used the shorthand description of "probabilities of each state",
#' but technically, the probabilities are the normalized probabilities D of the
#' data given each state at the internal nodes.
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
                          method = "odeint::runge_kutta_cash_karp54",
                          display_warning = TRUE,
                          use_normalization = TRUE,
                          return_root_state = FALSE) {
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
                use_normalization = use_normalization,
                return_root_state = return_root_state)
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
                              method = "odeint::runge_kutta_cash_karp54",
                              atol = 1e-8,
                              rtol = 1e-7,
                              display_warning = TRUE,
                              use_normalization = TRUE,
                              return_root_state = FALSE) {
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
                use_normalization = use_normalization,
                return_root_state = return_root_state)
}
