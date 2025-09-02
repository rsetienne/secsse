#' @title Test Likelihood for SecSSE model, using Rcpp
#' Loglikelihood calculation for the cla_SecSSE model given a set of parameters
#' and data using Rcpp
#' 
#' @inheritParams default_params_doc
#' 
#' @return The loglikelihood of the data given the parameters
#
#' @export
test_loglik <- function(parameter,
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
                        use_normalization = TRUE) {
  
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
  
  E <- nodeM[1:d]
  S <- nodeM[(2 * d + 1):(3 * d)]
  
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
  return(list("LL" = LL,
              "E" = E,
              "S" = S))
}