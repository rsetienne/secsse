#' @title Likelihood for SecSSE model
#' Loglikelihood calculation for the SecSSE model given a set of parameters and
#' data, calculated for a single branch
#' 
#' @inheritParams default_params_doc
#' @return The loglikelihood of the data given the parameter.
#' @export
secsse_single_branch_loglik <- function(parameter,
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
                                        display_warning = TRUE) {
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  q_matrix <- parameter[[3]]
  
  using_cla <- is.list(lambdas)
  
  num_modeled_traits <- ncol(q_matrix) / floor(num_concealed_states)
  
  traitStates = get_trait_states(parameter,
                                 num_concealed_states, display_warning)
  
  if (is.null(setting_calculation)) {

    check_root_state_weight(root_state_weight, traits)

    # make fake phy
    fake_phy <- ape::rphylo(n = 2, birth = 1, death = 0)
    fake_phy$edge.length[1:2] <- phy$edge.length[1]

    setting_calculation <- build_initStates_time(fake_phy,
                                                 c(traits, traits),
                                                 num_concealed_states,
                                                 sampling_fraction,
                                                 is_complete_tree,
                                                 mus,
                                                 num_modeled_traits,
                                                 traitStates = traitStates)
  } 
  
  states <- setting_calculation$states
  states <- states[-2, ]
  forTime <- setting_calculation$forTime
  forTime <- forTime[-2, ]
  
  d <- ncol(states) / 2
  
  if (!is.null(phy$root.edge)) {
      forTime[3] <- forTime[3] + phy$root.edge
  }
  
  RcppParallel::setThreadOptions(numThreads = num_threads)
  
  calcul <- calc_ll_single_branch_cpp(rhs = if (using_cla) "ode_cla" else "ode_standard",
                        states = states[1, ],
                        forTime = c(0, forTime[3]),
                        lambdas = lambdas,
                        mus = mus,
                        Q = q_matrix,
                        method = method,
                        atol = atol,
                        rtol = rtol,
                        see_states = see_ancestral_states)

  loglik <- calcul$loglik
  nodeM <- calcul$states
  mergeBranch <- calcul$merge_branch
  
  if (length(nodeM) > 2 * d) nodeM <- nodeM[1:(2 * d)]

  ## At the root
  weight_states <- get_weight_states(root_state_weight,
                                     num_concealed_states,
                                     mergeBranch,
                                     lambdas,
                                     nodeM,
                                     d,
                                     is_cla = using_cla)
  
  mergeBranch2 <- condition(cond,
                            mergeBranch,
                            weight_states,
                            lambdas,
                            nodeM,
                            is_root_edge = TRUE)
  
  wholeLike <- sum((mergeBranch2) * (weight_states))

  LL <- log(wholeLike) +
    loglik -
    penalty(pars = parameter, loglik_penalty = loglik_penalty)
  
  return(list("loglik" = LL,
              "nodeM" = nodeM,
              "merge_branch" = mergeBranch))
}