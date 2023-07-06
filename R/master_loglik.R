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
                          method = "odeint::bulirsch_stoer") {
  
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  q_matrix <- parameter[[3]]
  
  using_cla <- FALSE
  if (is.list(lambdas)) using_cla <- TRUE
  
  num_modeled_traits <- ncol(q_matrix) / floor(num_concealed_states)
  
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
                                                 num_modeled_traits)
  } else {
    # with a complete tree, we need to re-calculate the states every time we
    # run, because they are dependent on mu.
    if (is_complete_tree) {
      states <- build_states(phy = phy,
                             traits = traits,
                             num_concealed_states = num_concealed_states,
                             sampling_fraction = sampling_fraction,
                             is_complete_tree = is_complete_tree,
                             mus = mus)
    }
  }
  
  states <- setting_calculation$states
  forTime <- setting_calculation$forTime
  ances <- setting_calculation$ances
  
  d <- ncol(states) / 2
  
  if (see_ancestral_states == TRUE && num_threads != 1) {
    warning("see ancestral states only works with one thread, 
              setting to one thread")
    num_threads <- 1
  }
  
  calcul <- update_using_cpp(ances,
                             states,
                             forTime,
                             lambdas,
                             mus,
                             q_matrix,
                             method,
                             atol,
                             rtol,
                             is_complete_tree,
                             num_threads)
  
  loglik <- calcul$loglik
  nodeM <- calcul$nodeM
  mergeBranch <- calcul$mergeBranch
  states <- calcul$states
  
  if (length(nodeM) > 2 * d) nodeM <- nodeM[1:(2 * d)]
  
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
                                                      length(mergeBranch))
  
  mergeBranch2 <- condition(cond,
                            mergeBranch,
                            weight_states,
                            lambdas,
                            nodeM)
  
  wholeLike <- sum((mergeBranch2) * (weight_states))
  
  LL <- log(wholeLike) +
    loglik -
    penalty(pars = parameter, loglik_penalty = loglik_penalty)
  
  if (see_ancestral_states == TRUE) {
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
