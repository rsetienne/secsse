#' function to provide probabilities of observing states along branches. 
#' @title Likelihood for SecSSE model, using Rcpp
#' @param parameter list where the first is a table where lambdas across
#' different modes of speciation are shown, the second mus and the third
#'  transition rates.
#' @param phy phylogenetic tree of class phylo, ultrametric, fully-resolved,
#' rooted and with branch lengths.
#' @param traits vector with trait states, order of states must be the same as
#'  tree tips, for help, see vignette.
#' @param num_concealed_states number of concealed states, generally equivalent
#' to number of examined states.
#' @param ancestral_states ancestral states matrix provided by 
#' cla_secsse_loglik, this is used as starting points for manual integration
#' @param num_steps number of steps to integrate along a branch
#' @param cond condition on the existence of a node root: 'maddison_cond',
#' 'proper_cond'(default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:'maddison_weigh
#' ,'proper_weights'(default) or 'equal_weights'. It can also be specified the
#' root state:the vector c(1,0,0) indicates state 1 was the root state.
#' @param sampling_fraction vector that states the sampling proportion per trait
#' state. It must have as many elements as trait states.
#' @param setting_calculation argument used internally to speed up calculation.
#' It should be leave blank (default : setting_calculation = NULL)
#' @param loglik_penalty the size of the penalty for all parameters; default is
#' 0 (no penalty)
#' @param is_complete_tree whether or not a tree with all its extinct species is
#' provided
#' @param method integration method used, available are: 
#' "odeint::runge_kutta_cash_karp54", "odeint::runge_kutta_fehlberg78", 
#' "odeint::runge_kutta_dopri5", "odeint::bulirsch_stoer" and 
#' "odeint::runge_kutta4". Default method is:"odeint::bulirsch_stoer".
#' @param atol absolute tolerance of integration
#' @param rtol relative tolerance of integration
#' @return The loglikelihood of the data given the parameters
#' @description Using see_ancestral_states = TRUE in the function 
#' cla_secsse_loglik will provide posterior probabilities of the states of the 
#' model on the nodes of the tree, but will not give the values on the branches.
#' This function evaluates these probabilities at fixed time intervals dt. 
#' Because dt is fixed, this may lead to some inaccuracies, and dt is best 
#' chosen as small as possible. 
#' @export
cla_secsse_eval <- function(parameter,
                            phy,
                            traits,
                            num_concealed_states,
                            ancestral_states,
                            num_steps = 10,
                            cond = "proper_cond",
                            root_state_weight = "proper_weights",
                            sampling_fraction,
                            setting_calculation = NULL,
                            see_ancestral_states = FALSE,
                            loglik_penalty = 0,
                            is_complete_tree = FALSE,
                            method = "odeint::bulirsch_stoer",
                            atol = 1e-16,
                            rtol = 1e-16) {
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  Q <- parameter[[3]]  # nolint
  
  
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
                                                 mus)
  }
  
  forTime <- setting_calculation$forTime  # nolint
  ances <- setting_calculation$ances
  
  calcul <- c()
  ancescpp <- ances - 1
  forTimecpp <- forTime # nolint
  forTimecpp[, c(1, 2)] <- forTimecpp[, c(1, 2)] - 1 # nolint
  calcul <- cla_calThruNodes_store_cpp(ancescpp,
                                       ancestral_states,
                                       forTimecpp,
                                       lambdas,
                                       mus,
                                       Q,
                                       method,
                                       atol,
                                       rtol,
                                       num_steps,
                                       is_complete_tree)
  
  return(calcul)
}