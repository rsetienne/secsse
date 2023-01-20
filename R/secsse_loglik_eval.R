#' Logikelihood calculation for the SecSSE model given a set of parameters and
#' data, returning also the likelihoods along the branches
#' @title Likelihood for SecSSE model
#' @param parameter list where first vector represents lambdas, the second mus
#' and the third transition rates.
#' @param phy phylogenetic tree of class phylo, ultrametric, fully-resolved,
#' rooted and with branch lengths.
#' @param traits vector with trait states, order of states must be the same as
#' tree tips, for help, see vignette.
#' @param num_concealed_states number of concealed states, generally equivalent
#' to number of examined states.
#' @param ancestral_states ancestral states matrix provided by
#' secsse_loglik, this is used as starting points for the branch integration
#' @param cond condition on the existence of a node root: "maddison_cond",
#' "proper_cond"(default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:"maddison_weights",
#' "proper_weights"(default) or "equal_weights". It can also be specified the
#' root state:the vector c(1,0,0) indicates state 1 was the root state.
#' @param sampling_fraction vector that states the sampling proportion per
#' trait state. It must have as many elements as trait states.
#' @param setting_calculation argument used internally to speed up calculation.
#' It should be left blank (default : setting_calculation = NULL)
#' @param see_ancestral_states should the ancestral states be shown? Default
#' FALSE
#' @param loglik_penalty the size of the penalty for all parameters; default is
#' 0 (no penalty)
#' @param is_complete_tree whether or not a tree with all its extinct species
#' is provided
#' @param atol absolute tolerance of integration
#' @param rtol relative tolerance of integration
#' @param method integration method used, available are:
#' "odeint::runge_kutta_cash_karp54", "odeint::runge_kutta_fehlberg78",
#' "odeint::runge_kutta_dopri5", "odeint::bulirsch_stoer" and
#' "odeint::runge_kutta4". Default method is:"odeint::bulirsch_stoer".
#' @param num_steps number of substeps to show intermediate likelihoods
#' along a branch, if left to NULL, the intermediate likelihoods at every
#' integration evaluation are stored, which is more exact, but can lead to 
#' huge datasets / memory usage.
#' @return The loglikelihood of the data given the parameters
#' @examples
#' set.seed(5)
#' focal_tree <- ape::rphylo(n = 4, birth = 1, death = 0)
#' traits <- c(0, 1, 1, 0)
#' params <- secsse::id_paramPos(c(0, 1), 2)
#' params[[1]][] <- c(0.2, 0.2, 0.1, 0.1)
#' params[[2]][] <- 0.0
#' params[[3]][, ] <- 0.1
#' diag(params[[3]]) <- NA
#' #  Thus, we have for both, rates
#' # 0A, 1A, 0B and 1B. If we are interested in the posterior probability of trait 0,
#' # we have to provide a helper function that sums the probabilities of 0A and 0B, e.g.: 
#' helper_function <- function(x) {
#'   return(sum(x[c(5, 7)]) / sum(x)) # normalized by total sum, just in case.
#' }
#' 
#' secsse_loglik_eval(parameters = params,
#'                    phy = focal_tree,
#'                    traits = traits,
#'                    num_concealed_states = 2,
#'                    sampling_fraction = c(1, 1),
#'                    steps = 10)
#' @export
secsse_loglik_eval <- function(parameter,
                               phy,
                               traits,
                               num_concealed_states,
                               ancestral_states,
                               cond = "proper_cond",
                               root_state_weight = "proper_weights",
                               sampling_fraction,
                               setting_calculation = NULL,
                               loglik_penalty = 0,
                               is_complete_tree = FALSE,
                               atol = 1e-12,
                               rtol = 1e-12,
                               method = "odeint::bulirsch_stoer",
                               num_steps = NULL) {
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  q_matrix <- parameter[[3]]

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

  for_time <- setting_calculation$forTime
  ances <- setting_calculation$ances


  calcul <- calThruNodes_store_cpp(ances,
                                   ancestral_states,
                                   for_time,
                                   lambdas,
                                   mus,
                                   q_matrix,
                                   1,
                                   atol,
                                   rtol,
                                   method,
                                   is_complete_tree,
                                   ifelse(is.null(num_steps), 0, num_steps))
  # if the number of steps == NULL, pass a 0.
  return(calcul)
}
