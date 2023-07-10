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

  using_cla <- is.list(lambdas)

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
                        see_states = see_ancestral_states)
  loglik <- calcul$loglik
  nodeM <- calcul$node_M
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

#' Loglikelihood calculation for the SecSSE model given a set of parameters and
#' data
#' @title Likelihood for SecSSE model
#' @param parameter list where first vector represents lambdas, the second mus
#' and the third transition rates.
#' @param phy phylogenetic tree of class phylo, ultrametric, fully-resolved,
#' rooted and with branch lengths.
#' @param traits vector with trait states, order of states must be the same as
#' tree tips, for help, see vignette.
#' @param num_concealed_states number of concealed states, generally equivalent
#' to number of examined states.
#' @param cond condition on the existence of a node root: "maddison_cond",
#' "proper_cond"(default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:
#' "maddison_weights" or "proper_weights"(default).
#' It can also be specified the root state:the vector c(1, 0, 0)
#' indicates state 1 was the root state.
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
#' @param num_threads number of threads. Set to -1 to use all available threads.
#' Default is one thread.
#' @param atol absolute tolerance of integration
#' @param rtol relative tolerance of integration
#' @param method integration method used, available are:
#' "odeint::runge_kutta_cash_karp54", "odeint::runge_kutta_fehlberg78",
#' "odeint::runge_kutta_dopri5", "odeint::bulirsch_stoer" and
#' "odeint::runge_kutta4". Default method is:"odeint::bulirsch_stoer".
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
                          num_threads = 1,
                          atol = 1e-8,
                          rtol = 1e-7,
                          method = "odeint::bulirsch_stoer") {
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
                num_threads = num_threads,
                atol = atol,
                rtol = rtol,
                method = method)
}

#' Loglikelihood calculation for the cla_SecSSE model given a set of parameters
#' and data using Rcpp
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
#' @param cond condition on the existence of a node root: 'maddison_cond',
#' 'proper_cond'(default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:'maddison_weights
#'  or 'proper_weights'(default) oIt can also be specified the
#' root state:the vector c(1,0,0) indicates state 1 was the root state.
#' @param sampling_fraction vector that states the sampling proportion per trait
#' state. It must have as many elements as trait states.
#' @param setting_calculation argument used internally to speed up calculation.
#' It should be leave blank (default : setting_calculation = NULL)
#' @param see_ancestral_states should the ancestral states be shown? Deafault
#' FALSE
#' @param loglik_penalty the size of the penalty for all parameters; default is
#' 0 (no penalty)
#' @param is_complete_tree whether or not a tree with all its extinct species is
#' provided
#' @param num_threads number of threads to be used, default is 1. Set to -1 to
#' use all available threads.
#' @param method integration method used, available are:
#' "odeint::runge_kutta_cash_karp54", "odeint::runge_kutta_fehlberg78",
#' "odeint::runge_kutta_dopri5", "odeint::bulirsch_stoer" and
#' "odeint::runge_kutta4". Default method is:"odeint::bulirsch_stoer".
#' @param atol absolute tolerance of integration
#' @param rtol relative tolerance of integration
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
                              num_threads = 1,
                              method = "odeint::bulirsch_stoer",
                              atol = 1e-8,
                              rtol = 1e-7) {
  master_loglik(parameter,
                phy,
                traits,
                num_concealed_states,
                cond,
                root_state_weight,
                sampling_fraction,
                setting_calculation,
                see_ancestral_states,
                loglik_penalty,
                is_complete_tree,
                num_threads,
                atol,
                rtol,
                method)
}
