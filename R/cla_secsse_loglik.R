#' Logikelihood calculation for the cla_SecSSE model given a set of parameters
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
#' @param root_state_weight the method to weigh the states:'maddison_weigh
#' ,'proper_weights'(default) or 'equal_weights'. It can also be specified the
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
#' @note Multithreading might lead to a slightly reduced accuracy (in the order of 1e-8) and is therefore not enabled by default. Please use at your own discretion. 
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
                                  method = ifelse(num_threads == 1,
                                              "odeint::bulirsch_stoer",
                                              "odeint::runge_kutta_fehlberg78"),
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

  states <- setting_calculation$states
  forTime <- setting_calculation$forTime  # nolint
  ances <- setting_calculation$ances

  if (num_concealed_states != round(num_concealed_states)) {
    # for testing
    d <- ncol(states) / 2
    new_states <- states[, c(1:sqrt(d), (d + 1):((d + 1) + sqrt(d) - 1))]
    new_states <- states[, c(1, 2, 3, 10, 11, 12)]
    states <- new_states
  }

  loglik <- 0
  d <- ncol(states) / 2

  calcul <- c()
  if (num_threads == 1) {
    ancescpp <- ances - 1
    forTimecpp <- forTime # nolint
    forTimecpp[, c(1, 2)] <- forTimecpp[, c(1, 2)] - 1 # nolint
      calcul <- cla_calThruNodes_cpp(ancescpp,
                                     states,
                                     forTimecpp,
                                     lambdas,
                                     mus,
                                     Q,
                                     method,
                                     atol,
                                     rtol,
                                     is_complete_tree)
  } else {
    # because C++ indexes from 0, we need to adjust the indexing:
    ancescpp <- ances - 1
    forTimecpp <- forTime # nolint
    forTimecpp[, c(1, 2)] <- forTimecpp[, c(1, 2)] - 1 # nolint

    if (num_threads == -2) {
      calcul <- calc_cla_ll_threaded(ancescpp,
                                     states,
                                     forTimecpp,
                                     lambdas,
                                     mus,
                                     Q,
                                     1,
                                     method,
                                     is_complete_tree)
    } else {
      calcul <- calc_cla_ll_threaded(ancescpp,
                                     states,
                                     forTimecpp,
                                     lambdas,
                                     mus,
                                     Q,
                                     num_threads,
                                     method,
                                     is_complete_tree)
    }
  }

  mergeBranch <- calcul$mergeBranch # nolint
  nodeM <- calcul$nodeM  # nolint
  loglik <- calcul$loglik

  ## At the root
  mergeBranch2 <- mergeBranch # nolint
  lmb <- length(mergeBranch2)
  if (is.numeric(root_state_weight)) {
    weightStates <- rep(root_state_weight / num_concealed_states, # nolint
                        num_concealed_states)
  } else {
    if (root_state_weight == "maddison_weights") {
      weightStates <- mergeBranch / sum(mergeBranch2)
    }
    if (root_state_weight == "proper_weights") {
      numerator <- rep(NA, lmb)
      for (j in 1:lmb) {
        numerator[j] <- mergeBranch2[j] / sum(lambdas[[j]] *
                                  ((1 - nodeM[1:d]) %o% (1 - nodeM[1:d])))
      }
      weightStates <- numerator / sum(numerator) # nolint
    }
    if (root_state_weight == "equal_weights") {
      weightStates <- rep(1 / lmb, lmb) # nolint
    }
  }

  if (cond == "maddison_cond") {
    preCond <- rep(NA, lmb) # nolint
    for (j in 1:lmb) {
      if (root_state_weight == "equal_weights") {
        weightStates <- rep(1 / lmb, lmb) # nolint
      }
    }
  }

  if (cond == "maddison_cond") {
    preCond <- rep(NA, lmb) # nolint
    for (j in 1:lmb) {
      preCond[j] <- sum(weightStates[j] *
                        lambdas[[j]] *
                        (1 - nodeM[1:d][j]) ^ 2)
     }
    mergeBranch2 <- mergeBranch2 / sum(preCond) # nolint
  }

  if (is_complete_tree) {
    timeInte <- max(abs(ape::branching.times(phy))) # nolint
    y <- rep(0, lmb)

    nodeM <- ct_condition_cla(y, # nolint
                              timeInte,
                              lambdas,
                              mus,
                              Q,
                              "odeint::bulirsch_stoer",
                              1e-16,
                              1e-12)
    nodeM <- c(nodeM, y) # nolint
  }

  if (cond == "proper_cond") {
    preCond <- rep(NA, lmb) # nolint
    for (j in 1:lmb) {
      preCond[j] <- sum(lambdas[[j]] * ((1 - nodeM[1:d]) %o% (1 - nodeM[1:d])))
    }
    mergeBranch2 <- mergeBranch2 / preCond # nolint
  }

  wholeLike_atRoot <- sum(mergeBranch2 * weightStates) # nolint
  LL <- log(wholeLike_atRoot) + # nolint
        loglik -
        penalty(pars = parameter,
                loglik_penalty = loglik_penalty)

  if (see_ancestral_states == TRUE) {
    num_tips <- ape::Ntip(phy)
    ancestral_states <- states[(num_tips + 1):nrow(states), ]
    ancestral_states <- 
        ancestral_states[, -1 * (1:(ncol(ancestral_states) / 2))]
    rownames(ancestral_states) <- ances
    return(list(ancestral_states = ancestral_states, LL = LL))
  } else {
    return(LL)
  }
}
