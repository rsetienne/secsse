#' Logikelihood calculation for the cla_SecSSE model given a set of parameters and data using Rcpp
#' @title Likelihood for SecSSE model, using Rcpp
#' @param parameter list where the first is a table where lambdas across different modes of speciation are shown, the second mus and the third transition rates.
#' @param phy phylogenetic tree of class phylo, ultrametric, fully-resolved, rooted and with branch lengths.
#' @param traits vector with trait states, order of states must be the same as tree tips, for help, see vignette.
#' @param num_concealed_states number of concealed states, generally equivalent to number of examined states.
#' @param use_fortran Should the Fortran code for numerical integration be called? Default is TRUE.
#' @param methode Solver for the set of equations, default is 'ode45'.
#' @param cond condition on the existence of a node root: 'maddison_cond','proper_cond'(default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:'maddison_weights','proper_weights'(default) or 'equal_weights'. It can also be specified the root state:the vector c(1,0,0) indicates state 1 was the root state.
#' @param sampling_fraction vector that states the sampling proportion per trait state. It must have as many elements as trait states.
#' @param run_parallel should the routine to run in parallel be called?
#' @param setting_calculation argument used internally to speed up calculation. It should be leave blank (default : setting_calculation = NULL)
#' @param setting_parallel argument used internally to set a parallel calculation. It should be left blank (default : setting_parallel = NULL)
#' @param see_ancestral_states should the ancestral states be shown? Deafault FALSE
#' @param loglik_penalty the size of the penalty for all parameters; default is 0 (no penalty)
#' @param is_complete_tree whether or not a tree with all its extinct species is provided
#' @param func function to be used in solving the ODE system. Currently only for testing purposes.
#' @note To run in parallel it is needed to load the following libraries when windows: apTreeshape, doparallel and foreach. When unix, it requires: apTreeshape, doparallel, foreach and doMC
#' @return The loglikelihood of the data given the parameters
#' @examples
#'rm(list=ls(all=TRUE))
#'library(secsse)
#'library(DDD)
#'library(deSolve)
#'#library(diversitree)
#'library(apTreeshape)
#'library(foreach)
#'set.seed(13)
#'phylotree <- ape::rcoal(12, tip.label = 1:12)
#'traits <- sample(c(0,1,2),ape::Ntip(phylotree),replace=TRUE)
#'num_concealed_states <- 3
#'sampling_fraction <- c(1,1,1)
#'methode <- 'ode45'
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
#'# which is a list with 9 matrices, corresponding to the 9 states (0A,1A,2A,0B,etc)
#'# if we want to calculate a single likelihood:
#'parameter <- idparlist
#'lambd_and_modeSpe <- parameter$lambdas
#'lambd_and_modeSpe[1,] <- c(0.2,0.2,0.2,0.4,0.4,0.4,0.01,0.01,0.01)
#'parameter[[1]] <- prepare_full_lambdas(traits,num_concealed_states,lambd_and_modeSpe)
#'parameter[[2]] <- rep(0,9)
#'masterBlock <- matrix(0.07,ncol=3,nrow=3,byrow=TRUE)
#'diag(masterBlock) <- NA
#'parameter [[3]] <- q_doubletrans(traits,masterBlock,diff.conceal = FALSE)
#'cla_secsse_loglik_cpp(parameter, phy, traits, num_concealed_states,
#'                  use_fortran = FALSE, methode = 'ode45', cond = 'maddison_cond',
#'                  root_state_weight = 'maddison_weights', sampling_fraction,
#'                  run_parallel = FALSE, setting_calculation = NULL,
#'                  setting_parallel = NULL, see_ancestral_states = FALSE,
#'                  loglik_penalty = 0)
#'# LL = -37.8741
#' @export
cla_secsse_loglik_cpp <- function(parameter,
                              phy,
                              traits,
                              num_concealed_states,
                              use_fortran = FALSE,
                              methode = "odeint::runge_kutta_dopri5",
                              cond = "proper_cond",
                              root_state_weight = "proper_weights",
                              sampling_fraction,
                              run_parallel = FALSE,
                              setting_calculation = NULL,
                              setting_parallel = NULL,
                              see_ancestral_states = FALSE,
                              loglik_penalty = 0,
                              is_complete_tree = FALSE,
                              func =
                                ifelse(is_complete_tree, "cla_secsse_runmod_ct",
                                       ifelse(use_fortran == FALSE,
                                              cla_secsse_loglik_rhs,
                                              "cla_secsse_runmod"))) {
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  Q <- parameter[[3]]


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
    forTime <- setting_calculation$forTime
    ances <- setting_calculation$ances

    if (num_concealed_states != round(num_concealed_states)) {
      # for testing
      d <- ncol(states)/2
      new_states <- states[, c(1:sqrt(d), (d + 1):((d + 1) + sqrt(d) - 1))]
      new_states <- states[, c(1, 2, 3, 10, 11, 12)]
      states <- new_states
    }

    loglik <- 0
    ly <- ncol(states)
    d <- ncol(states)/2

    calcul <- cla_calThruNodes_cpp(ances,
                                   states,
                                   forTime,
                                   lambdas,
                                   mus,
                                   Q,
                                   methode)

  mergeBranch <- calcul$mergeBranch
  nodeM <- calcul$nodeM
  loglik <- calcul$loglik

  ## At the root
  mergeBranch2 <- (mergeBranch)
  if (is.numeric(root_state_weight)) {
    giveWeights <- root_state_weight/num_concealed_states
    weightStates <- rep(giveWeights, num_concealed_states)

  } else {
    if (root_state_weight == "maddison_weights") {
      weightStates <- (mergeBranch2)/sum((mergeBranch2))
    }

    if (root_state_weight == "proper_weights") {
      numerator <- NULL
      for (j in 1:length(mergeBranch2)) {
        numerator <- c(numerator,
                       (mergeBranch2[j] /
                          (sum(lambdas[[j]] * (1 - nodeM[1:d][j])^2))))
      }
      denomin <- NULL
      for (j in 1:length(mergeBranch2)) {
        denomin <- c(denomin,
                     (mergeBranch2[j] /
                        (sum(lambdas[[j]] * (1 - nodeM[1:d][j])^2))))
      }

      weightStates <- numerator/sum(denomin)
    }
    if (root_state_weight == "equal_weights") {
      weightStates <- rep(1 / length(mergeBranch2), length(mergeBranch2))
    }
  }

  if (cond == "maddison_cond") {
    preCond <- NULL
    for (j in 1:length(weightStates)) {
      preCond <- c(preCond,
                   sum(weightStates[j] *
                         lambdas[[j]] *
                         (1 - nodeM[1:d][j])^2))
    }
    mergeBranch2 <- mergeBranch2/(sum(preCond))
  }

  if (cond == "proper_cond") {
    preCond <- NULL
    for (j in 1:length(mergeBranch2)) {
      preCond <- c(preCond,
                   sum((lambdas[[j]] * (1 - nodeM[1:d][j])^2)))
    }
    mergeBranch2 <- mergeBranch2/preCond
  }

  atRoot <- ((mergeBranch2) * (weightStates))

  wholeLike <- sum(atRoot)
  LL <- log(wholeLike) + loglik - penalty(pars = parameter,
                                          loglik_penalty = loglik_penalty)

  if (see_ancestral_states == TRUE) {
    num_tips <- ape::Ntip(phy)
    ancestral_states <- states[(num_tips + 1):nrow(states), ]
    ancestral_states <- ancestral_states[, -(1:(ncol(ancestral_states)/2))]
    rownames(ancestral_states) <- ances
    return(list(ancestral_states = ancestral_states, LL = LL))
  } else {
    return(LL)
  }
}
