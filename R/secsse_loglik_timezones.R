#' Logikelihood calculation for the SecSSE model given a set of parameters and 
#' data, using different parameters before and after a critical timepoint.
#' @title Likelihood for SecSSE model
#' @param parameter list where the first vector represents lambdas before,
#' the second mus before and the third transition rates before. The fourth 
#' vector represents lambdas after, the fifth vector mus after and the sixth
#' vector transition rates after.
#' @param phy phylogenetic tree of class phylo, ultrametric, fully-resolved, 
#' rooted and with branch lengths.
#' @param traits vector with trait states, order of states must be the same as 
#' tree tips, for help, see vignette.
#' @param num_concealed_states number of concealed states, generally equivalent 
#' to number of examined states.
#' @param critical_t the critical time point at which all the rates change.
#' @param cond condition on the existence of a node root: "maddison_cond",
#' "proper_cond"(default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:"maddison_weights","proper_weights"(default) or "equal_weights". It can also be specified the 
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
#' @param num_threads number of threads. Set to -1 to use all available threads. 
#' Default is one thread.
#' @param atol absolute tolerance of integration
#' @param rtol relative tolerance of integration
#' @param method integration method used, available are: 
#' "odeint::runge_kutta_cash_karp54", "odeint::runge_kutta_fehlberg78", 
#' "odeint::runge_kutta_dopri5", "odeint::bulirsch_stoer" and 
#' "odeint::runge_kutta4". Default method is:"odeint::bulirsch_stoer".
#' @return The loglikelihood of the data given the parameters
#' @note Multithreading might lead to a slightly reduced accuracy (in the order of 1e-10) and is therefore not enabled by default. Please use at your own discretion. 
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
#' run_parallel <- FALSE
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
secsse_loglik_timezones <- function(parameter,
                                    phy,
                                    traits,
                                    num_concealed_states,
                                    critical_t,
                                    cond = "proper_cond",
                                    root_state_weight = "proper_weights",
                                    sampling_fraction,
                                    setting_calculation = NULL,
                                    see_ancestral_states = FALSE,
                                    loglik_penalty = 0,
                                    is_complete_tree = FALSE,
                                    num_threads = 1,
                                    atol = 1e-12,
                                    rtol = 1e-12,
                                    method = "odeint::bulirsch_stoer") {
  lambdas1 <- parameter[[1]]
  mus1 <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  Q1 <- parameter[[3]]
  
  lambdas2 <- parameter[[1]]
  mus2 <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  Q2 <- parameter[[3]]
  
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
                                                 mus1)
  }
  
  states <- setting_calculation$states
  
  forTime <- setting_calculation$forTime
  ances <- setting_calculation$ances
  
  if (num_concealed_states != round(num_concealed_states)) { # for testing
    d <- ncol(states) / 2
    new_states <- states[,c(1:sqrt(d),(d + 1):((d + 1) + sqrt(d) - 1))]
    new_states <- states[,c(1,2,3,10,11,12)]
    states <- new_states
  }
  
  ly <- ncol(states)
  d <- ncol(states) / 2
  
  
  calcul <- calThruNodes_timezones_cpp(ances,
                                       states,
                                       forTime,
                                       lambdas1,
                                       mus1,
                                       Q1,
                                       lambdas2,
                                       mus2,
                                       Q2,
                                       critical_t,
                                       1,
                                       atol,
                                       rtol,
                                       method,
                                       is_complete_tree)
  
  loglik <- calcul$loglik
  nodeM <- calcul$nodeM
  mergeBranch <- calcul$mergeBranch
  
  if (length(nodeM) > 2 * d) nodeM <- nodeM[1:(2 * d)]
  
  ## At the root
  mergeBranch2 <- (mergeBranch)
  if (is.numeric(root_state_weight)) {
    weightStates <- rep(root_state_weight / num_concealed_states, 
                        num_concealed_states)
  } else {
    if (root_state_weight == "maddison_weights") {
      weightStates <- (mergeBranch2)/sum((mergeBranch2))
    }
    
    if (root_state_weight == "proper_weights") {
      weightStates <- (mergeBranch2 / 
                         (lambdas * (1 - nodeM[1:d]) ^ 2)) / 
        sum((mergeBranch2 / (lambdas * (1 - nodeM[1:d]) ^ 2)))
    }
    
    if (root_state_weight == "equal_weights") {
      weightStates <- rep(1 / length(mergeBranch2),
                          length(mergeBranch2))
    }
  }
  
  if (is_complete_tree) {
    warning("Complete tree conditioning not supported for multiple timezones")
  }
  
  if (cond == "maddison_cond") {
    mergeBranch2 <- 
      mergeBranch2 / sum(weightStates * lambdas * (1 - nodeM[1:d]) ^ 2)
  }
  
  if (cond == "proper_cond") {
    mergeBranch2 <- mergeBranch2 / (lambdas * (1 - nodeM[1:d]) ^ 2)
  }
  
  
  wholeLike <- sum((mergeBranch2) * (weightStates))
  LL <- log(wholeLike) + loglik - penalty(pars = parameter,
                                          loglik_penalty = loglik_penalty)
  
  if (see_ancestral_states == TRUE) {
    num_tips <- ape::Ntip(phy)
    ancestral_states <- states[(num_tips + 1):nrow(states), ]
    ancestral_states <- ancestral_states[, -(1:(ncol(ancestral_states) / 2))]
    rownames(ancestral_states) <- ances
    return(list(ancestral_states = ancestral_states,LL = LL))
  } else {
    return(LL)
  }
}
  