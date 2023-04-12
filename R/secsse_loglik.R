#' Logikelihood calculation for the SecSSE model given a set of parameters and
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
#' "maddison_weights","proper_weights"(default) or "equal_weights".
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
#' @note Multithreading might lead to a slightly reduced accuracy
#' (in the order of 1e-10) and is therefore not enabled by default.
#' Please use at your own discretion.
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
                          atol = 1e-12,
                          rtol = 1e-12,
                          method = "odeint::bulirsch_stoer") {
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

  states <- setting_calculation$states

  if (num_concealed_states != round(num_concealed_states)) { # for test case
    d <- ncol(states) / 2
    new_states <- states[, c(1:sqrt(d), (d + 1):((d + 1) + sqrt(d) - 1))]
    new_states <- states[, c(1, 2, 3, 10, 11, 12)]
    states <- new_states
  }


  if (is_complete_tree) {
    states <- build_states(phy = phy,
                           traits = traits,
                           num_concealed_states = num_concealed_states,
                           sampling_fraction = sampling_fraction,
                           is_complete_tree = is_complete_tree,
                           mus = mus)
  }

  forTime <- setting_calculation$forTime
  ances <- setting_calculation$ances

  d <- ncol(states) / 2

  if (see_ancestral_states == TRUE) {
    if (num_threads != 1) {
      warning("see ancestral states only works with one thread, 
              setting to one thread")
      num_threads <- 1
    }
  }
  calcul <- c()
  if (num_threads == 1) {
    calcul <- calThruNodes_cpp(ances,
                               states,
                               forTime,
                               lambdas,
                               mus,
                               q_matrix,
                               1,
                               atol,
                               rtol,
                               method,
                               is_complete_tree)
  } else {
    ancescpp <- ances - 1
    forTimecpp <- forTime
    forTimecpp[, c(1, 2)] <- forTimecpp[, c(1, 2)] - 1

    if (num_threads == -2) {
      calcul <- calc_ll_threaded(lambdas,
                                 mus,
                                 q_matrix,
                                 ancescpp,
                                 forTimecpp,
                                 states,
                                 1,
                                 method,
                                 is_complete_tree)
    } else {
      calcul <- calc_ll_threaded(lambdas,
                                 mus,
                                 q_matrix,
                                 ancescpp,
                                 forTimecpp,
                                 states,
                                 num_threads,
                                 method,
                                 is_complete_tree)
    }
  }

  loglik <- calcul$loglik
  nodeM <- calcul$nodeM
  mergeBranch <- calcul$mergeBranch
  states <- calcul$states

  if (length(nodeM) > 2 * d) nodeM <- nodeM[1:(2 * d)]

  ## At the root
  mergeBranch2 <- (mergeBranch)

  weightStates <- get_weight_states(root_state_weight,
                                    num_concealed_states,
                                    mergeBranch,
                                    lambdas,
                                    nodeM,
                                    d,
                                    is_cla = FALSE)

  if (is_complete_tree) {
    time_inte <- max(abs(ape::branching.times(phy))) # nolint
    y <- rep(0, 2 * length(mergeBranch2))

    nodeM <- ct_condition(y, # nolint
                          time_inte,
                          lambdas,
                          mus,
                          q_matrix,
                          method,
                          atol,
                          rtol)
  }

  if (cond == "maddison_cond") {
    mergeBranch2 <-
      mergeBranch2 / sum(weightStates * lambdas * (1 - nodeM[1:d]) ^ 2)
  }

  if (cond == "proper_cond") {
    mergeBranch2 <- mergeBranch2 / (lambdas * (1 - nodeM[1:d]) ^ 2)
  }

  wholeLike <- sum((mergeBranch2) * (weightStates))
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

#' @keywords internal
check_tree <- function(phy, is_complete_tree) {
  if (ape::is.rooted(phy) == FALSE) {
    stop("The tree needs to be rooted.")
  }

  if (ape::is.binary(phy) == FALSE) {
    stop("The tree needs to be fully resolved.")
  }
  if (ape::is.ultrametric(phy) == FALSE && is_complete_tree == FALSE) {
    stop("The tree needs to be ultrametric.")
  }
  if (any(phy$edge.length == 0)) {
    stop("The tree must have internode distancs that are all larger than 0.")
  }
}

check_traits <- function(traits, sampling_fraction) {
  if (is.matrix(traits)) {
    if (length(sampling_fraction) != length(sort(unique(traits[, 1])))) {
      stop("Sampling_fraction must have as many elements 
           as the number of traits.")
    }

    if (all(sort(unique(as.vector(traits))) == sort(unique(traits[, 1]))) ==
        FALSE) {
      stop(
        "Check your trait argument; if you have more than one column,
        make sure all your states are included in the first column."
      )
    }
  } else {
    if (length(sampling_fraction) != length(sort(unique(traits)))) {
      stop("Sampling_fraction must have as many elements as 
           the number of traits.")
    }
  }

  if (length(sort(unique(as.vector(traits)))) < 2) {
    stop("The trait has only one state.")
  }
}

check_root_state_weight <- function(root_state_weight, traits) {
  if (is.numeric(root_state_weight)) {
    if (length(root_state_weight) != length(sort(unique(traits)))) {
      stop("There need to be as many elements in root_state_weight 
           as there are traits.")
    }
    if (length(which(root_state_weight == 1)) != 1) {
      stop("The root_state_weight needs only one 1.")
    }
  } else {
    if (any(root_state_weight == "maddison_weights" |
            root_state_weight == "equal_weights" |
            root_state_weight == "proper_weights") == FALSE) {
      stop("The root_state_weight must be any of 
           maddison_weights, equal_weights, or proper_weights.")
    }
  }
}

check_input <- function(traits,
                        phy,
                        sampling_fraction,
                        root_state_weight,
                        is_complete_tree) {
  check_root_state_weight(root_state_weight, sampling_fraction)

  check_tree(phy, is_complete_tree)

  check_traits(traits, sampling_fraction)
}

create_states <- function(usetraits,
                          states,
                          sampling_fraction,
                          num_concealed_states,
                          d,
                          traitStates,
                          is_complete_tree,
                          phy,
                          ly,
                          mus,
                          nb_tip) {
  if (anyNA(usetraits)) {
    nas <- which(is.na(traits))
    for (iii in seq_along(nas)) {
      states[nas[iii], ] <- c(1 - rep(sampling_fraction,
                                      num_concealed_states),
                              rep(sampling_fraction, num_concealed_states))
    }
  }

  for (iii in seq_along(traitStates)) { # Initial state probabilities
    StatesPresents <- d + iii
    toPlaceOnes <- StatesPresents +
      length(traitStates) * (0:(num_concealed_states - 1))
    tipSampling <- 1 * sampling_fraction
    states[which(usetraits ==
                   traitStates[iii]), toPlaceOnes] <- tipSampling[iii]
  }

  if (is_complete_tree) {
    extinct_species <- geiger::is.extinct(phy)
    if (!is.null(extinct_species)) {
      for (i in seq_along(extinct_species)) {
        states[which(phy$tip.label == extinct_species[i]), (d + 1):ly] <-
          mus * states[which(phy$tip.label == extinct_species[i]), (d + 1):ly]
      }
    }
    for (iii in 1:nb_tip) {
      states[iii, 1:d] <- 0
    }
  } else {
    for (iii in 1:nb_tip) {
      states[iii, 1:d] <- rep(1 - sampling_fraction, num_concealed_states)
    }
  }

  return(states)
}


build_states <- function(phy,
                         traits,
                         num_concealed_states,
                         sampling_fraction,
                         is_complete_tree = FALSE,
                         mus = NULL) {
  if (!is.matrix(traits)) {
    traits <- matrix(traits, nrow = length(traits), ncol = 1, byrow = FALSE)
  }

  if (length(phy$tip.label) != nrow(traits)) {
    stop("Number of species in the tree must be the same as in the trait file")
  }
  traitStates <- sort(unique(traits[, 1]))

  nb_tip <- ape::Ntip(phy)
  nb_node <- phy$Nnode
  ly <- length(traitStates) * 2 * num_concealed_states
  states <- matrix(ncol = ly, nrow = nb_tip + nb_node)
  d <- ly / 2
  ## In a example of 3 states, the names of the colums would be like:
  ##
  ## colnames(states) <- c("E0A","E1A","E2A","E0B","E1B","E2B",
  ##                   "D0A","D1A","D2A","D0B","D1B","D2B")
  states[1:nb_tip, ] <- 0
  ## I repeat the process of state assignment as many times as columns I have
  for (iv in seq_len(ncol(traits))) {
    states <- create_states(traits[, iv],
                            states,
                            sampling_fraction,
                            num_concealed_states,
                            d,
                            traitStates,
                            is_complete_tree,
                            phy,
                            ly,
                            mus,
                            nb_tip)
  }
  return(states)
}

build_initStates_time <- function(phy,
                                  traits,
                                  num_concealed_states,
                                  sampling_fraction,
                                  is_complete_tree = FALSE,
                                  mus = NULL) {
  states <- build_states(phy,
                         traits,
                         num_concealed_states,
                         sampling_fraction,
                         is_complete_tree,
                         mus)
  phy$node.label <- NULL
  split_times <- sort(event_times(phy), decreasing = FALSE)
  ances <- as.numeric(names(split_times))
  forTime <- matrix(NA, ncol = 3, nrow = nrow(phy$edge))
  forTime[, 1:2] <- phy$edge

  for (ab in seq_along(ances)) {
    focalTime <- ances[ab]
    desRows <- which(phy$edge[, 1] == focalTime)
    desNodes <- phy$edge[desRows, 2]
    rootward_age <- split_times[which(names(split_times) == focalTime)]
    for (desIndex in 1:2) {
      ## to Find the time which the integration will be done in
      if (any(desNodes[desIndex] == names(split_times))) {
        tipward_age <- split_times[which(names(split_times) ==
                                         desNodes[desIndex])]
        timeInterv <- c(tipward_age, rootward_age)
      } else {
        timeInterv <- c(0, rootward_age)
      }
      forTime[which(forTime[, 2] ==
                    desNodes[desIndex]), 3] <- timeInterv[2] - timeInterv[1]
    }
  }
  return(list(
    states = states,
    ances = ances,
    forTime = forTime
  ))
}


get_weight_states <- function(root_state_weight,
                              num_concealed_states,
                              mergeBranch,
                              lambdas,
                              nodeM,
                              d,
                              is_cla = FALSE) {

  if (is.numeric(root_state_weight)) {
    weight_states <- rep(root_state_weight / num_concealed_states,
                         num_concealed_states)
  } else {
    if (root_state_weight == "maddison_weights") {
      weight_states <- (mergeBranch) / sum((mergeBranch))
    }

    if (root_state_weight == "proper_weights") {
      if (is_cla) {
        lmb <- length(mergeBranch)
        numerator <- rep(NA, lmb)
        for (j in 1:lmb) {
          numerator[j] <- mergeBranch[j] / sum(lambdas[[j]] *
                                        ((1 - nodeM[1:d]) %o% (1 - nodeM[1:d])))
        }
        weight_states <- numerator / sum(numerator) # nolint
      } else {
      weight_states <- (mergeBranch /
                          (lambdas * (1 - nodeM[1:d]) ^ 2)) /
        sum((mergeBranch / (lambdas * (1 - nodeM[1:d]) ^ 2)))
      }
    }

    if (root_state_weight == "equal_weights") {
      weight_states <- rep(1 / length(mergeBranch), length(mergeBranch))
    }
  }

  return(weight_states)
}
