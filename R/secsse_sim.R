#' Function to simulate a tree, conditional on observing all states.
#' @param lambdas speciation rates, in the form of a list of matrices
#' @param mus extinction rates, in the form of a vector
#' @param qs The Q matrix, for example the result of function q_doubletrans, but
#' generally in the form of a matrix.
#' @param num_concealed_states number of concealed states
#' @param crown_age crown age of the tree, tree will be simulated conditional
#' on non-extinction and this crown age.
#' @param pool_init_states pool of initial states at the crown, in case this is
#' different from all available states, otherwise leave at NULL
#' @param maxSpec Maximum number of species in the tree (please note that the
#' tree is not conditioned on this number, but that this is a safeguard against
#' generating extremely large trees).
#' @param conditioning can be 'obs_states', 'true_states' or 'none', the tree is
#' simulated until one is generated that contains all observed states 
#' ('obs_states'), all true states (e.g. all combinations of obs and hidden
#' states), or is always returned ('none').
#' @param non_extinction should the tree be conditioned on non-extinction of the
#' crown lineages? Default is TRUE.
#' @param verbose provide intermediate output.
#' @param max_tries maximum number of simulations to try to obtain a tree.
#' @param drop_extinct should extinct species be dropped from the tree? default
#' is TRUE.
#' @return a list with four properties: phy: reconstructed phylogeny,
#' true_traits: the true traits in order of tip label, obs_traits: observed
#' traits, ignoring hidden traits and lastly:
#' initialState, delineating the initial state at the root used.
#' @description By default, secsse_sim assumes CLA-secsse simulation, e.g.
#' inheritance of traits at speciation need not be symmetrical, and can be
#' specified through usage of lambda-matrices. Hence, the input for lambdas
#' is typically a list of matrices.
#'
#' Simulation is performed with a randomly
#' sampled initial trait at the crown - if you, however - want a specific,
#' single, trait used at the crown, you can reduce the possible traits by
#' modifying 'pool_init_states'.
#'
#' By default, the algorithm keeps simulating until it generates a tree where
#' both crown lineages survive to the present - this is to ensure that the tree
#' has a crown age that matches the used crown age. You can modify
#' 'non-extinction' to deviate from this behaviour.
#' @export
secsse_sim <- function(lambdas,
                       mus,
                       qs,
                       crown_age,
                       num_concealed_states,
                       pool_init_states = NULL,
                       maxSpec = 1e5,
                       conditioning = "none",
                       non_extinction = TRUE,
                       verbose = FALSE,
                       max_tries = 1e6,
                       drop_extinct = TRUE) {

  if (is.matrix(lambdas)) {
    # need to be converted
    lambdas <- prepare_full_lambdas(names(mus),
                                   num_concealed_states = num_concealed_states,
                                   lambdas)
  }

  if (length(lambdas) != length(mus)) {
    stop("Every state must have a single rate of speciation and extinction")
  }

  if (nrow(qs) != length(lambdas)) {
    stop("Incorrect number of transition rates")
  }
  diag(qs) <- 0

  if (is.null(pool_init_states)) {
    pool_init_states <- -1 + (seq_along(mus))
  } else {
    if (is.numeric(pool_init_states)) {
      stop("pool of initial states needs to be characters,
           e.g. c('0A', '1B') etc")
    }

    # now we have to match the indices
    all_states <- names(mus)
    indices <- which(all_states %in% pool_init_states)
    pool_init_states <- -1 + indices
  }

  if (!conditioning %in% c("none", "true_states", "obs_states")) {
    stop("unknown conditioning, please pick from
         'none', 'obs_states', 'true_states'")
  }

  res <- secsse_sim_cpp(mus,
                        lambdas,
                        qs,
                        crown_age,
                        maxSpec,
                        pool_init_states,
                        conditioning,
                        num_concealed_states,
                        non_extinction,
                        verbose,
                        max_tries)

  if (length(res$traits) < 1) {
    warning("crown lineages died out")
    return(list(phy = "ds",
                traits = 0,
                extinct = res$tracker[2],
                overshoot = res$tracker[3],
                conditioning = res$tracker[4]))
  }

  Ltable        <- res$ltable

  speciesID     <- res$traits[seq(2, length(res$traits), by = 2)]
  initialState  <- res$initial_state
  Ltable[, 1]   <- crown_age - Ltable[, 1] # simulation starts at 0,
                                           # not at crown age
  notmin1 = which(Ltable[, 4] != -1)
  Ltable[notmin1, 4] = crown_age - c(Ltable[notmin1, 4])
  Ltable[which(Ltable[, 4] == crown_age + 1), 4] = -1          

  indices       <- seq(1, length(res$traits), by = 2)
  speciesTraits <- 1 + res$traits[indices]

  phy <- DDD::L2phylo(Ltable, dropextinct = drop_extinct)

  true_traits <- sortingtraits(data.frame(cbind(paste0("t", abs(speciesID)),
                                             speciesTraits),
                                       row.names = NULL),
                            phy)

  true_traits <- names(mus)[true_traits]
  obs_traits <- c()
  for (i in 1:length(true_traits)) {
    obs_traits[i] <- stringr::str_sub(true_traits[i], 1, -2)
  }
  #obs_traits <- as.numeric(gsub("[^0-9.-]", "", true_traits))

  if (sum(Ltable[, 4] < 0)) {
      return(list(phy = phy,
                true_traits = true_traits,
                obs_traits = obs_traits,
                initialState = initialState,
                extinct = res$tracker[2],
                overshoot = res$tracker[3],
                conditioning = res$tracker[4],
                event_counter = res$event_counter,
                extinct_draw = res$extinct_draw))
  } else {
    warning("simulation did not meet minimal requirements")
    return(list(phy = "ds",
                traits = 0,
                extinct = res$tracker[2],
                overshoot = res$tracker[3],
                conditioning = res$tracker[4],
                event_counter = res$event_counter,
                extinct_draw = res$extinct_draw))
  }
}
