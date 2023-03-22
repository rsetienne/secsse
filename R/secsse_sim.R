#' Function to simulate a tree, conditional on observing all states.
#' @param lambdas speciation rates, in the form of a list of matrices
#' @param mus extinction rates, in the form of a vector
#' @param qs The Q matrix, for example the result of function q_doubletrans, but
#' generally in the form of a matrix.
#' @param crown_age crown age of the tree, tree will be simulated conditional
#' on non-extinction and this crown age.
#' @param pool_init_states pool of initial states at the crown, in case this is
#' different from all available states, otherwise leave at NULL
#' @param maxSpec Maximum number of species in the tree (please note that the
#' tree is not conditioned on this number, but that this is a safeguard against
#' generating extremely large trees).
#' @param conditioning can be 'obs_states', 'true_states' or 'none', the tree is
#' simulated until one is generated that contains all observed states, all
#' true states (e.g. obs x hidden states), or is always returned.
#' @param non_extinction should the tree be conditioned on non-extinction of the
#' crown lineages? Default is TRUE.
#' @param verbose provide intermediate output.
#' @param max_tries maximum number of simulations to try to obtain a tree.
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
                       pool_init_states = NULL,
                       maxSpec = 1e5,
                       conditioning = "none",
                       non_extinction = TRUE,
                       verbose = FALSE,
                       max_tries = 1e6) {

  if (is.matrix(lambdas)) {
    hidden_traits <- unique(gsub("[[:digit:]]+", "", names(mus)))
    # need to be converted
    lambdas <- prepare_full_lambdas(
                                names(mus),
                                num_concealed_states = length(hidden_traits),
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
      stop("pool of initial states needs to be characters, e.g. c('0A', '1B') etc")
    }
    
    # now we have to match the indices
    all_states <- names(mus)
    indices <- which(all_states %in% pool_init_states)
    pool_init_states <- -1 + indices
  }

  conditioning_vec <- c(-1)
  if (conditioning == "true_states") {
    conditioning_vec <- -1 + seq_along(mus)
  }
  if (conditioning == "obs_states") {
    obs_traits <- as.numeric(gsub("[^0-9.-]", "", names(mus)))
    conditioning_vec <- sort(unique(obs_traits))
  }

  res <- secsse_sim_cpp(mus,
                        lambdas,
                        qs,
                        crown_age,
                        maxSpec,
                        pool_init_states,
                        conditioning_vec,
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

  num_alive_species <- length(which(Ltable[, 4] == -1))

  accept_simulation <- FALSE

  if (non_extinction) {
    if (num_alive_species <= maxSpec &&
        Ltable[1, 4] == -1 &&
        Ltable[2, 4] == -1) {
      accept_simulation <- TRUE
    }
  } else {
    if (num_alive_species <= maxSpec) {
      accept_simulation <- TRUE
    }
  }

  if (accept_simulation) {
    speciesID     <- res$traits[seq(2, length(res$traits), by = 2)]
    initialState  <- res$initial_state
    Ltable[, 1]   <- crown_age - Ltable[, 1] # simulation starts at 0,
                                             # not at crown age

    indices       <- seq(1, length(res$traits), by = 2)
    speciesTraits <- 1 + res$traits[indices]
    if (verbose) {
    #  cat("making phylogeny\n")
    }
    phy <- DDD::L2phylo(Ltable, dropextinct = TRUE)

    true_traits <- sortingtraits(data.frame(cbind(paste0("t", abs(speciesID)),
                                             speciesTraits),
                                       row.names = NULL),
                            phy)
    
    true_traits <- names(mus)[true_traits]
    obs_traits <- as.numeric(gsub("[^0-9.-]", "", true_traits))
    
    return(list(phy = phy,
                true_traits = true_traits,
                obs_traits = obs_traits,
                initialState = initialState,
                extinct = res$tracker[2],
                overshoot = res$tracker[3],
                conditioning = res$tracker[4]))
  } else {
    warning("simulation did not meet minimal requirements")
    return(list(phy = "ds",
                traits = 0,
                extinct = res$tracker[2],
                overshoot = res$tracker[3],
                conditioning = res$tracker[4]))
  }
}
