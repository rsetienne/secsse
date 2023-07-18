#' Function to simulate a tree, conditional on observing all states.
#' 
#' @inheritParams default_params_doc
#'
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
#' modifying `pool_init_states`.
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
                       conditioning = "obs_states",
                       non_extinction = TRUE,
                       verbose = FALSE,
                       max_tries = 1e6,
                       drop_extinct = TRUE,
                       seed = NULL) {

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

  if (is.null(seed)) seed <- -1

  condition_vec <- c()
  if (is.vector(conditioning)) {
    condition_vec <- conditioning
    conditioning <- "custom"
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
                        max_tries,
                        seed,
                        condition_vec)

  Ltable        <- res$ltable
  
  if (sum(Ltable[, 4] == -1) < 2) {
    warning("crown lineages died out")
    return(list(phy = "ds",
                traits = 0,
                extinct = res$tracker[2],
                overshoot = res$tracker[3],
                conditioning = res$tracker[4]))
  }

  if (sum(res$tracker) >= max_tries) {
    warning("Couldn't simulate a tree in enough tries,
            try increasing max_tries")
    return(list(phy = "ds",
                traits = 0,
                extinct = res$tracker[2],
                overshoot = res$tracker[3],
                conditioning = res$tracker[4]))
  }

  

  initialState  <- res$initial_state
  Ltable[, 1]   <- crown_age - Ltable[, 1] # simulation starts at 0,
                                           # not at crown age
  notmin1 <- which(Ltable[, 4] != -1)
  Ltable[notmin1, 4] <- crown_age - c(Ltable[notmin1, 4])
  Ltable[which(Ltable[, 4] == crown_age + 1), 4] <- -1

  # indices       <- seq(1, length(res$traits), by = 2)
  speciesTraits <- 1 + Ltable[, 5]
  used_id <- abs(Ltable[, 3])

  phy <- DDD::L2phylo(Ltable, dropextinct = drop_extinct)
  
 
  if (drop_extinct) {
    to_drop <- which(Ltable[, 4] != -1)
    if (length(to_drop) > 0) {
      used_id <- used_id[-to_drop]
      speciesTraits <- speciesTraits[-to_drop]
    }
  }

  true_traits <- sortingtraits(data.frame(cbind(paste0("t", used_id),
                                             speciesTraits),
                                          row.names = NULL),
                               phy)

  true_traits <- names(mus)[true_traits]
  obs_traits <- c()
  for (i in seq_along(true_traits)) {
    obs_traits[i] <- substr(true_traits[i], 1, (nchar(-2) - 1))
  }
  
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
