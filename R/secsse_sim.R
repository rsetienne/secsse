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
                       sampling_fraction = NULL,
                       max_spec = 1e5,
                       min_spec = 2,
                       max_species_extant = TRUE,
                       tree_size_hist = FALSE,
                       conditioning = "none",
                       non_extinction = TRUE,
                       verbose = FALSE,
                       max_tries = 1e6,
                       drop_extinct = TRUE,
                       start_at_crown = TRUE,
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
    if (length(indices) < 1) {
      # most likely states without hidden labels have been provided
      letters_conceal <- LETTERS[1:num_concealed_states]
      all_states2 <- c()
      for (i in 1:length(pool_init_states)) {
        for (j in 1:length(letters_conceal)) {
          all_states2 <- c(all_states2, paste0(pool_init_states[i], letters_conceal[j]))
        }
      }
      indices <- which(all_states %in% all_states2)
    }
    pool_init_states <- -1 + indices
  }

  if (is.null(seed)) seed <- -1

  condition_vec <- vector()
  if (length(conditioning) > 1) {
    condition_vec <- conditioning
    conditioning <- "custom"
    true_traits <- names(mus)
    all_states <- c()
    for (i in seq_along(true_traits)) {
      all_states[i] <- substr(true_traits[i], 1, (nchar(true_traits[i]) - 1))
    }
    all_states <- unique(all_states)
    
    indices <- which(all_states %in% condition_vec)
    condition_vec <- -1 + indices
  } 
  
  res <- generate_phy(mus,
                      lambdas,
                      qs,
                      crown_age,
                      max_spec,
                      max_species_extant,
                      min_spec,
                      pool_init_states,
                      conditioning,
                      num_concealed_states,
                      non_extinction,
                      verbose,
                      max_tries,
                      seed,
                      condition_vec,
                      tree_size_hist,
                      start_at_crown,
                      drop_extinct) 
  
  if (res$status == "success" || res$status == "single_species_tree") {
    if (sum(sampling_fraction) == length(sampling_fraction) ||
        is.null(sampling_fraction)) {
      return(res) # sampling fraction = 1
    }
    
    # if not, we have to subsample
    all_traits <- names(mus)
    for (i in 1:length(all_traits)) {
      all_traits[i] <- substr(all_traits[i], 1, (nchar(all_traits[i]) - 1))
    }
    all_traits <- unique(all_traits) # these are now all observed traits
    tips_to_remove <- c()
    for (i in seq_along(sampling_fraction)) {
      trait_tips <- which(res$obs_traits == all_traits[i])
      to_keep <- stats::rbinom(length(trait_tips), 1, sampling_fraction[i]) # sampling fraction is survival probability
      trait_tips_for_removal <- trait_tips[to_keep == 0]
      tips_to_remove <- c(tips_to_remove, trait_tips_for_removal)
    }
    # now we need to remove all the tips
    if (length(tips_to_remove) > 0) {
      res$obs_traits <- res$obs_traits[-tips_to_remove]
      res$true_traits <- res$true_traits[-tips_to_remove]
      res$phy <- ape::drop.tip(res$phy, tips_to_remove)
    }
  }
  return(res)
}

#' @keywords internal
generate_phy <- function(mus,
                         lambdas,
                         qs,
                         crown_age,
                         max_spec,
                         max_species_extant,
                         min_spec,
                         pool_init_states,
                         conditioning,
                         num_concealed_states,
                         non_extinction,
                         verbose,
                         max_tries,
                         seed,
                         condition_vec,
                         tree_size_hist,
                         start_at_crown,
                         drop_extinct) {
  res <- secsse_sim_cpp(mus,
                        lambdas,
                        qs,
                        crown_age,
                        max_spec,
                        max_species_extant,
                        min_spec,
                        pool_init_states,
                        conditioning,
                        num_concealed_states,
                        non_extinction,
                        verbose,
                        max_tries,
                        seed,
                        condition_vec,
                        tree_size_hist,
                        start_at_crown)
  
  if (length(res) < 1) { # this happens upon a throw
    return(list(phy = "ds",
                traits = 0,
                status = "error"))
  }
  
  
  
  Ltable        <- res$ltable
  
  out_hist <- 0
  if (tree_size_hist == TRUE) out_hist <- res$hist_tree_size
  
  if (sum(res$tracker) >= max_tries) {
    warning("Couldn't simulate a tree in enough tries,
            try increasing max_tries")
    
    return(list(phy = "ds",
                traits = 0,
                extinct = res$tracker[2],
                overshoot = res$tracker[3],
                conditioning = res$tracker[4],
                small = res$tracker[6],
                size_hist = out_hist,
                status = "not enough tries"))
  }
  
  if (start_at_crown == FALSE && sum(Ltable[, 4] == -1) == 1) {
    # fake phy
    phy <- ape::rphylo(n = 2, birth = 0.2, death = 0)
    phy$edge.length[-2]
    phy$tip.label <- phy$tip.label[-2]
    phy$edge <- phy$edge[-2, ]
    phy$edge <- matrix(data = phy$edge, nrow = 1) # important!
    phy$edge.length <- crown_age # this is now root age
    
    speciesTraits <- 1 + Ltable[, 5]
    
    true_traits <- names(mus)[speciesTraits]
    obs_traits <- c()
    for (i in seq_along(true_traits)) {
      obs_traits[i] <- substr(true_traits[i], 1, (nchar(true_traits[i]) - 1))
    }
    initialState  <- names(mus)[1 + res$initial_state]
    return(list(phy = phy,
                true_traits = true_traits,
                obs_traits = obs_traits,
                initialState = initialState,
                extinct = res$tracker[2],
                overshoot = res$tracker[3],
                conditioning = res$tracker[4],
                small = res$tracker[6],
                size_hist = out_hist,
                status = "single_species_tree"))
    
    
  } else if (sum(Ltable[, 4] == -1) < 2) {
    warning("crown lineages died out")
    return(list(phy = "ds",
                traits = 0,
                extinct = res$tracker[2],
                overshoot = res$tracker[3],
                conditioning = res$tracker[4],
                small = res$tracker[6],
                size_hist = out_hist,
                status = "extinction"))
  }
  
  
  
  
  
  initialState  <- names(mus)[1 + res$initial_state]
  Ltable[, 1]   <- crown_age - Ltable[, 1] # simulation starts at 0,
  # not at crown age
  notmin1 <- which(Ltable[, 4] != -1)
  Ltable[notmin1, 4] <- crown_age - c(Ltable[notmin1, 4])
  Ltable[which(Ltable[, 4] == crown_age + 1), 4] <- -1
  
  # indices       <- seq(1, length(res$traits), by = 2)
  speciesTraits <- 1 + Ltable[, 5]
  used_id <- abs(Ltable[, 3])
  
  #phy <- DDD::L2phylo(Ltable, dropextinct = drop_extinct)
  phy <- treestats::l_to_phylo(Ltable, drop_extinct = drop_extinct)
  
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
    obs_traits[i] <- substr(true_traits[i], 1, (nchar(true_traits[i]) - 1))
  }
  
  if (sum(Ltable[, 4] < 0) > 0) {
    return(list(phy = phy,
                true_traits = true_traits,
                obs_traits = obs_traits,
                initialState = initialState,
                extinct = res$tracker[2],
                overshoot = res$tracker[3],
                conditioning = res$tracker[4],
                small = res$tracker[6],
                size_hist = out_hist,
                status = "success"))
  } else {
    warning("simulation did not meet minimal requirements")
    return(list(phy = "ds",
                traits = 0,
                extinct = res$tracker[2],
                overshoot = res$tracker[3],
                conditioning = res$tracker[4],
                small = res$tracker[6],
                size_hist = out_hist,
                status = "requirements not met"))
  }
}