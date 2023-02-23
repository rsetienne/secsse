#' Function to simulate a tree, conditional on observing all states.
#' @param lambdas speciation rates, in the form of a list of matrices
#' @param mus extinction rates, in the form of a vector
#' @param qs The Q matrix, for example the result of function q_doubletrans, but
#' generally in the form of a matrix.
#' @param crown_age crown age of the tree, tree will be simulated conditional
#' on non-extinction and this crown age.
#' @param pool_init_states pool of initial states, in case this is different
#' from all available states, otherwise leave at NULL
#' @param maxSpec Maximum number of species in the tree (please note that the
#' tree is not conditioned on this number, but that this is a safeguard against
#' generating extremely large trees).
#' @param conditioning can be 'obs_states', 'true_states' or 'none', the tree is
#' simulated until one is generated that contains all observed states, all 
#' true states (e.g. obs x hidden states), or is always returned.
#' @return a list with three properties: phy: reconstructed phylogeny,
#' traits: the traits in order of tip label and thirdly: 
#' initialState, delineating the initial state at the root used. 
#' @export
secsse_sim <- function(lambdas,
                       mus,
                       qs,
                       crown_age,
                       pool_init_states = NULL,
                       maxSpec,
                       conditioning = 'none') {
  if (length(lambdas) != length(mus)) {
    stop("Every state must have a single rate of speciation and extinction")
  }
  
  if (nrow(qs) != length(lambdas)) {
    stop("Incorrect number of transition rates")
  }
  diag(qs) <- 0
  
  if (is.null(pool_init_states)) {
    pool_init_states <- -1 + (1:length(mus))
  }
  
  conditioning_vec <- c(-1)
  if (conditioning == "true_states") {
    conditioning_vec <- -1 + 1:length(mus)
  } 
  if (conditioning == "obs_states") {
    temp_vec <- c()
    for (i in seq_along(names(mus))) {
      local_num <- substr(names(mus)[[i]], 1, 1)
      temp_vec[i] <- as.numeric(local_num)
    }
    conditioning_vec <- sort(unique(temp_vec))
  }

  res <- secsse::secsse_sim_cpp(mus,
                                lambdas,
                                qs,
                                crown_age,
                                maxSpec,
                                pool_init_states,
                                conditioning_vec)
  
  Ltable <- res$ltable
  speciesTraits <- 1 + res$traits[seq(1, length(res$traits), by = 2)]
  speciesID     <- res$traits[seq(2, length(res$traits), by = 2)]
  initialState <- res$initial_state

  Ltable[, 1] <- crown_age - Ltable[, 1] 
  
  if (length(speciesID) <= maxSpec &&  
      Ltable[1, 4] == -1 &&
      Ltable[2, 4] == -1 ) {
  
    phy <- DDD::L2phylo(Ltable, dropextinct = TRUE)
    
    traits <- sortingtraits(data.frame(cbind(paste0("t", abs(speciesID)), 
                                               speciesTraits),
                                         row.names = NULL), 
                              phy)
    return(list(phy = phy,
                traits = traits,
                initialState = initialState))
    } else {
      warning("crown lineages died out")
      return(list(phy = "ds",
                  traits = 0))
    }
}
