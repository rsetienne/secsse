#' Function to simulate a tree, conditional on observing all states.
#' @param states states
#' @param lambdas speciation rates
#' @param mus extinction rates
#' @param timeSimul Crown age of the tree
#' @param qs The Q matrix, for example the result of function q_doubletrans.
#' @param pool_init_states pool of initial states, in case this is different
#' from all available states.
#' @param maxSpec Maximum number of species in the tree (please note that the
#' tree is not conditioned on this number, but that this is a safeguard against
#' generating extremely large trees).
#' @return a list with three properties: phy: Phylogenies,
#' traits: the traits and thirdly: initialState, delineating the
#' initial state.
#' @export
secsse_sim <- function(lambdas,
                       mus,
                       timeSimul,
                       qs,
                       pool_init_states = NULL,
                       maxSpec) {
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
  
  res <- secsse::secsse_sim_cpp(mus,
                                lambdas,
                                qs,
                                timeSimul,
                                maxSpec,
                                pool_init_states)
  
  Ltable <- res$ltable
  speciesTraits <- 1 + res$traits[seq(1, length(res$traits), by = 2)]
  speciesID     <- res$traits[seq(2, length(res$traits), by = 2)]
  initialState <- res$initial_state

  Ltable[, 1] <- timeSimul - Ltable[, 1] 
  
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
