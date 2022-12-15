#' Function to simulate a tree under the cla secsse model
#' @param states states
#' @param lambdas speciation rates
#' @param mus extinction rates
#' @param timeSimul Crown age of the tree
#' @param qs The Q matrix, for example the result of function q_doubletrans.
#' @param pool_init_states I don't know what this is?
#' @param maxSpec Maximum number of species in the tree (please note that the
#' tree is not conditioned on this number, but that this is a safeguard against
#' generating extremely large trees).
#' @param the_traits_to_sim the traits to simulate
#' @return a list with three properties: phy: Phylogenies,
#' traits: the traits and thirdly: initialState, delineating the
#' initial state.
#' @export
cla_secsse_sim <- function(states,
                           lambdas,
                           mus,
                           timeSimul,
                           qs,
                           pool_init_states,
                           maxSpec,
                           the_traits_to_sim){
  if (length(lambdas) != length(mus) ||
      length(lambdas) != length(states) ||
      length(mus) != length(states) ) {
    stop("Every state must have a single rate of speciation and extinction")
  }
  
  if (nrow(qs) != length(states)) {
    stop("Incorrect number of transition rates")
  }
  
  speciesTraits <- c(1, 2) # to initialize the while loop  
  while (length(unique(speciesTraits)) != length(states)) {
    initialState <- sample(pool_init_states, 1)
    speciesTraits <- c(initialState, initialState)
    preTree <- secsse_sim(timeSimul,
                          states,
                          mus,
                          lambdas,
                          qs,
                          speciesTraits,
                          maxSpec)
    Ltable <- preTree$Ltable
    speciesTraits <- preTree$speciesTraits
    speciesID <- preTree$speciesID
  }
  
  if (length(speciesID) < maxSpec &  
      Ltable[1, 4] == 0 & 
      Ltable[2, 4] == 0 ) {
    age <- timeSimul
    Ltable[which(Ltable[, 4] == 0), 4] <- -1
    Ltable[,1] = age - c(Ltable[, 1])
    notmin1 = which(Ltable[,4] != -1)
    Ltable[notmin1,4] = age - c(Ltable[notmin1, 4])
    Ltable[which(Ltable[,4] == age + 1), 4] = -1
    phy <- DDD::L2phylo(Ltable, dropextinct = T)
    
    num_traits <- length(the_traits_to_sim)
    
    for (ii in 1:num_traits) {
      toConceal <- NULL
      for (jj in 1:num_traits) {
        toConceal <- c(toConceal,
                       which(speciesTraits == 
                               states[ii + (num_traits * (jj - 1))]))
        
      }
      speciesTraits[toConceal] <- ii
    }
    
    speciesTraits <- as.numeric(speciesTraits)
    traits <- sortingtraits(data.frame(cbind(paste0("t", speciesID), 
                                             speciesTraits)), 
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
