#' @keywords internal
get_rates <- function(speciesID,
                      states,
                      speciesTraits,
                      mus,
                      lambdas,
                      qs) {
  species_mus <- NULL
  for (i in 1:length(speciesID)) {
    species_mus <- c(species_mus,
                     mus[which(states == speciesTraits[i])])
  }
  
  species_lambdas <- NULL
  for (i in 1:length(speciesID)) {
    species_lambdas <- c(species_lambdas,
                         sum(lambdas[[which(states == speciesTraits[i])]]))
    
  }
  
  shiftprob <- NULL
  for (i in 1:length(speciesID)) {
    shiftprob <- c(shiftprob,
                   sum(qs[which(speciesTraits[i] == states), ], na.rm = T))
  }
  
  return(list(species_mus = species_mus,
              species_lambdas = species_lambdas,
              shiftprob = shiftprob))
}

#' @keywords internal
event_extinction <- function(species_mus,
                             speciesTraits,
                             states,
                             Ltable,
                             speciesID,
                             mus,
                             timeStep) {
  dying <- sample(speciesID, 1, prob = species_mus) 
  Ltable[which(Ltable[,3] == dying), 4] <- timeStep
  speciesTraits <- speciesTraits[-which(speciesID == dying)]
  speciesID <- speciesID[-which(speciesID == dying)]
  return(list(speciesID = speciesID,
              speciesTraits = speciesTraits,
              Ltable = Ltable))
}


#' @keywords internal
event_speciation <- function(species_lambdas,
                             states,
                             lambdas,
                             Ltable,
                             speciesTraits,
                             speciesID,
                             timeStep) {
  if (length(speciesID) == 1) {
    mother <- speciesID
  } else {
    mother <- sample(speciesID, 1, prob = species_lambdas)  
  }
  
  mother_trait <- speciesTraits[which(mother == speciesID)]
  mother_lambdas <- lambdas[[which(mother_trait == states)]]
  
  matri_positio <- matrix(1:(length(states)^2),
                          ncol = length(states),
                          nrow = length(states),
                          byrow = FALSE)
  all_lambdas_mother <- as.vector(mother_lambdas)
  
  picked_speciation <- sample(1:(length(states)^2),
                              1, 
                              prob = all_lambdas_mother)
  
  picked_speciation_cell <- which(matri_positio == picked_speciation, 
                                  arr.ind = TRUE)
  
  state_to_parent <- picked_speciation_cell[2] # parents are in colums
  state_to_parent <- states[state_to_parent]
  state_to_daugther <- picked_speciation_cell[1]
  state_to_daugther <- states[state_to_daugther]
  speciesTraits[which(mother == speciesID)] <- state_to_parent
  speciesTraits <- c(speciesTraits,state_to_daugther)
  Ltable <- rbind(Ltable,
                  matrix(c(timeStep, mother, max(Ltable[, 3]) + 1, 0), 
                         ncol = 4))
  speciesID <- c(speciesID, max(Ltable[, 3]))
  return(list(Ltable = Ltable,
              speciesTraits = speciesTraits,
              speciesID = speciesID))
  
}

#' @keywords internal
event_traitshift <- function(shiftprob,
                             speciesTraits,
                             states,
                             speciesID,
                             qs) {
  if (length(speciesID) == 1) {
    ID_chosenspecies <- speciesID
  } else {
    ID_chosenspecies <- sample(speciesID, 1, prob = shiftprob) 
  }
  trait_chosen_species <- speciesTraits[which(speciesID == ID_chosenspecies)]
  state_prob_to_shift <- qs[which(trait_chosen_species == states), ]
  state_prob_to_shift[which(is.na(state_prob_to_shift))] <- 0
  shift_to <- sample(states, 1, prob = state_prob_to_shift )
  speciesTraits[which(speciesID == ID_chosenspecies)] <- shift_to
  return(speciesTraits = speciesTraits)
}


#' function to simulate a tree under the SecSSE model.
#' @param timeSimul Crown age of the tree
#' @param states states
#' @param lambdas speciation rates
#' @param mus extinction rates
#' @param qs transition rates in the form of a Q matrix (see the output of 
#' q_doubletrans)
#' @param speciesTraits initial species' traits.
#' @param maxSpec Maximum number of species in the tree (please note that the
#' tree is not conditioned on this number, but that this is a safeguard against
#' generating extremely large trees).
#' @return a list with four components: 1) phy: the resulting phylogeny,
#' 2) ltable the associated Ltable, 3) speciesID a list of species
#' labels and 4) speciesTraits a list of species' traits.
secsse_sim <- function(timeSimul,
                       states,
                       lambdas,
                       mus,
                       qs,
                       speciesTraits,
                       maxSpec) {
  
  timeStep <- 0
  Ltable <- matrix(c(0, 0, 1, 0, 
                     0, 1, 2, 0),
                   ncol = 4,
                   nrow = 2,
                   byrow = TRUE)
  speciesID <- c(1, 2)
  
  while (length(speciesID != 0)  && 
         timeStep <= timeSimul) {
    
    rates <- get_rates(speciesID, states, speciesTraits, mus, lambdas, qs)
    species_mus <- rates$species_mus
    species_lambdas <- rates$species_lambdas
    shiftprob <- rates$shiftprob
    totalRate <- sum(species_mus, species_lambdas, shiftprob)
    timeElapsed <- rexp(1, rate = totalRate)
    timeStep <- timeStep + timeElapsed
    
    #cat("the simulated time is:", timeStep, "\n")
    if (timeStep >= timeSimul) {
      break  
    }
    
    if (length(speciesID) >= maxSpec) {
      stop("too many species")
    }
    ##Events
    ## 1=trait shift
    ## 2=speciation
    ## 3=extinction
    
    event <- sample(c(1, 2, 3), 1,
                    prob = c(sum(shiftprob),
                             sum(species_lambdas),
                             sum(species_mus)))
    
    if (length(speciesID) == 1 && event == 3) {
      warning("Clade Extinction")
      break
    }
    
    
    if (event == 1) {
      speciesTraits <- event_traitshift(shiftprob,
                                        speciesTraits,
                                        states,
                                        speciesID,
                                        qs)
    }
    
    if (event == 2) {
      eventSpec <- event_speciation(species_lambdas,
                                    states,
                                    lambdas,
                                    Ltable,
                                    speciesTraits,
                                    speciesID,
                                    timeStep)
      Ltable <- eventSpec$Ltable
      speciesTraits <- eventSpec$speciesTraits
      speciesID <- eventSpec$speciesID
    }
    
    if (event == 3) {
      eventExt <- event_extinction(species_mus,
                                   speciesTraits,
                                   states,
                                   Ltable,
                                   speciesID,
                                   mus,
                                   timeStep)
      
      speciesID <- eventExt$speciesID
      speciesTraits <- eventExt$speciesTraits
      Ltable <- eventExt$Ltable
    }
    
  }
  return(list(phy = DDD::L2phylo(Ltable),
              Ltable = Ltable,
              speciesID = speciesID,
              speciesTraits = speciesTraits))
}
