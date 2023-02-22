#' @keywords internal
get_rates <- function(speciesID,
                      states,
                      speciesTraits,
                      mus,
                      lambdas,
                      qs) {
  species_mus <- rep(NA, length(speciesID))
  species_lambdas <- rep(NA, length(speciesID))
  shiftprob <- rep(NA, length(speciesID))
  for (i in seq_along(speciesID)) {
    matches <- which(states == speciesTraits[i])
    species_mus[i] <- mus[matches]
    species_lambdas[i] <- sum(lambdas[[matches]])
    shiftprob[i] <- sum(qs[matches, ], na.rm = TRUE)
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
  all_lambdas_mother <- as.vector(mother_lambdas)
  picked_speciation <- sample(1:(length(states)^2), 1, prob = all_lambdas_mother)
  matri_positio <- matrix(1:(length(states)^2),
                          ncol = length(states),
                          nrow = length(states),
                          byrow = FALSE)
  picked_speciation_cell <- which(matri_positio == picked_speciation,
                                  arr.ind = TRUE)
  
  state_to_parent <- picked_speciation_cell[2] # parents are in colums
  state_to_parent <- states[state_to_parent]
  state_to_daugther <- picked_speciation_cell[1]
  state_to_daugther <- states[state_to_daugther]
  speciesTraits[which(mother == speciesID)] <- state_to_parent
  speciesTraits <- c(speciesTraits,state_to_daugther)
  
  newL = nrow(Ltable)
  newL = newL + 1
  Ltable <- rbind(Ltable,
                  matrix(c(timeStep, mother, sign(mother) * newL, -1),
                         ncol = 4)) 
  speciesID <- c(speciesID, sign(mother) * newL)
  
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
#' @export
secsse_sim_R <- function(timeSimul,
                       states,
                       lambdas,
                       mus,
                       qs,
                       speciesTraits,
                       maxSpec) {

  timeStep <- 0
  Ltable <- matrix(c(0, 0, -1, -1,
                     0, -1, 2, -1),
                   ncol = 4,
                   nrow = 2,
                   byrow = TRUE)
  speciesID <- c(-1, 2)

  while (length(speciesID != 0) && timeStep <= timeSimul) {
    rates <- get_rates(speciesID, states, speciesTraits, mus, lambdas, qs)
    species_mus <- rates$species_mus
    species_lambdas <- rates$species_lambdas
    shiftprob <- rates$shiftprob
    totalRate <- sum(species_mus, species_lambdas, shiftprob)
    timeElapsed <- rexp(1, rate = totalRate)
    timeStep <- timeStep + timeElapsed

    if (timeStep > timeSimul) {
      break 
    }

    if (length(speciesID) > maxSpec) {
    #  warning("Too many species")
      break
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
     # print("Clade Extinction")
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
  
  #correcting Ltable to use DDD::L2phylo
  new_Ltable <- Ltable
  notmin1 <- which(Ltable[,4] != -1)
  Ltable[notmin1, 4] <- timeSimul - c(Ltable[notmin1, 4])

  time_change <- timeSimul - (Ltable[,1])
  new_Ltable <- cbind(time_change, Ltable[,2:4])

  #getting only examined traits
  examTraits <- c()
  for (i in 1:length(speciesTraits)) {
    focal_thing <- speciesTraits[i]
    local_num <- substr(focal_thing, 1, 1)
    examTraits[i] <- as.numeric(local_num)
  }

  # building phylogeny and matching order of traits and tips
  resulting_phylogeny <- NA
  if (length(speciesID) != 1) {
    resulting_phylogeny <- DDD::L2phylo(as.matrix(new_Ltable),
                                        dropextinct = TRUE)
    examTraits <- secsse::sortingtraits(
                            data.frame(cbind(paste0("t", abs(speciesID)), 
                                             examTraits), row.names = NULL),
                                             phy = resulting_phylogeny)
  }

  return(list(phy = resulting_phylogeny,
              new_Ltable = new_Ltable,
              speciesID = speciesID,
              speciesTraits = speciesTraits,
              examTraits = examTraits))
}
