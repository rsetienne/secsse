check_input <- function(traits,
                        phy,
                        sampling_fraction,
                        root_state_weight,
                        is_complete_tree){
  if(is.numeric(root_state_weight)){
    if(length(root_state_weight) != length(sort(unique(traits)))){
      stop("There need to be as many elements in root_state_weight as there are traits.")
    }
    if(length(which(root_state_weight == 1)) != 1){
      stop("The root_state_weight needs only one 1.")
    }
  } else {
    if(any(root_state_weight == "maddison_weights" |
           root_state_weight == "equal_weights" |
           root_state_weight == "proper_weights") == FALSE){
      stop("The root_state_weight must be any of maddison_weights, equal_weights, or proper_weights.")
    }
  }
  
  if (ape::is.rooted(phy) == FALSE) {
    stop("The tree needs to be rooted.")
  }
  
  if (ape::is.binary(phy) == FALSE) {
    stop("The tree needs to be fully resolved.")
  }
  if (ape::is.ultrametric(phy) == FALSE & is_complete_tree == FALSE){
    stop("The tree needs to be ultrametric.")
  }
  if(any(phy$edge.length == 0)){
    stop('The tree must have internode distancs that are all larger than 0.')
  }
  
  if (is.matrix(traits)) {
    if(length(sampling_fraction) != length(sort(unique(traits[, 1])))){
      stop("Sampling_fraction must have as many elements as the number of traits.")
    }
    
    if (all(sort(unique(as.vector(traits))) == sort(unique(traits[, 1]))) == 
       FALSE){
      stop(
        "Check your trait argument; if you have more than one column, make sure all your states are included in the first column."
      )
    }
  } else{
    if (length(sampling_fraction) != length(sort(unique(traits)))){
      stop("Sampling_fraction must have as many elements as the number of traits.")
    }
  }
  
  if(length(sort(unique(as.vector(traits)))) < 2)
  {
    stop("The trait has only one state.")
  }
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
  traitStates <- sort(unique(traits[,1]))
  
  nb_tip <- ape::Ntip(phy)
  nb_node <- phy$Nnode
  ly <- length(traitStates) * 2 * num_concealed_states
  states <- matrix(ncol = ly, nrow = nb_tip + nb_node)
  d <- ly / 2
  ## In a example of 3 states, the names of the colums would be like:
  ##
  ## colnames(states) <- c("E0A","E1A","E2A","E0B","E1B","E2B",
  ##                   "D0A","D1A","D2A","D0B","D1B","D2B")
  states[1:nb_tip,] <- 0
  ## I repeat the process of state assignment as many times as columns I have
  for (iv in 1:ncol(traits)) {
    usetraits <- traits[,iv]
    if (anyNA(usetraits)) {
      nas <- which(is.na(traits))
      for (iii in 1:length(nas)) {
        states[nas[iii],] <- c(1 - rep(sampling_fraction, num_concealed_states),
                               rep(sampling_fraction, num_concealed_states))
      }
    }
    
    for (iii in 1:length(traitStates)) { # Initial state probabilities
      StatesPresents <- d + iii
      toPlaceOnes <- StatesPresents + 
                     length(traitStates) * (0:(num_concealed_states - 1))
      tipSampling <- 1 * sampling_fraction
      states[which(usetraits == traitStates[iii]), toPlaceOnes] <- tipSampling[iii]
    }
    if (is_complete_tree) {
      extinct_species <- geiger::is.extinct(phy)
      if (!is.null(extinct_species)) {
        for (i in 1:length(extinct_species)) {
          states[which(phy$tip.label == extinct_species[i]), (d + 1):ly] <- 
            mus * states[which(phy$tip.label == extinct_species[i]), (d + 1):ly]
        }
      }
      for (iii in 1:nb_tip) {
        states[iii,1:d] <- 0
      }
    } else {
      for (iii in 1:nb_tip) {
        states[iii,1:d] <- rep(1 - sampling_fraction, num_concealed_states)
      }
    }
  }
  return(states)
}

build_initStates_time <- function(phy,
                                  traits,
                                  num_concealed_states,
                                  sampling_fraction,
                                  is_complete_tree = FALSE,
                                  mus = NULL){ 
  states <- build_states(phy, 
                         traits,
                         num_concealed_states,
                         sampling_fraction,
                         is_complete_tree,
                         mus)
  phy$node.label <- NULL
  split_times <- sort(event_times(phy), decreasing = F)
  ances <- as.numeric(names(split_times))
  forTime <- matrix(NA,ncol = 3,nrow = nrow(phy$edge))
  forTime[,1:2] <- phy$edge
  
  for (ab in 1:length(ances)) {
    focalTime <- ances[ab]
    desRows <- which(phy$edge[, 1] == focalTime)
    desNodes <- phy$edge[desRows, 2]
    rootward_age <- split_times[which(names(split_times) == focalTime)]
    for (desIndex in 1:2) {
      ## to Find the time which the integration will be done in
      if (any(desNodes[desIndex] == names(split_times))) {
        tipward_age <- split_times[which(names(split_times) == desNodes[desIndex])]
        timeInterv <- c(tipward_age, rootward_age)
      } else {
        timeInterv <- c(0, rootward_age)
      }  
      forTime[which(forTime[,2] == desNodes[desIndex]),3] <- timeInterv[2] - timeInterv[1]
    }
  }
  return(list(
    states = states,
    ances = ances,
    forTime = forTime
  ))
}