check_input <- function(traits,
                        phy,
                        sampling_fraction,
                        root_state_weight,
                        is_complete_tree) {
  if (is.numeric(root_state_weight)) {
    if (length(root_state_weight) != length(sort(unique(traits)))) {
      stop("There need to be as many elements in root_state_weight as there are traits.")
    }
    if (length(which(root_state_weight == 1)) != 1) {
      stop("The root_state_weight needs only one 1.")
    }
  } else {
    if (any(root_state_weight == "maddison_weights" |
            root_state_weight == "equal_weights" |
            root_state_weight == "proper_weights") == FALSE) {
      stop("The root_state_weight must be any of maddison_weights, equal_weights, or proper_weights.")
    }
  }

  if (ape::is.rooted(phy) == FALSE) {
    stop("The tree needs to be rooted.")
  }

  if (ape::is.binary(phy) == FALSE) {
    stop("The tree needs to be fully resolved.")
  }
  if (ape::is.ultrametric(phy) == FALSE & is_complete_tree == FALSE) {
    stop("The tree needs to be ultrametric.")
  }
  if (any(phy$edge.length == 0)) {
    stop('The tree must have internode distancs that are all larger than 0.')
  }

  if (is.matrix(traits)) {
    if (length(sampling_fraction) != length(sort(unique(traits[, 1])))) {
      stop("Sampling_fraction must have as many elements as the number of traits.")
    }

    if (all(sort(unique(as.vector(traits))) == sort(unique(traits[, 1]))) ==
        FALSE) {
      stop(
        "Check your trait argument; if you have more than one column, make sure all your states are included in the first column."
      )
    }
  } else {
    if (length(sampling_fraction) != length(sort(unique(traits)))) {
      stop("Sampling_fraction must have as many elements as the number of traits.")
    }
  }

  if (length(sort(unique(as.vector(traits)))) < 2) {
    stop("The trait has only one state.")
  }
}

secsse_loglik_rhs <- function(t,
                              y,
                              parameter) {
  ly <- length(y)
  d <- ly/2
  Es <- y[1:d]
  Ds <- y[(d + 1):ly]
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  Q <- parameter[[3]]
  diag(Q) <- 0
  dE <- mus - (lambdas + mus + Q %*% (rep(1, d))) * Es +
       lambdas * Es * Es +
       (Q %*% Es )
  dD <- -(lambdas + mus + Q %*% (rep(1,d))) * Ds +
         2 * lambdas * Es * Ds +
         (Q %*% Ds )
  return(list(c(dE,dD)))
}

ode_FORTRAN <- function(
  y,
  times,
  func = "secsse_runmod",
  parms,
  method,
  ...
)
{
  n_vars <- length(y)
  parms <- as.numeric(unlist(parms))
  n_pars <- length(parms)
  probs <- deSolve::ode(y = y, parms = c(n_vars + 0.), rpar = parms,
                        times = times, func = func, initfunc = "secsse_initmod",
                        ynames = c("SV"), dimens = n_pars, nout = 1,
                        dllname = "secsse", method = method, ...
  )[,1:(n_vars + 1)]
  #print(as.numeric(c(probs[1,1],probs[2,c(1,6:9)])));
  return(probs)
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
  ##                   "D0A","D1A","D2B","D0B","D1B","D2B")
  states[1:nb_tip,] <- 0
  #if(is.matrix(traits)){ ## I repeat the process of state assignment as many times as columns I have
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
      toPlaceOnes <- StatesPresents + length(traitStates) * (0:(num_concealed_states - 1))
      tipSampling <- 1 * sampling_fraction
      states[which(usetraits == traitStates[iii]),toPlaceOnes] <- tipSampling[iii]
    }
    if (is_complete_tree) {
      extinct_species <- geiger::is.extinct(phy)
      if (!is.null(extinct_species)) {
        for (i in 1:length(extinct_species)) {
          states[which(phy$tip.label == extinct_species[i]),(d + 1):ly] <- mus
        }
      }
      for (iii in 1:nb_tip) {
        states[iii, 1:d] <- 0
      }
    } else {
      for (iii in 1:nb_tip) {
        states[iii, 1:d] <- rep(1 - sampling_fraction,num_concealed_states)
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

calThruNodes <- function(ances,
                         states,
                         loglik,
                         forTime,
                         parameter,
                         use_fortran,
                         methode,
                         phy,
                         func) {
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  Q <- parameter[[3]]
  nb_node <- phy$Nnode
  reltol <- 1e-12
  abstol <- 1e-16
  ly <- ncol(states)
  d <- ncol(states) / 2
  focal <- ances
  desRows <- which(phy$edge[, 1] == focal)
  desNodes <- phy$edge[desRows, 2]
  nodeM <- numeric()
  nodeN <- numeric()

  for (desIndex in 1:2) {
    y <- states[desNodes[desIndex],]
    timeInte <- forTime[which(forTime[, 2] == desNodes[desIndex]), 3]
    ##  To do the calculation in both lineages

    nodeMN <- deSolve::ode(y = y,
                           func = func,
                           times = c(0, timeInte),
                           parms = parameter,
                           rtol = reltol,
                           atol = abstol,
                           method = methode)

    if (desIndex == 1) {
      nodeN <- nodeMN
    }
    if (desIndex == 2) {
      nodeM <- nodeMN
    }
  }
  ## At the node
  nodeM <- as.numeric(nodeM[2, -1])
  nodeN <- as.numeric(nodeN[2, -1])

  ff <- normalize_loglik(nodeM[(1:d) + d], loglik);
  nodeM[(1:d) + d] <- ff$probs;
  loglik <- ff$loglik

  ff <- normalize_loglik(nodeN[(1:d) + d], loglik);
  nodeN[(1:d) + d] <- ff$probs;
  loglik <- ff$loglik

  mergeBranch <- nodeM[(1:d) + d] * nodeN[(1:d) + d] * lambdas[(1:d)]

  ff <- normalize_loglik(mergeBranch, loglik);
  mergeBranch <- ff$probs;
  loglik <- ff$loglik

  newstate <- nodeM[1:d] ## extinction probabilities
  newstate <- c(newstate, mergeBranch)
  states[focal,] <- newstate
  return(list(states = states,
              loglik = loglik,
              mergeBranch = mergeBranch,
              nodeM = nodeM))
}

build_initStates_time_bigtree <- function(phy,
                                          traits,
                                          num_concealed_states,
                                          sampling_fraction,
                                          is_complete_tree = FALSE,
                                          mus = NULL) {
  initStates_list <- build_initStates_time(phy,
                                           traits,
                                           num_concealed_states,
                                           sampling_fraction,
                                           is_complete_tree,
                                           mus)
  states <- initStates_list$states
  ances <- initStates_list$ances
  forTime <- initStates_list$forTime

  phySplit <- phy
  phySplit$node.label <- NULL
  nspp <- length(phySplit$tip.label)

  split_times <-
    sort(event_times(phySplit), decreasing = T)
  interNodes <- as.numeric(names(split_times))

  formattedtree <- apTreeshape::as.treeshape(phySplit)
  smaller <- apTreeshape::smaller.clade.spectrum(formattedtree)
  smaller <- cbind(smaller, (nspp + 1):(nrow(phySplit$edge) + 1))

  optSplit <- NULL
  for (I in 1:length(interNodes)) {
    optSplit <-
      rbind(optSplit, smaller[which(smaller[, 3] == interNodes[I]),])
  }

  optSplit <- optSplit[-1,] # To prevent the root to be chosen
  NodeOptimSplit <-
    optSplit[which(optSplit[, 2] == max(optSplit[, 2]))[1], 3]

  phy2 <- phylobase::phylo4(phy)

  sub1 <- as.numeric(phylobase::children(phy2,NodeOptimSplit)[1])
  sub2 <- as.numeric(phylobase::children(phy2,NodeOptimSplit)[2])
  jointSubs <- NodeOptimSplit

  cat("The best split is into two subtrees with",
      length(geiger::tips(phy, sub1)),"and",length(geiger::tips(phy, sub2)),"tips \n")

  descenSub1 <- sort(phylobase::descendants(phy2, sub1, type = "all"))
  descenSub1 <-
    as.numeric(descenSub1[(length(geiger::tips(phy, sub1)) + 1):length(descenSub1)])

  descenSub2 <- sort(phylobase::descendants(phy2, sub2, type = "all"))
  descenSub2 <-
    as.numeric(descenSub2[(length(geiger::tips(phy, sub2)) + 1):length(descenSub2)])

  if (length(geiger::tips(phy, sub2)) == length(descenSub2)) {
    # the siblng node has no descendent nodes
    descenSub2 <- sub2
  }

  ancesSub1 <- NULL
  for (i in 1:length(ances)) {
    if (any(ances[i] == descenSub1)) {
      ancesSub1 <- c(ancesSub1, ances[i])
    }
  }

  ancesSub2 <- NULL
  for (i in 1:length(ances)) {
    if (any(ances[i] == descenSub2)) {
      ancesSub2 <- c(ancesSub2, ances[i])
    }
  }

  ancesRest <- NULL
  for (i in 1:length(ances)) {
    if (any(any(ances[i] == descenSub1) |
            any(ances[i] == descenSub2)) == FALSE) {
      ancesRest <- c(ancesRest, ances[i])
    }
  }

  return(
    list(
      states = states,
      ancesSub1 = ancesSub1,
      ancesSub2 = ancesSub2,
      ancesRest = ancesRest,
      forTime = forTime
    )
  )
}

#' Logikelihood calculation for the SecSSE model given a set of parameters and data
#' @title Likelihood for SecSSE model
#' @param parameter list where first vector represents lambdas, the second mus and the third transition rates.
#' @param phy phylogenetic tree of class phylo, ultrametric, fully-resolved, rooted and with branch lengths.
#' @param traits vector with trait states, order of states must be the same as tree tips, for help, see vignette.
#' @param num_concealed_states number of concealed states, generally equivalent to number of examined states.
#' @param use_fortran Should the Fortran code for numerical integration be called? Default is TRUE.
#' @param methode Solver for the set of equations, default is "ode45".
#' @param cond condition on the existence of a node root: "maddison_cond","proper_cond"(default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:"maddison_weights","proper_weights"(default) or "equal_weights". It can also be specified the root state:the vector c(1,0,0) indicates state 1 was the root state.
#' @param sampling_fraction vector that states the sampling proportion per trait state. It must have as many elements as trait states.
#' @param run_parallel should the routine to run in parallel be called?
#' @param setting_calculation argument used internally to speed up calculation. It should be left blank (default : setting_calculation = NULL)
#' @param setting_parallel argument used internally to set a parallel calculation. It should be left blank (default : setting_parallel = NULL)
#' @param see_ancestral_states should the ancestral states be shown? Deafault FALSE
#' @param loglik_penalty the size of the penalty for all parameters; default is 0 (no penalty)
#' @param is_complete_tree whether or not a tree with all its extinct species is provided
#' @param func function to be used in solving the ODE system Currently only for testing purposes.
#' @note To run in parallel it is needed to load the following libraries when windows: apTreeshape, doparallel and foreach. When unix, it requires: apTreeshape, doparallel, foreach and doMC
#' @return The loglikelihood of the data given the parameters
#' @examples
#' rm(list = ls(all = TRUE))
#' library(secsse)
#' library(DDD)
#' library(deSolve)
#' library(apTreeshape)
#' library(foreach)
#' set.seed(13)
#' phylotree <- ape::rcoal(31, tip.label = 1:31)
#' traits <- sample(c(0,1,2),ape::Ntip(phylotree),replace = TRUE)
#' num_concealed_states <- 2
#' use_fortran <- TRUE
#' cond <- "proper_cond"
#' methode <- "ode45"
#' root_state_weight <- "proper_weights"
#' sampling_fraction <- c(1,1,1)
#' run_parallel <- FALSE
#' drill <- id_paramPos(traits,num_concealed_states)
#' drill[[1]][] <- c(0.12,0.01,0.2,0.21,0.31,0.23)
#' drill[[2]][] <- 0
#' drill[[3]][,] <- 0.1
#' diag(drill[[3]]) <- NA
#' secsse_loglik(parameter = drill,phylotree,traits,num_concealed_states,
#'    use_fortran,methode,cond,root_state_weight,sampling_fraction,see_ancestral_states = FALSE)
#'
#' #[1] -113.1018
#' @export
secsse_loglik <- function(parameter,
                          phy,
                          traits,
                          num_concealed_states,
                          use_fortran = TRUE,
                          methode = "ode45",
                          cond = "proper_cond",
                          root_state_weight = "proper_weights",
                          sampling_fraction,
                          run_parallel = FALSE,
                          setting_calculation = NULL,
                          setting_parallel = NULL,
                          see_ancestral_states = FALSE,
                          loglik_penalty = 0,
                          is_complete_tree = FALSE,
                          func = ifelse(is_complete_tree,
                                        "secsse_runmod_ct",
                                        ifelse(use_fortran == FALSE,
                                               secsse_loglik_rhs,
                                               "secsse_runmod2")
                          )){
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  Q <- parameter[[3]]

  if (is.null(setting_calculation)) {
    check_input(traits,phy,sampling_fraction,root_state_weight,is_complete_tree)
    setting_calculation <- build_initStates_time(phy,traits,num_concealed_states,sampling_fraction,is_complete_tree,mus)
  }

  states <- setting_calculation$states
  forTime <- setting_calculation$forTime
  ances <- setting_calculation$ances

  if (num_concealed_states != round(num_concealed_states)) { # for testing
    d <- ncol(states) / 2
    new_states <- states[,c(1:sqrt(d),(d + 1):((d + 1) + sqrt(d) - 1))]
    new_states <- states[,c(1,2,3,10,11,12)]
    states <- new_states
  }


  ly <- ncol(states)
  d <- ncol(states) / 2

  logliks <- c()
  loglik <- 0
  for (i in 1:length(ances)) {

    calcul <- calThruNodes(ances[i],
                           states,
                           loglik,
                           forTime,
                           parameter,
                           use_fortran = use_fortran,
                           methode = methode,
                           phy = phy,
                           func = func)
    states <- calcul$states
    loglik <- calcul$loglik
    if (is.na(loglik)) {stop('NA encountered in loglik')}
    logliks[i] <- loglik
    nodeN <- calcul$nodeN
  }


  mergeBranch <- calcul$mergeBranch
  nodeM <- calcul$nodeM
  ## At the root
  mergeBranch2 <- (mergeBranch)
  if (is.numeric(root_state_weight)) {
    giveWeights <- root_state_weight/num_concealed_states
    weightStates <- rep(giveWeights,num_concealed_states)
  } else {
    if (root_state_weight == "maddison_weights") {
      weightStates <- (mergeBranch2)/sum((mergeBranch2))
    }

    if (root_state_weight == "proper_weights") {
      weightStates <- (mergeBranch2/(lambdas * (1 - nodeM[1:d]) ^ 2))/sum((mergeBranch2/(lambdas * (1 - nodeM[1:d]) ^ 2)))
    }

    if (root_state_weight == "equal_weights") {
      weightStates <- rep(1/length(mergeBranch2),length(mergeBranch2))
    }
  }

  if (cond == "maddison_cond") {
    mergeBranch2 <-
      mergeBranch2 / sum(weightStates * lambdas * (1 - nodeM[1:d]) ^ 2)
  }

  if (cond == "proper_cond") {
    mergeBranch2 <- mergeBranch2 / (lambdas * (1 - nodeM[1:d]) ^ 2)
  }

  atRoot <- ((mergeBranch2) * (weightStates))

  wholeLike <- sum(atRoot)
  LL <- log(wholeLike) + loglik - penalty(pars = parameter,loglik_penalty = loglik_penalty)

  if (see_ancestral_states == TRUE) {
    num_tips <- ape::Ntip(phy)
    ancestral_states <- states[(num_tips + 1):nrow(states),]
    ancestral_states <- ancestral_states[,-(1:(ncol(ancestral_states)/2))]
    rownames(ancestral_states) <- ances
    return(list(ancestral_states = ancestral_states,LL = LL))
  } else {
    return(LL)
  }
}
