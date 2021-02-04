cla_secsse_loglik_rhs <- function(t,y,parameter){
  ly <- length(y)
  d <- ly/2
  Es <- y[1:d]
  Ds <- y[(d+1):ly]
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  Q <- parameter[[3]]
  diag(Q) <- 0
  
  all_states <- cbind(Ds,Es)
  a <- cbind(all_states[,2],all_states[,1])
  b <- t(all_states)
  cross_D_E <- a%*%b

  dD <- -((unlist(lapply(lambdas,sum))) + mus + Q %*% (rep(1,d)))  * Ds +( Q %*% Ds ) + unlist(lapply(lapply(lambdas,"*",cross_D_E),sum))
  dE <- -((unlist(lapply(lambdas,sum))) + mus + Q %*% (rep(1,d)))  * Es + ( Q %*% Es ) + mus + unlist(lapply(lapply(lambdas,"*",Es%*%t(Es)),sum))
  
  return(list(c(dE,dD)))
}


cla_calThruNodes <- function(
  ances,
  states,
  loglik,
  forTime,
  parameter,
  use_fortran,
  methode,
  phy,
  func
){
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  Q <- parameter[[3]]
  nb_node <- phy$Nnode
  reltol <- 1e-12
  abstol <- 1e-16
  hmax <- NULL
  ly <- ncol(states)
  d <- ncol(states)/2
  
  #ances <- ances[1] #########################  REMOVE!!!!!!!!!!!!
  #desIndex <- 1
  
  focal <- ances
  desRows <- which(phy$edge[, 1] == focal)
  desNodes <- phy$edge[desRows, 2]
  
  nodeM <- numeric()
  nodeN <- numeric()

  for(desIndex in 1:2){
    y <- states[desNodes[desIndex],]
    timeInte <- forTime[which(forTime[,2] == desNodes[desIndex]),3]
    ##  To make the calculation in both lineages
    
    if(use_fortran == FALSE){
       nodeMN <- deSolve::ode(y = y,
                              func = cla_secsse_loglik_rhs,
                              times = c(0,timeInte),
                              parms = parameter,
                              rtol = reltol,
                              atol = abstol,
                              hmax = NULL,
                              method = methode)
    } else {
       nodeMN <- ode_FORTRAN(y = y,
                             func = func,
                             times = c(0,timeInte),
                             parms = parameter,
                             rtol = reltol,
                             atol = abstol,
                             method = methode)
    }
    if(desIndex == 1){
      nodeN <- nodeMN
    }
    if(desIndex == 2){
      nodeM <- nodeMN 
    }
   # print(nodeMN)
  }
  ## At the node
  nodeM <- as.numeric(nodeM[2,-1])
  nodeN <- as.numeric(nodeN[2,-1])
  ff <- normalize_loglik(nodeM[(1:d) + d],loglik); nodeM[(1:d) + d] <- ff$probs; loglik <- ff$loglik
  ff <- normalize_loglik(nodeN[(1:d) + d],loglik); nodeN[(1:d) + d] <- ff$probs; loglik <- ff$loglik
  # cat(rbind(nodeM[(d+1):length(nodeM)],nodeN[(d+1):length(nodeN)]),"\n")
  
  all_states <- cbind(nodeM[(d + 1):length(nodeM)],nodeN[(d + 1):length(nodeN)])
  a <- cbind(all_states[,2],all_states[,1])
  b <- t(all_states)
  cross_M_N <- a%*%b
  
  # mergeBranch <- NULL
  #for(iii in 1:d){
  #combProb <- 0.5*sum(lapply(lambdas,"*",cross_M_N)[[1]]) ## multiplication of probabilities of both branches
  #mergeBranch <- c(mergeBranch,combProb)
  mergeBranch <- 0.5 * (unlist(lapply(lapply(lambdas,"*",cross_M_N),sum)))
  #}
  ff <- normalize_loglik(mergeBranch,loglik); mergeBranch <- ff$probs; loglik <- ff$loglik
  #sumD <- sum(mergeBranch)
  #mergeBranch <- mergeBranch/sumD
  #loglik <- loglik + log(sumD)
  #cat(mergeBranch,"\n")
  newstate <- nodeM[1:d] ## extinction probabilities
  newstate <- c(newstate,mergeBranch)
  states[focal,] <- newstate
  #print(parameter); print(loglik)
  return(list(states = states,loglik = loglik,mergeBranch = mergeBranch,nodeM = nodeM))
}

#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
#' @importFrom doMC registerDoMC

cla_doParalThing <- function(take_ancesSub,
                         states,
                         loglik,
                         forTime,
                         parameter,
                         use_fortran,
                         methode,
                         phy,
                         func
){
  #cl <- makeCluster(2)
  #registerDoParallel(cl)
  
  
  
  #.packages=c("foreach"),
  ii <- NULL
  rm(ii)
  statesNEW <- foreach::foreach (ii = 1:2,
                                 .packages = c(
                                   "secsse",
                                   #"diversitree",
                                   "deSolve",
                                   "phylobase",
                                   "foreach",
                                   "doParallel",
                                   "geiger",
                                   "apTreeshape"),
                                 .export = c(
                                   "cla_secsse_loglik",
                                   "ode_FORTRAN",
                                   "cla_calThruNodes"
                                   )) %dopar% { 
                                     ancesSub <- take_ancesSub[[ii]]
                                     for(i in 1:length(ancesSub)){
                                       calcul <- 
                                         cla_calThruNodes(ancesSub[i],
                                                          states,
                                                          loglik,
                                                          forTime,
                                                          parameter,
                                                          use_fortran = use_fortran,
                                                          methode = methode,
                                                          phy = phy,
                                                          func = func)
                                       loglik <- calcul$loglik
                                       states <- calcul$states
                                     }
                                     list(states, loglik)
                                   }
  return(statesNEW)
}

#' Logikelihood calculation for the cla_SecSSE model given a set of parameters and data
#' @title Likelihood for SecSSE model
#' @param parameter list where the first is a table where lambdas across different modes of speciation are shown, the second mus and the third transition rates.
#' @param phy phylogenetic tree of class phylo, ultrametric, fully-resolved, rooted and with branch lengths.
#' @param traits vector with trait states, order of states must be the same as tree tips, for help, see vignette.
#' @param num_concealed_states number of concealed states, generally equivalent to number of examined states.
#' @param use_fortran Should the Fortran code for numerical integration be called? Default is TRUE.
#' @param methode Solver for the set of equations, default is "ode45".
#' @param cond condition on the existence of a node root: "maddison_cond","proper_cond"(default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:"maddison_weights","proper_weights"(default) or "equal_weights". It can also be specified the root state:the vector c(1,0,0) indicates state 1 was the root state.
#' @param sampling_fraction vector that states the sampling proportion per trait state. It must have as many elements as trait states.
#' @param run_parallel should the routine to run in parallel be called?
#' @param setting_calculation argument used internally to speed up calculation. It should be leave blank (default : setting_calculation = NULL)
#' @param setting_parallel argument used internally to set a parallel calculation. It should be left blank (default : setting_parallel = NULL)
#' @param see_ancestral_states should the ancestral states be shown? Deafault FALSE
#' @param loglik_penalty the size of the penalty for all parameters; default is 0 (no penalty)
#' @param is_complete_tree whether or not a tree with all its extinct species is provided
#' @param func function to be used in solving the ODE system. Currently only for testing purposes.
#' @note To run in parallel it is needed to load the following libraries when windows: apTreeshape, doparallel and foreach. When unix, it requires: apTreeshape, doparallel, foreach and doMC
#' @return The loglikelihood of the data given the parameters
#' @examples
#'rm(list=ls(all=TRUE))
#'library(secsse)
#'library(DDD)
#'library(deSolve)
#'#library(diversitree)
#'library(apTreeshape)
#'library(foreach)
#'set.seed(13)
#'phylotree <- ape::rcoal(12, tip.label = 1:12)
#'traits <- sample(c(0,1,2),ape::Ntip(phylotree),replace=TRUE)
#'num_concealed_states <- 3
#'sampling_fraction <- c(1,1,1)
#'methode <- "ode45"
#'phy <- phylotree

#'# the idparlist for a ETD model (dual state inheritance model of evolution)
#'# would be set like this:
#'idparlist <- cla_id_paramPos(traits,num_concealed_states)
#'lambd_and_modeSpe <- idparlist$lambdas
#'lambd_and_modeSpe[1,] <- c(1,1,1,2,2,2,3,3,3)
#'idparlist[[1]] <- lambd_and_modeSpe
#'idparlist[[2]][] <- 0
#'masterBlock <- matrix(4,ncol=3,nrow=3,byrow=TRUE) 
#'diag(masterBlock) <- NA
#'idparlist [[3]] <- q_doubletrans(traits,masterBlock,diff.conceal = FALSE)
#'# Now, internally, clasecsse sorts the lambda matrices, so they look like:
#'prepare_full_lambdas(traits,num_concealed_states,idparlist[[1]]) 
#'# which is a list with 9 matrices, corresponding to the 9 states (0A,1A,2A,0B,etc)
#'# if we want to calculate a single likelihood:
#'parameter <- idparlist
#'lambd_and_modeSpe <- parameter$lambdas
#'lambd_and_modeSpe[1,] <- c(0.2,0.2,0.2,0.4,0.4,0.4,0.01,0.01,0.01)
#'parameter[[1]] <- prepare_full_lambdas(traits,num_concealed_states,lambd_and_modeSpe) 
#'parameter[[2]] <- rep(0,9)
#'masterBlock <- matrix(0.07,ncol=3,nrow=3,byrow=TRUE) 
#'diag(masterBlock) <- NA
#'parameter [[3]] <- q_doubletrans(traits,masterBlock,diff.conceal = FALSE)
#'cla_secsse_loglik(parameter, phy, traits, num_concealed_states,
#'                  use_fortran = FALSE, methode = "ode45", cond = "maddison_cond",
#'                  root_state_weight = "maddison_weights", sampling_fraction,
#'                  run_parallel = FALSE, setting_calculation = NULL,
#'                  setting_parallel = NULL, see_ancestral_states = FALSE,
#'                  loglik_penalty = 0)
#'# LL = -37.8741
#' @export
cla_secsse_loglik <- function(parameter,
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
                                            "cla_secsse_runmod_ct",
                                            ifelse(use_fortran == FALSE,
                                                   cla_secsse_loglik_rhs,
                                                   "cla_secsse_runmod")
                                            )
                              ){
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  Q <- parameter[[3]]
  
  if(run_parallel == TRUE){ 
    if(is.null(setting_calculation)){
      check_input(traits,phy,sampling_fraction,root_state_weight,is_complete_tree)
      setting_calculation <- 
        build_initStates_time_bigtree(phy, traits, num_concealed_states, sampling_fraction, is_complete_tree, mus)
    }
    
    states <- setting_calculation$states
    forTime <- setting_calculation$forTime
    ancesSub1 <- setting_calculation$ancesSub1
    ancesSub2 <- setting_calculation$ancesSub2
    ancesRest <- setting_calculation$ancesRest
    
    if(num_concealed_states != round(num_concealed_states)){ # for testing 
      d <- ncol(states) / 2 
      new_states <- states[,c(1:sqrt(d),(d + 1):((d + 1) + sqrt(d) - 1))]
      new_states <- states[,c(1,2,3,10,11,12)]
      states <- new_states
    }
        
    loglik <- 0
    ly <- ncol(states)
    d <- ncol(states) / 2
    take_ancesSub <- list(ancesSub1, ancesSub2)
    
    if(is.null(setting_parallel)){
      if(.Platform$OS.type == "windows"){
        cl <- parallel::makeCluster(2)
        doParallel::registerDoParallel(cl)
        # pass libPath to workers
        # see https://stackoverflow.com/questions/6412459/how-to-specify-the-location-of-r-packages-in-foreach-packages-pkg-do
        # https://gitlab.com/CarlBrunius/MUVR/-/issues/11
        parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
        on.exit(parallel::stopCluster(cl))
      }
      if(.Platform$OS.type == "unix"){
        doMC::registerDoMC(2)
      } 
    }
    statesNEW <- cla_doParalThing(take_ancesSub,
                              states,
                              loglik,
                              forTime,
                              parameter,
                              use_fortran,
                              methode,
                              phy,
                              func
    )
    
    comingfromSub1 <- statesNEW[[1]][[1]]
    comingfromSub2 <- statesNEW[[2]][[1]]
    loglik <- statesNEW[[1]][[2]] + statesNEW[[2]][[2]]
    
    thoseCalculated <- 
      which(comingfromSub2[, ncol(comingfromSub2)] > 0 &
              comingfromSub2[, ncol(comingfromSub2)] < 1 &
              (is.na(comingfromSub2[, ncol(comingfromSub2)]) == FALSE))
    
    comingfromSub1[thoseCalculated, ] <- comingfromSub2[thoseCalculated, ]
    states <- comingfromSub1
    
    for(i in 1:length(ancesRest)){
      calcul <- 
        cla_calThruNodes(ancesRest[i], states, loglik, forTime, parameter, use_fortran = use_fortran,methode = methode, phy = phy, func = func)
      states <- calcul$states
      loglik <- calcul$loglik
      
    }
  } else {
    if(is.null(setting_calculation)){
      check_input(traits,phy,sampling_fraction,root_state_weight,is_complete_tree)
      setting_calculation <- build_initStates_time(phy,traits,num_concealed_states,sampling_fraction,is_complete_tree,mus)
    } 
    
    states <- setting_calculation$states
    forTime <- setting_calculation$forTime
    ances <- setting_calculation$ances
    
    if(num_concealed_states != round(num_concealed_states)){ # for testing 
      d <- ncol(states) / 2 
      new_states <- states[,c(1:sqrt(d),(d + 1):((d + 1) + sqrt(d) - 1))]
      new_states <- states[,c(1,2,3,10,11,12)]
      states <- new_states
    }
    
    loglik <- 0
    ly <- ncol(states)
    d <- ncol(states)/2
    
    for(i in 1:length(ances)){
      calcul <- cla_calThruNodes(ances[i],
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
      nodeN <- calcul$nodeN
    }
  }
  
  mergeBranch <- calcul$mergeBranch
  nodeM <- calcul$nodeM
  
  ## At the root
  
  mergeBranch2 <- (mergeBranch)
  if(is.numeric(root_state_weight)){
    giveWeights <- root_state_weight/num_concealed_states
    weightStates <- rep(giveWeights,num_concealed_states)
    
  } else {
    if(root_state_weight == "maddison_weights"){  
      weightStates <- (mergeBranch2)/sum((mergeBranch2))
    }
    
    if(root_state_weight == "proper_weights"){
      numerator <- NULL
      for(j in 1:length(mergeBranch2)){
        numerator <- c(numerator,
                      (mergeBranch2[j]/(sum(lambdas[[j]] * (1 - nodeM[1:d][j])^2))))
      }
      denomin <- NULL
      for(j in 1:length(mergeBranch2)){
        denomin <- c(denomin,(mergeBranch2[j]/(sum(lambdas[[j]] * (1 -nodeM[1:d][j]) ^2))))
      }
      
      weightStates <- numerator/sum(denomin)
    }
    #     if(root_state_weight == "proper_weights"){
    #   weightStates <- NULL
    #   for(j in 1:length(mergeBranch2)){
    #     weightStates <- c(weightStates,
    #                      (mergeBranch2[j]/(sum(lambdas[[j]] * (1 - nodeM[1:d][j]) ^2)))/
    #                        sum((mergeBranch2[j]/(sum(lambdas[[j]] * (1 -nodeM[1:d][j]) ^2)))))
    #   }
    # }
    
    if(root_state_weight == "equal_weights"){  
      weightStates <- rep(1/length(mergeBranch2),length(mergeBranch2))
    }
  }  
  
  if(cond == "maddison_cond"){
    preCond <- NULL
    for(j in 1:length(weightStates)){
      preCond <- c(preCond,
                  sum(weightStates[j] * lambdas[[j]] *  (1 - nodeM[1:d][j]) ^ 2)
      )
    }
    mergeBranch2 <- 
      mergeBranch2/(sum(preCond))
  }
  
  if(cond == "proper_cond"){
    preCond <- NULL
    for(j in 1:length(mergeBranch2)){
      preCond <- c(preCond,
                  sum((lambdas[[j]] *  (1 - nodeM[1:d][j]) ^ 2)))
    }
    mergeBranch2 <- 
      mergeBranch2/preCond
  }
  
  atRoot <- ((mergeBranch2) * (weightStates))
  
  wholeLike <- sum(atRoot)
  LL <- log(wholeLike) + loglik - penalty(pars = parameter,loglik_penalty = loglik_penalty)
  #print(unique(unlist(parameter[[1]]))); print(LL);  
  if(see_ancestral_states == TRUE){
    num_tips <- ape::Ntip(phy)
    ancestral_states <- states[(num_tips+1):nrow(states),]
    ancestral_states <- ancestral_states[,-(1:(ncol(ancestral_states)/2))]
    rownames(ancestral_states) <- ances
    return(list(ancestral_states=ancestral_states,LL=LL))
    } else {
    return(LL)
  }
}
