cla_secsse_loglik_rhs <- function(t, y, parameter) {
  ly <- length(y)
  d <- ly/2
  Es <- y[1:d]
  Ds <- y[(d + 1):ly]
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  Q <- parameter[[3]]
  diag(Q) <- 0
  
  all_states <- cbind(Ds, Es)
  a <- cbind(all_states[, 2], all_states[, 1])
  b <- t(all_states)
  cross_D_E <- a %*% b
  
  dD <- -((unlist(lapply(lambdas, sum))) +
            mus + Q %*% (rep(1, d))) *
    Ds + (Q %*% Ds) + unlist(lapply(lapply(lambdas, "*", cross_D_E), sum))
  dE <- -((unlist(lapply(lambdas, sum))) + mus + Q %*% (rep(1, d))) *
    Es +
    (Q %*% Es) +
    mus +
    unlist(lapply(lapply(lambdas, "*", Es %*% t(Es)), sum))
  
  return(list(c(dE, dD)))
}

cla_secsse_runmod_ct_e_R <- function(t, y, parameter) {
  ly <- length(y)
  d <- ly/2
  
  Es <- y[1:d]
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  Q <- parameter[[3]]
  diag(Q) <- 0
  
  dC <- 0
  for (i in 1:d) {
    dC[i] = mus[i] * (1 - Es[i])
    for (j in 1:d) {
      dC[i] <- dC[i] + Q[i, j] * (Es[j] - Es[i])
      for (k in 1:d) {
        if (lambdas[[i]][j, k] != 0) {
          dC[i] <- dC[i] + lambdas[[i]][j, k] * (Es[j] * Es[k] - Es[i])
        }
      }
    }
  }
  
  return(list(c(dC, rep(0, d))))
}


cla_calThruNodes <- function(ances,
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
  hmax <- NULL
  ly <- ncol(states)
  d <- ncol(states)/2
  
  focal <- ances
  desRows <- which(phy$edge[, 1] == focal)
  desNodes <- phy$edge[desRows, 2]
  
  nodeM <- numeric()
  nodeN <- numeric()
  
  for (desIndex in 1:2) {
    y <- states[desNodes[desIndex], ]
    timeInte <- forTime[which(forTime[, 2] == desNodes[ desIndex]), 3]
    ## To make the calculation in both lineages
    
    if (use_fortran == FALSE) {
      nodeMN <- deSolve::ode(y = y,
                             func = cla_secsse_loglik_rhs,
                             times = c(0, timeInte),
                             parms = parameter,
                             rtol = reltol,
                             atol = abstol,
                             hmax = NULL,
                             method = methode)
    } else {
      stop("FORTRAN was removed from this version")
    }
    if (desIndex == 1) {
      nodeN <- nodeMN
    }
    if (desIndex == 2) {
      nodeM <- nodeMN
    }
  }
  ## At the node
  nodeM <- as.numeric(nodeM[2, -1])
  nodeN <- as.numeric(nodeN[2, -1]) # nodeN = c(E, D)
  ff <- normalize_loglik(nodeM[(1:d) + d], loglik)
  nodeM[(1:d) + d] <- ff$probs
  loglik <- ff$loglik
  ff <- normalize_loglik(nodeN[(1:d) + d], loglik)
  nodeN[(1:d) + d] <- ff$probs
  loglik <- ff$loglik
  
  all_states <- cbind(nodeM[(d + 1):length(nodeM)],
                      nodeN[(d + 1):length(nodeN)])
  a <- cbind(all_states[, 2], all_states[, 1])
  b <- t(all_states)
  cross_M_N <- a %*% b
  
  # probabilities of both branches mergeBranch <- c(mergeBranch,combProb)
  mergeBranch <- 0.5 * (unlist(lapply(lapply(lambdas, "*", cross_M_N), sum)))
  # }
  ff <- normalize_loglik(mergeBranch, loglik)
  mergeBranch <- ff$probs
  loglik <- ff$loglik
  newstate <- nodeM[1:d]  ## extinction probabilities
  newstate <- c(newstate, mergeBranch)
  states[focal, ] <- newstate
  # print(parameter); print(loglik)
  return(list(states = states,
              loglik = loglik,
              mergeBranch = mergeBranch,
              nodeM = nodeM))
}