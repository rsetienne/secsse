cla_secsse_loglik_rhs <- function(t, y, parameter) {
  ly <- length(y)
  d <- ly / 2
  Es <- y[1:d] # nolint
  Ds <- y[(d + 1):ly] # nolint
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  Q <- parameter[[3]] # nolint
  diag(Q) <- 0

  all_states <- cbind(Ds, Es)
  a <- cbind(all_states[, 2], all_states[, 1])
  b <- t(all_states)
  cross_D_E <- a %*% b # nolint

  dD <- -((unlist(lapply(lambdas, sum))) + # nolint
            mus + Q %*% (rep(1, d))) *
    Ds + (Q %*% Ds) + unlist(lapply(lapply(lambdas, "*", cross_D_E), sum)) 
  dE <- -((unlist(lapply(lambdas, sum))) + mus + Q %*% (rep(1, d))) * # nolint
    Es +
    (Q %*% Es) +
    mus +
    unlist(lapply(lapply(lambdas, "*", Es %*% t(Es)), sum))

  return(list(c(dE, dD)))
}

cla_secsse_runmod_ct_e_R <- function(t, y, parameter) { # nolint
  ly <- length(y)
  d <- ly / 2

  Es <- y[1:d] # nolint
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  Q <- parameter[[3]] # nolint
  diag(Q) <- 0

  dC <- 0 # nolint
  for (i in 1:d) {
    dC[i] <- mus[i] * (1 - Es[i])
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
