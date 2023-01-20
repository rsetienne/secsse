## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width = 6)
knitr::opts_chunk$set(fig.height = 6)
library(secsse)

## ----starting_conditions------------------------------------------------------
set.seed(5)
focal_tree <- ape::rphylo(n = 4, birth = 1, death = 0)
traits <- c(0, 1, 1, 0)

plot(focal_tree)

## ----simple likelihood--------------------------------------------------------
params <- secsse::id_paramPos(c(0, 1), 2)
params[[1]][] <- c(0.2, 0.2, 0.1, 0.1)
params[[2]][] <- 0.0
params[[3]][, ] <- 0.1
diag(params[[3]]) <- NA


ll <- secsse::secsse_loglik(parameter = params,
                             phy = focal_tree,
                             traits = traits,
                             num_concealed_states = 2,
                             see_ancestral_states = TRUE,
                             sampling_fraction = c(1, 1))
ll

## ----states-------------------------------------------------------------------
ll$states

## ----helper function----------------------------------------------------------
helper_function <- function(x) {
  return(sum(x[c(5, 7)]) / sum(x)) # normalized by total sum, just in case.
}

## ----exact--------------------------------------------------------------------
secsse::plot_state_exact(parameters = params,
                 focal_tree = focal_tree,
                 traits = traits,
                 num_concealed_states = 2,
                 sampling_fraction = c(1, 1),
                 prob_func = helper_function)

secsse::plot_state_exact(parameters = params,
                 focal_tree = focal_tree,
                 traits = traits,
                 num_concealed_states = 2,
                 sampling_fraction = c(1, 1),
                 steps = 10,
                 prob_func = helper_function)

secsse::plot_state_exact(parameters = params,
                 focal_tree = focal_tree,
                 traits = traits,
                 num_concealed_states = 2,
                 sampling_fraction = c(1, 1),
                 steps = 100,
                 prob_func = helper_function)

## ----cla secsse---------------------------------------------------------------
set.seed(13)
phylotree <- ape::rcoal(12, tip.label = 1:12)
traits <- sample(c(0, 1, 2), ape::Ntip(phylotree), replace = TRUE)
num_concealed_states <- 3
sampling_fraction <- c(1, 1, 1)
phy <- phylotree
# the idparlist for a ETD model (dual state inheritance model of evolution)
# would be set like this:
idparlist <- secsse::cla_id_paramPos(traits, num_concealed_states)
lambd_and_modeSpe <- idparlist$lambdas
lambd_and_modeSpe[1, ] <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
idparlist[[1]] <- lambd_and_modeSpe
idparlist[[2]][] <- 0
masterBlock <- matrix(4, ncol = 3, nrow = 3, byrow = TRUE)
diag(masterBlock) <- NA
idparlist[[3]] <- q_doubletrans(traits, masterBlock, diff.conceal = FALSE)
# Now, internally, clasecsse sorts the lambda matrices, so they look like
#  a list with 9 matrices, corresponding to the 9 states
# (0A,1A,2A,0B, etc)

parameter <- idparlist
lambda_and_modeSpe <- parameter$lambdas
lambda_and_modeSpe[1, ] <- c(0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 0.01, 0.01, 0.01)
parameter[[1]] <- prepare_full_lambdas(traits, num_concealed_states,
lambda_and_modeSpe)
parameter[[2]] <- rep(0, 9)
masterBlock <- matrix(0.07, ncol = 3, nrow = 3, byrow = TRUE)
diag(masterBlock) <- NA
parameter[[3]] <- q_doubletrans(traits, masterBlock, diff.conceal = FALSE)

## ----helper function cla------------------------------------------------------
helper_function <- function(x) {
  return(sum(x[c(10, 13, 16)]) / sum(x)) # normalized by total sum, just in case
}

## ----plot cla-----------------------------------------------------------------
secsse::plot_state_exact_cla(parameters = parameter,
                             focal_tree = phy,
                             traits = traits,
                             num_concealed_states = 3,
                             sampling_fraction = sampling_fraction,
                             cond = "maddison_cond",
                             root_state_weight = "maddison_weights",
                             is_complete_tree = FALSE,
                             prob_func = helper_function,
                             steps = 10)
