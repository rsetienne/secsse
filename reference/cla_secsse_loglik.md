# Likelihood for SecSSE model, using Rcpp Loglikelihood calculation for the cla_SecSSE model given a set of parameters and data using Rcpp

Likelihood for SecSSE model, using Rcpp Loglikelihood calculation for
the cla_SecSSE model given a set of parameters and data using Rcpp

## Usage

``` r
cla_secsse_loglik(
  parameter,
  phy,
  traits,
  num_concealed_states,
  cond = "proper_cond",
  root_state_weight = "proper_weights",
  sampling_fraction,
  setting_calculation = NULL,
  see_ancestral_states = FALSE,
  loglik_penalty = 0,
  is_complete_tree = FALSE,
  take_into_account_root_edge = FALSE,
  num_threads = 1,
  method = "odeint::runge_kutta_cash_karp54",
  atol = 1e-08,
  rtol = 1e-07,
  display_warning = TRUE,
  use_normalization = TRUE,
  return_root_state = FALSE
)
```

## Arguments

- parameter:

  list where first vector represents lambdas, the second mus and the
  third transition rates.

- phy:

  phylogenetic tree of class `phylo`, rooted and with branch lengths.
  Alternatively, multiple phylogenetic trees can be provided as the
  `multiPhylo` class.

- traits:

  vector with trait states for each tip in the phylogeny. The order of
  the states must be the same as the tree tips. For help, see
  [`vignette("starting_secsse", package = "secsse")`](https://rsetienne.github.io/secsse/articles/starting_secsse.md).
  When providing a `multiPhylo` set of multiple phylognies, traits
  should be a list where each entry in the list corresponds to the
  matching phylogeny on that position.

- num_concealed_states:

  number of concealed states, generally equivalent to the number of
  examined states in the dataset.

- cond:

  condition on the existence of a node root: `"maddison_cond"`,
  `"proper_cond"` (default). For details, see vignette.

- root_state_weight:

  the method to weigh the states: `"maddison_weights"`,
  `"proper_weights"` (default) or `"equal_weights"`. It can also be
  specified for the root state: the vector `c(1, 0, 0)` indicates state
  1 was the root state. When using a `multiPhylo` object,
  root_state_weight should be list where each entry in the list
  corresponds to the root_state_weight for each tree.

- sampling_fraction:

  vector that states the sampling proportion per trait state. It must
  have as many elements as there are trait states. When using a
  `multiPhylo` object, sampling fraction should be list where each entry
  in the list corresponds to the sampling proportion for each tree.

- setting_calculation:

  argument used internally to speed up calculation. This should be left
  blank (default : `setting_calculation = NULL`).

- see_ancestral_states:

  Boolean for whether the ancestral states for each of the internal
  nodes should be output. Defaults to `FALSE`.

- loglik_penalty:

  the size of the penalty for all parameters; default is 0 (no penalty).

- is_complete_tree:

  logical specifying whether or not a tree with all its extinct species
  is provided. If set to `TRUE`, it also assumes that all *all* extinct
  lineages are present on the tree. Defaults to `FALSE`.

- take_into_account_root_edge:

  if TRUE, the LL integration is continued along the root edge. This
  also affects conditioning (as now, conditioning no longer needs to
  assume a speciation event at the start of the tree)

- num_threads:

  number of threads to be used. Default is one thread.

- method:

  ODE integration method. Choose from:
  `"odeint::runge_kutta_cash_karp54"`,
  `"odeint::runge_kutta_fehlberg78"`, `"odeint::runge_kutta_dopri5"`,
  `"odeint::bulirsch_stoer"` and `"odeint::runge_kutta4"`. Default
  method is: `"odeint::runge_kutta_cash_karp54"`.

- atol:

  A numeric specifying the absolute tolerance of integration.

- rtol:

  A numeric specifying the relative tolerance of integration.

- display_warning:

  display a warning if necessary

- use_normalization:

  normalize the density vector during integration, more accurate but
  slower (default = TRUE)

- return_root_state:

  if TRUE, returns the state of the system at the root, this can be
  useful to use as the starting point of a simulation. When used in ML,
  after finishing the ML optimization, the found optimum is evaluated
  one more time to retrieve the root state (to avoid having to store the
  root state every ML evaluation).

## Value

A List with property LL: The loglikelihood of the data given the
parameters, and potentially the root state.

## Examples

``` r
rm(list=ls(all=TRUE))
library(secsse)
set.seed(13)
phylotree <- ape::rcoal(12, tip.label = 1:12)
traits <- sample(c(0,1,2),ape::Ntip(phylotree),replace=TRUE)
num_concealed_states <- 3
sampling_fraction <- c(1,1,1)
phy <- phylotree
# the idparlist for a ETD model (dual state inheritance model of evolution)
# would be set like this:
idparlist <- cla_id_paramPos(traits,num_concealed_states)
lambd_and_modeSpe <- idparlist$lambdas
lambd_and_modeSpe[1,] <- c(1,1,1,2,2,2,3,3,3)
idparlist[[1]] <- lambd_and_modeSpe
idparlist[[2]][] <- 0
masterBlock <- matrix(4,ncol=3,nrow=3,byrow=TRUE)
diag(masterBlock) <- NA
idparlist [[3]] <- q_doubletrans(traits,masterBlock,diff.conceal = FALSE)
# Now, internally, clasecsse sorts the lambda matrices, so they look like:
prepare_full_lambdas(traits,num_concealed_states,idparlist[[1]])
#> [[1]]
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    1    0    0    0    0    0    0    0    0
#>  [2,]    0    0    0    0    0    0    0    0    0
#>  [3,]    0    0    0    0    0    0    0    0    0
#>  [4,]    0    0    0    0    0    0    0    0    0
#>  [5,]    0    0    0    0    0    0    0    0    0
#>  [6,]    0    0    0    0    0    0    0    0    0
#>  [7,]    0    0    0    0    0    0    0    0    0
#>  [8,]    0    0    0    0    0    0    0    0    0
#>  [9,]    0    0    0    0    0    0    0    0    0
#> 
#> [[2]]
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    0    0    0    0    0    0    0    0    0
#>  [2,]    0    1    0    0    0    0    0    0    0
#>  [3,]    0    0    0    0    0    0    0    0    0
#>  [4,]    0    0    0    0    0    0    0    0    0
#>  [5,]    0    0    0    0    0    0    0    0    0
#>  [6,]    0    0    0    0    0    0    0    0    0
#>  [7,]    0    0    0    0    0    0    0    0    0
#>  [8,]    0    0    0    0    0    0    0    0    0
#>  [9,]    0    0    0    0    0    0    0    0    0
#> 
#> [[3]]
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    0    0    0    0    0    0    0    0    0
#>  [2,]    0    0    0    0    0    0    0    0    0
#>  [3,]    0    0    1    0    0    0    0    0    0
#>  [4,]    0    0    0    0    0    0    0    0    0
#>  [5,]    0    0    0    0    0    0    0    0    0
#>  [6,]    0    0    0    0    0    0    0    0    0
#>  [7,]    0    0    0    0    0    0    0    0    0
#>  [8,]    0    0    0    0    0    0    0    0    0
#>  [9,]    0    0    0    0    0    0    0    0    0
#> 
#> [[4]]
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    0    0    0    0    0    0    0    0    0
#>  [2,]    0    0    0    0    0    0    0    0    0
#>  [3,]    0    0    0    0    0    0    0    0    0
#>  [4,]    0    0    0    2    0    0    0    0    0
#>  [5,]    0    0    0    0    0    0    0    0    0
#>  [6,]    0    0    0    0    0    0    0    0    0
#>  [7,]    0    0    0    0    0    0    0    0    0
#>  [8,]    0    0    0    0    0    0    0    0    0
#>  [9,]    0    0    0    0    0    0    0    0    0
#> 
#> [[5]]
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    0    0    0    0    0    0    0    0    0
#>  [2,]    0    0    0    0    0    0    0    0    0
#>  [3,]    0    0    0    0    0    0    0    0    0
#>  [4,]    0    0    0    0    0    0    0    0    0
#>  [5,]    0    0    0    0    2    0    0    0    0
#>  [6,]    0    0    0    0    0    0    0    0    0
#>  [7,]    0    0    0    0    0    0    0    0    0
#>  [8,]    0    0    0    0    0    0    0    0    0
#>  [9,]    0    0    0    0    0    0    0    0    0
#> 
#> [[6]]
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    0    0    0    0    0    0    0    0    0
#>  [2,]    0    0    0    0    0    0    0    0    0
#>  [3,]    0    0    0    0    0    0    0    0    0
#>  [4,]    0    0    0    0    0    0    0    0    0
#>  [5,]    0    0    0    0    0    0    0    0    0
#>  [6,]    0    0    0    0    0    2    0    0    0
#>  [7,]    0    0    0    0    0    0    0    0    0
#>  [8,]    0    0    0    0    0    0    0    0    0
#>  [9,]    0    0    0    0    0    0    0    0    0
#> 
#> [[7]]
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    0    0    0    0    0    0    0    0    0
#>  [2,]    0    0    0    0    0    0    0    0    0
#>  [3,]    0    0    0    0    0    0    0    0    0
#>  [4,]    0    0    0    0    0    0    0    0    0
#>  [5,]    0    0    0    0    0    0    0    0    0
#>  [6,]    0    0    0    0    0    0    0    0    0
#>  [7,]    0    0    0    0    0    0    3    0    0
#>  [8,]    0    0    0    0    0    0    0    0    0
#>  [9,]    0    0    0    0    0    0    0    0    0
#> 
#> [[8]]
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    0    0    0    0    0    0    0    0    0
#>  [2,]    0    0    0    0    0    0    0    0    0
#>  [3,]    0    0    0    0    0    0    0    0    0
#>  [4,]    0    0    0    0    0    0    0    0    0
#>  [5,]    0    0    0    0    0    0    0    0    0
#>  [6,]    0    0    0    0    0    0    0    0    0
#>  [7,]    0    0    0    0    0    0    0    0    0
#>  [8,]    0    0    0    0    0    0    0    3    0
#>  [9,]    0    0    0    0    0    0    0    0    0
#> 
#> [[9]]
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    0    0    0    0    0    0    0    0    0
#>  [2,]    0    0    0    0    0    0    0    0    0
#>  [3,]    0    0    0    0    0    0    0    0    0
#>  [4,]    0    0    0    0    0    0    0    0    0
#>  [5,]    0    0    0    0    0    0    0    0    0
#>  [6,]    0    0    0    0    0    0    0    0    0
#>  [7,]    0    0    0    0    0    0    0    0    0
#>  [8,]    0    0    0    0    0    0    0    0    0
#>  [9,]    0    0    0    0    0    0    0    0    3
#> 
# which is a list with 9 matrices, corresponding to the 9 states
# (0A,1A,2A,0B,etc)
# if we want to calculate a single likelihood:
parameter <- idparlist
lambda_and_modeSpe <- parameter$lambdas
lambda_and_modeSpe[1,] <- c(0.2,0.2,0.2,0.4,0.4,0.4,0.01,0.01,0.01)
parameter[[1]] <- prepare_full_lambdas(traits,num_concealed_states,
lambda_and_modeSpe)
parameter[[2]] <- rep(0,9)
masterBlock <- matrix(0.07, ncol=3, nrow=3, byrow=TRUE)
diag(masterBlock) <- NA
parameter [[3]] <- q_doubletrans(traits,masterBlock,diff.conceal = FALSE)
cla_secsse_loglik(parameter, phy, traits, num_concealed_states,
                 cond = 'maddison_cond',
                 root_state_weight = 'maddison_weights', sampling_fraction,
                 setting_calculation = NULL,
                 see_ancestral_states = FALSE,
                 loglik_penalty = 0)
#> Warning: Deduced names and order of used states to be: 0, 1, 2
#> if this is incorrect, consider passing states as matching numeric 
#>   ordering, e.g. 1 for the first state, 2 for the second etc.
#> [1] -42.18407
# LL = -42.18407
```
