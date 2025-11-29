# Plot the local probability along a tree

Plot the local probability along the tree, including the branches

## Usage

``` r
plot_state_exact(
  parameters,
  phy,
  traits,
  num_concealed_states,
  sampling_fraction,
  cond = "proper_cond",
  root_state_weight = "proper_weights",
  is_complete_tree = FALSE,
  method = "odeint::runge_kutta_cash_karp54",
  atol = 1e-16,
  rtol = 1e-16,
  num_steps = 100,
  prob_func = NULL,
  verbose = FALSE
)
```

## Arguments

- parameters:

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

- sampling_fraction:

  vector that states the sampling proportion per trait state. It must
  have as many elements as there are trait states. When using a
  `multiPhylo` object, sampling fraction should be list where each entry
  in the list corresponds to the sampling proportion for each tree.

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

- is_complete_tree:

  logical specifying whether or not a tree with all its extinct species
  is provided. If set to `TRUE`, it also assumes that all *all* extinct
  lineages are present on the tree. Defaults to `FALSE`.

- method:

  integration method used, available are:
  `"odeint::runge_kutta_cash_karp54"`,
  `"odeint::runge_kutta_fehlberg78"`, `"odeint::runge_kutta_dopri5"`,
  `"odeint::bulirsch_stoer"` and `"odeint::runge_kutta4"`. Default
  method is: `"odeint::runge_kutta_cash_karp54"`.

- atol:

  A numeric specifying the absolute tolerance of integration.

- rtol:

  A numeric specifying the relative tolerance of integration.

- num_steps:

  number of substeps to show intermediate likelihoods along a branch.

- prob_func:

  a function to calculate the probability of interest, see description.

- verbose:

  sets verbose output; default is `TRUE` when `optimmethod` is
  `"simplex"`. If `optimmethod` is set to `"simplex"`, then even if set
  to `FALSE`, optimizer output will be shown.

## Value

ggplot2 object

## Details

This function will evaluate the log likelihood locally along all
branches and plot the result. When `num_steps` is left to `NULL`, all
likelihood evaluations during integration are used for plotting. This
may work for not too large trees, but may become very memory heavy for
larger trees. Instead, the user can indicate a number of steps, which
causes the probabilities to be evaluated at a distinct amount of steps
along each branch (and the probabilities to be properly integrated in
between these steps). This provides an approximation, but generally
results look very similar to using the full evaluation. The function
used for `prob_func` will be highly dependent on your system. for
instance, for a 3 observed, 2 hidden states model, the probability of
state A is `prob[1] + prob[2] + prob[3]`, normalized by the row sum.
`prob_func` will be applied to each row of the 'states' matrix (you can
thus test your function on the states matrix returned when
`'see_ancestral_states = TRUE'`). Please note that the first N columns
of the states matrix are the extinction rates, and the `(N+1):2N`
columns belong to the speciation rates, where
`N = num_obs_states * num_concealed_states`. A typical `prob_func`
function will look like:

    my_prob_func <- function(x) {
      return(sum(x[5:8]) / sum(x))
    }

## Examples

``` r
set.seed(5)
phy <- ape::rphylo(n = 4, birth = 1, death = 0)
traits <- c(0, 1, 1, 0)
params <- secsse::id_paramPos(c(0, 1), 2)
params[[1]][] <- c(0.2, 0.2, 0.1, 0.1)
params[[2]][] <- 0.0
params[[3]][, ] <- 0.1
diag(params[[3]]) <- NA
#  Thus, we have for both, rates
# 0A, 1A, 0B and 1B. If we are interested in the posterior probability of
# trait 0,we have to provide a helper function that sums the probabilities of
# 0A and 0B, e.g.:
helper_function <- function(x) {
  return(sum(x[c(5, 7)]) / sum(x)) # normalized by total sum, just in case.
}

out_plot <- plot_state_exact(parameters = params,
                             phy = phy,
                             traits = traits,
                             num_concealed_states = 2,
                             sampling_fraction = c(1, 1),
                             num_steps = 10,
                             prob_func = helper_function)
#> Warning: Deduced names and order of used states to be: 0, 1
#> if this is incorrect, consider passing states as matching numeric 
#>   ordering, e.g. 1 for the first state, 2 for the second etc.
```
