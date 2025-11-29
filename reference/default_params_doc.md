# Default parameter documentation

This function's purpose is to list all parameter documentation to be
inherited by the relevant functions.

## Usage

``` r
default_params_doc(
  phy,
  traits,
  num_concealed_states,
  idparslist,
  initparsopt,
  idparsfix,
  idparsopt,
  idfactorsopt,
  parsfix,
  cond,
  root_state_weight,
  sampling_fraction,
  tol,
  maxiter,
  optimethod,
  num_cycles,
  loglik_penalty,
  is_complete_tree,
  verbose,
  num_threads,
  atol,
  rtol,
  method,
  parameter,
  setting_calculation,
  num_steps,
  see_ancestral_states,
  lambdas,
  mus,
  qs,
  crown_age,
  init_state_probs,
  maxSpec,
  conditioning,
  non_extinction,
  max_tries,
  drop_extinct,
  seed,
  prob_func,
  parameters,
  masterBlock,
  diff.conceal,
  trait_info,
  lambd_and_modeSpe,
  initloglik,
  initfactors,
  idparsfuncdefpar,
  functions_defining_params,
  state_names,
  transition_matrix,
  model,
  shift_matrix,
  q_matrix,
  lambda_list,
  object,
  params,
  param_posit,
  ml_pars,
  mu_vector,
  max_spec,
  min_spec,
  max_species_extant,
  tree_size_hist,
  start_at_crown,
  optimmethod,
  display_warning,
  take_into_account_root_edge,
  use_normalization,
  return_root_state
)
```

## Arguments

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

- idparslist:

  overview of parameters and their values.

- initparsopt:

  a numeric vector with the initial guess of the parameters to be
  estimated.

- idparsfix:

  a numeric vector with the ID of the fixed parameters.

- idparsopt:

  a numeric vector with the ID of parameters to be estimated.

- idfactorsopt:

  id of the factors that will be optimized. There are not fixed factors,
  so use a constant within `functions_defining_params`.

- parsfix:

  a numeric vector with the value of the fixed parameters.

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

- tol:

  A numeric vector with the maximum tolerance of the optimization
  algorithm. Default is `c(1e-04, 1e-05, 1e-05)`.

- maxiter:

  max number of iterations. Default is
  `1000 * round((1.25) ^ length(idparsopt))`.

- num_cycles:

  Number of cycles of the optimization. When set to `Inf`, the
  optimization will be repeated until the result is, within the
  tolerance, equal to the starting values, with a maximum of 10 cycles.

- loglik_penalty:

  the size of the penalty for all parameters; default is 0 (no penalty).

- is_complete_tree:

  logical specifying whether or not a tree with all its extinct species
  is provided. If set to `TRUE`, it also assumes that all *all* extinct
  lineages are present on the tree. Defaults to `FALSE`.

- verbose:

  sets verbose output; default is `TRUE`.

- num_threads:

  number of threads to be used. Default is one thread.

- atol:

  A numeric specifying the absolute tolerance of integration.

- rtol:

  A numeric specifying the relative tolerance of integration.

- method:

  ODE integration method. Choose from:
  `"odeint::runge_kutta_cash_karp54"`,
  `"odeint::runge_kutta_fehlberg78"`, `"odeint::runge_kutta_dopri5"`,
  `"odeint::bulirsch_stoer"` and `"odeint::runge_kutta4"`. Default
  method is: `"odeint::runge_kutta_cash_karp54"`.

- parameter:

  list where first vector represents lambdas, the second mus and the
  third transition rates.

- setting_calculation:

  argument used internally to speed up calculation. This should be left
  blank (default : `setting_calculation = NULL`).

- num_steps:

  number of substeps to show intermediate likelihoods along a branch.

- see_ancestral_states:

  Boolean for whether the ancestral states for each of the internal
  nodes should be output. Defaults to `FALSE`.

- lambdas:

  speciation rates, in the form of a list of matrices.

- mus:

  extinction rates, in the form of a vector.

- qs:

  The Q matrix, for example the result of function q_doubletrans, but
  generally in the form of a matrix.

- crown_age:

  crown age of the tree, tree will be simulated conditional on
  non-extinction and this crown age.

- init_state_probs:

  The user can provide a vector with probabilities of observing each
  state at the root, the root state is then drawn from this
  distribution. Alternatively, the user can provide a vector of
  characters representing the full names (including the concealed state)
  of states used to initialize the root (e.g. '0A', not '0').

- conditioning:

  can be `"obs_states"`, `"true_states"` or `"none"`, the tree is
  simulated until one is generated that contains all observed states
  (`"obs_states"`), all true states (e.g. all combinations of obs and
  hidden states), or is always returned (`"none"`). Alternatively, a
  vector with the names of required observed states can be provided,
  e.g. c("S", "N").

- non_extinction:

  boolean stating if the tree should be conditioned on non-extinction of
  the crown lineages. Defaults to `TRUE`.

- max_tries:

  maximum number of simulations to try to obtain a tree.

- drop_extinct:

  boolean stating if extinct species should be dropped from the tree.
  Defaults to `TRUE`.

- seed:

  pseudo-random number generator seed.

- prob_func:

  a function to calculate the probability of interest, see description.

- parameters:

  list where first vector represents lambdas, the second mus and the
  third transition rates.

- masterBlock:

  matrix of transitions among only examined states, `NA` in the main
  diagonal, used to build the full transition rates matrix.

- diff.conceal:

  Boolean stating if the concealed states should be different. E.g. that
  the transition rates for the concealed states are different from the
  transition rates for the examined states. Normally it should be
  `FALSE` in order to avoid having a huge number of parameters.

- trait_info:

  data frame where first column has species ids and the second one is
  the trait associated information.

- lambd_and_modeSpe:

  a matrix with the 4 models of speciation possible.

- initloglik:

  A numeric with the value of loglikehood obtained prior to
  optimisation. Only used internally.

- initfactors:

  the initial guess for a factor (it should be set to `NULL` when no
  factors).

- idparsfuncdefpar:

  id of the parameters which will be a function of optimized and/or
  fixed parameters. The order of id should match
  `functions_defining_params`.

- functions_defining_params:

  a list of functions. Each element will be a function which defines a
  parameter e.g. `id_3 <- (id_1 + id_2) / 2`. See example.

- state_names:

  vector of names of all observed states.

- transition_matrix:

  a matrix containing a description of all speciation events, where the
  first column indicates the source state, the second and third column
  indicate the two daughter states, and the fourth column gives the rate
  indicator used. E.g.: `["SA", "S", "A", 1]` for a trait state `"SA"`
  which upon speciation generates two daughter species with traits `"S"`
  and `"A"`, where the number 1 is used as indicator for optimization of
  the likelihood.

- model:

  used model, choice of `"ETD"` (Examined Traits Diversification),
  `"CTD"` (Concealed Traits Diversification) or `"CR"` (Constant Rate).

- shift_matrix:

  matrix of shifts, indicating in order:

  1.  starting state (typically the column in the transition matrix)

  2.  ending state (typically the row in the transition matrix)

  3.  associated rate indicator.

- q_matrix:

  `q_matrix` with only transitions between observed states.

- lambda_list:

  previously generated list of lambda matrices, used to infer the rate
  number to start with.

- object:

  lambda matrices, `q_matrix` or mu vector.

- params:

  parameters in order, where each value reflects the value of the
  parameter at that position, e.g. `c(0.3, 0.2, 0.1)` will fill out the
  value 0.3 for the parameter with rate identifier 1, 0.2 for the
  parameter with rate identifier 2 and 0.1 for the parameter with rate
  identifier 3.

- param_posit:

  initial parameter structure, consisting of a list with three entries:

  1.  lambda matrices

  2.  mus

  3.  Q matrix

  In each entry, integers numbers (1-n) indicate the parameter to be
  optimized.

- ml_pars:

  resulting parameter estimates as returned by for instance
  [`cla_secsse_ml()`](https://rsetienne.github.io/secsse/reference/cla_secsse_ml.md),
  having the same structure as `param_post`.

- mu_vector:

  previously defined mus - used to choose indicator number.

- max_spec:

  Maximum number of species in the tree (please note that the tree is
  not conditioned on this number, but that this is a safeguard against
  generating extremely large trees).

- min_spec:

  Minimum number of species in the tree.

- max_species_extant:

  Should the maximum number of species be counted in the reconstructed
  tree (if TRUE) or in the complete tree (if FALSE).

- tree_size_hist:

  if TRUE, returns a vector of all found tree sizes.

- start_at_crown:

  if FALSE, the simulation starts with one species instead of the two
  assumed by default by secsse (also in ML), and the resulting crown age
  will be lower than the set crown age. This allows for direct
  comparison with BiSSE and facilitates implementing speciation effects
  at the crown.

- optimmethod:

  A string with method used for optimization. Default is `"simplex"`.
  Alternative is `"subplex"`. Both are called from
  [`DDD::optimizer()`](https://rsetienne.github.io/DDD/reference/optimizer.html),
  simplex is implemented natively in DDD, while subplex is ultimately
  called from
  [`subplex::subplex()`](https://rdrr.io/pkg/subplex/man/subplex.html).

- display_warning:

  display a warning if necessary

- take_into_account_root_edge:

  if TRUE, the LL integration is continued along the root edge. This
  also affects conditioning (as now, conditioning no longer needs to
  assume a speciation event at the start of the tree)

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

Nothing
