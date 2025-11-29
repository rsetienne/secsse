# Function to simulate a tree, conditional on observing all states.

By default, secsse_sim assumes CLA-secsse simulation, e.g. inheritance
of traits at speciation need not be symmetrical, and can be specified
through usage of lambda-matrices. Hence, the input for lambdas is
typically a list of matrices.

Simulation is performed with a randomly sampled initial trait at the
crown - if you, however - want a specific, single, trait used at the
crown, you can reduce the possible traits by modifying
`init_state_probs`.

By default, the algorithm keeps simulating until it generates a tree
where both crown lineages survive to the present - this is to ensure
that the tree has a crown age that matches the used crown age. You can
modify 'non-extinction' to deviate from this behaviour.

## Usage

``` r
secsse_sim(
  lambdas,
  mus,
  qs,
  crown_age,
  num_concealed_states,
  init_state_probs = NULL,
  sampling_fraction = NULL,
  max_spec = 1e+05,
  min_spec = 2,
  max_species_extant = TRUE,
  tree_size_hist = FALSE,
  conditioning = "none",
  non_extinction = TRUE,
  verbose = FALSE,
  max_tries = 1e+06,
  drop_extinct = TRUE,
  start_at_crown = TRUE,
  seed = NULL
)
```

## Arguments

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

- num_concealed_states:

  number of concealed states, generally equivalent to the number of
  examined states in the dataset.

- init_state_probs:

  The user can provide a vector with probabilities of observing each
  state at the root, the root state is then drawn from this
  distribution. Alternatively, the user can provide a vector of
  characters representing the full names (including the concealed state)
  of states used to initialize the root (e.g. '0A', not '0').

- sampling_fraction:

  vector that states the sampling proportion per trait state. It must
  have as many elements as there are trait states. When using a
  `multiPhylo` object, sampling fraction should be list where each entry
  in the list corresponds to the sampling proportion for each tree.

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

- verbose:

  sets verbose output; default is `TRUE`.

- max_tries:

  maximum number of simulations to try to obtain a tree.

- drop_extinct:

  boolean stating if extinct species should be dropped from the tree.
  Defaults to `TRUE`.

- start_at_crown:

  if FALSE, the simulation starts with one species instead of the two
  assumed by default by secsse (also in ML), and the resulting crown age
  will be lower than the set crown age. This allows for direct
  comparison with BiSSE and facilitates implementing speciation effects
  at the crown.

- seed:

  pseudo-random number generator seed.

## Value

a list with four properties: phy: reconstructed phylogeny, true_traits:
the true traits in order of tip label, obs_traits: observed traits,
ignoring hidden traits and lastly: initialState, delineating the initial
state at the root used.
