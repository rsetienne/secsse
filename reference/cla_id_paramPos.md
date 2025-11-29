# Parameter structure setting for cla_secsse It sets the parameters (speciation, extinction and transition) IDs. Needed for ML calculation with cladogenetic options (cla_secsse_ml)

Parameter structure setting for cla_secsse It sets the parameters
(speciation, extinction and transition) IDs. Needed for ML calculation
with cladogenetic options (cla_secsse_ml)

## Usage

``` r
cla_id_paramPos(traits, num_concealed_states)
```

## Arguments

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

## Value

A list that includes the ids of the parameters for ML analysis.

## Examples

``` r
traits <- sample(c(0,1,2), 45,replace = TRUE) #get some traits
num_concealed_states <- 3
param_posit <- cla_id_paramPos(traits, num_concealed_states)
```
