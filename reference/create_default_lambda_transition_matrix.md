# Helper function to create a default lambda list

This function generates a generic lambda list, assuming no transitions
between states, e.g. a species of observed state 0 generates daughter
species with state 0 as well.

## Usage

``` r
create_default_lambda_transition_matrix(
  state_names = c("0", "1"),
  model = "ETD"
)
```

## Arguments

- state_names:

  vector of names of all observed states.

- model:

  used model, choice of `"ETD"` (Examined Traits Diversification),
  `"CTD"` (Concealed Traits Diversification) or `"CR"` (Constant Rate).

## Examples

``` r
lambda_matrix <-
     create_default_lambda_transition_matrix(state_names = c(0, 1),
                                             model = "ETD")
lambda_list <- create_lambda_list(state_names = c(0, 1),
                                  num_concealed_states = 2,
                                  transition_matrix = lambda_matrix,
                                  model = "ETD")
```
