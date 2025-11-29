# Helper function to automatically create lambda matrices, based on input. When choosing the CTD model, rates associated with observed states are now re-distributed to concealed states. This implicitly assumes that the number of observed and concealed states is identical.

Helper function to automatically create lambda matrices, based on input.
When choosing the CTD model, rates associated with observed states are
now re-distributed to concealed states. This implicitly assumes that the
number of observed and concealed states is identical.

## Usage

``` r
create_lambda_list(
  state_names = c(0, 1),
  num_concealed_states = 2,
  transition_matrix,
  model = "ETD"
)
```

## Arguments

- state_names:

  vector of names of all observed states.

- num_concealed_states:

  number of concealed states, generally equivalent to the number of
  examined states in the dataset.

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

## Examples

``` r
trans_matrix <- c(0, 0, 0, 1)
trans_matrix <- rbind(trans_matrix, c(1, 1, 1, 2))
lambda_list <- create_lambda_list(state_names = c(0, 1),
                                  num_concealed_states = 2,
                                  transition_matrix = trans_matrix,
                                  model = "ETD")
```
