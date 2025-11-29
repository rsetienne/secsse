# Helper function to create a default `shift_matrix` list

This function generates a generic shift matrix to be used with the
function
[`create_q_matrix()`](https://rsetienne.github.io/secsse/reference/create_q_matrix.md).

## Usage

``` r
create_default_shift_matrix(
  state_names = c("0", "1"),
  num_concealed_states = 2,
  mu_vector = NULL
)
```

## Arguments

- state_names:

  vector of names of all observed states.

- num_concealed_states:

  number of concealed states, generally equivalent to the number of
  examined states in the dataset.

- mu_vector:

  previously defined mus - used to choose indicator number.

## Examples

``` r
shift_matrix <- create_default_shift_matrix(state_names = c(0, 1),
                                            num_concealed_states = 2,
                                            mu_vector = c(1, 2, 1, 2))
q_matrix <- create_q_matrix(state_names = c(0, 1),
                            num_concealed_states = 2,
                            shift_matrix = shift_matrix,
                            diff.conceal = FALSE)
```
