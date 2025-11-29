# Helper function to neatly setup a Q matrix, without transitions to concealed states (only observed transitions shown)

Helper function to neatly setup a Q matrix, without transitions to
concealed states (only observed transitions shown)

## Usage

``` r
create_q_matrix(
  state_names,
  num_concealed_states,
  shift_matrix,
  diff.conceal = FALSE
)
```

## Arguments

- state_names:

  vector of names of all observed states.

- num_concealed_states:

  number of concealed states, generally equivalent to the number of
  examined states in the dataset.

- shift_matrix:

  matrix of shifts, indicating in order:

  1.  starting state (typically the column in the transition matrix)

  2.  ending state (typically the row in the transition matrix)

  3.  associated rate indicator.

- diff.conceal:

  Boolean stating if the concealed states should be different. E.g. that
  the transition rates for the concealed states are different from the
  transition rates for the examined states. Normally it should be
  `FALSE` in order to avoid having a huge number of parameters.

## Value

transition matrix

## Examples

``` r
shift_matrix <- c(0, 1, 5)
shift_matrix <- rbind(shift_matrix, c(1, 0, 6))
q_matrix <- secsse::create_q_matrix(state_names = c(0, 1),
                                    num_concealed_states = 2,
                                    shift_matrix = shift_matrix,
                                    diff.conceal = TRUE)
```
