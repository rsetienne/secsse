# Function to expand an existing q_matrix to a number of concealed states

Function to expand an existing q_matrix to a number of concealed states

## Usage

``` r
expand_q_matrix(q_matrix, num_concealed_states, diff.conceal = FALSE)
```

## Arguments

- q_matrix:

  `q_matrix` with only transitions between observed states.

- num_concealed_states:

  number of concealed states, generally equivalent to the number of
  examined states in the dataset.

- diff.conceal:

  Boolean stating if the concealed states should be different. E.g. that
  the transition rates for the concealed states are different from the
  transition rates for the examined states. Normally it should be
  `FALSE` in order to avoid having a huge number of parameters.

## Value

updated q matrix

## Note

This is highly similar to
[`q_doubletrans()`](https://rsetienne.github.io/secsse/reference/q_doubletrans.md).
