# Helper function to enter parameter value on their right place

Helper function to enter parameter value on their right place

## Usage

``` r
fill_in(object, params)
```

## Arguments

- object:

  lambda matrices, `q_matrix` or mu vector.

- params:

  parameters in order, where each value reflects the value of the
  parameter at that position, e.g. `c(0.3, 0.2, 0.1)` will fill out the
  value 0.3 for the parameter with rate identifier 1, 0.2 for the
  parameter with rate identifier 2 and 0.1 for the parameter with rate
  identifier 3.

## Value

lambda matrices, `q_matrix` or mu vector with the correct values in
their right place.
