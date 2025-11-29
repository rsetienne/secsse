# Extract parameter values out of the result of a maximum likelihood inference run

Extract parameter values out of the result of a maximum likelihood
inference run

## Usage

``` r
extract_par_vals(param_posit, ml_pars)
```

## Arguments

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

## Value

Vector of parameter estimates.
