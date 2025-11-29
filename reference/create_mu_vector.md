# Generate mus vector

Generate mus vector

## Usage

``` r
create_mu_vector(state_names, num_concealed_states, model = "CR", lambda_list)
```

## Arguments

- state_names:

  vector of names of all observed states.

- num_concealed_states:

  number of concealed states, generally equivalent to the number of
  examined states in the dataset.

- model:

  used model, choice of `"ETD"` (Examined Traits Diversification),
  `"CTD"` (Concealed Traits Diversification) or `"CR"` (Constant Rate).

- lambda_list:

  previously generated list of lambda matrices, used to infer the rate
  number to start with.

## Value

mu vector
