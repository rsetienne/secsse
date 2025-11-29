# Package index

## All functions

- [`cla_id_paramPos()`](https://rsetienne.github.io/secsse/reference/cla_id_paramPos.md)
  : Parameter structure setting for cla_secsse It sets the parameters
  (speciation, extinction and transition) IDs. Needed for ML calculation
  with cladogenetic options (cla_secsse_ml)

- [`cla_secsse_loglik()`](https://rsetienne.github.io/secsse/reference/cla_secsse_loglik.md)
  : Likelihood for SecSSE model, using Rcpp Loglikelihood calculation
  for the cla_SecSSE model given a set of parameters and data using Rcpp

- [`cla_secsse_ml()`](https://rsetienne.github.io/secsse/reference/cla_secsse_ml.md)
  : Maximum likehood estimation for (SecSSE)

- [`cla_secsse_ml_func_def_pars()`](https://rsetienne.github.io/secsse/reference/cla_secsse_ml_func_def_pars.md)
  : Maximum likehood estimation for (SecSSE) with parameter as complex
  functions. Cladogenetic version

- [`create_default_lambda_transition_matrix()`](https://rsetienne.github.io/secsse/reference/create_default_lambda_transition_matrix.md)
  : Helper function to create a default lambda list

- [`create_default_shift_matrix()`](https://rsetienne.github.io/secsse/reference/create_default_shift_matrix.md)
  :

  Helper function to create a default `shift_matrix` list

- [`create_lambda_list()`](https://rsetienne.github.io/secsse/reference/create_lambda_list.md)
  : Helper function to automatically create lambda matrices, based on
  input. When choosing the CTD model, rates associated with observed
  states are now re-distributed to concealed states. This implicitly
  assumes that the number of observed and concealed states is identical.

- [`create_mu_vector()`](https://rsetienne.github.io/secsse/reference/create_mu_vector.md)
  : Generate mus vector

- [`create_q_matrix()`](https://rsetienne.github.io/secsse/reference/create_q_matrix.md)
  : Helper function to neatly setup a Q matrix, without transitions to
  concealed states (only observed transitions shown)

- [`event_times()`](https://rsetienne.github.io/secsse/reference/event_times.md)
  : Event times of a (possibly non-ultrametric) phylogenetic tree

- [`example_phy_GeoSSE`](https://rsetienne.github.io/secsse/reference/example_phy_GeoSSE.md)
  : A phylogeny with traits at the tips

- [`expand_q_matrix()`](https://rsetienne.github.io/secsse/reference/expand_q_matrix.md)
  : Function to expand an existing q_matrix to a number of concealed
  states

- [`extract_par_vals()`](https://rsetienne.github.io/secsse/reference/extract_par_vals.md)
  : Extract parameter values out of the result of a maximum likelihood
  inference run

- [`fill_in()`](https://rsetienne.github.io/secsse/reference/fill_in.md)
  : Helper function to enter parameter value on their right place

- [`id_paramPos()`](https://rsetienne.github.io/secsse/reference/id_paramPos.md)
  :

  Parameter structure setting Sets the parameters (speciation,
  extinction and transition) ids. Needed for ML calculation
  ([`secsse_ml()`](https://rsetienne.github.io/secsse/reference/secsse_ml.md)).

- [`phylo_vignette`](https://rsetienne.github.io/secsse/reference/phylo_vignette.md)
  : A phylogenetic reconstuction to run the vignette

- [`plot_idparslist()`](https://rsetienne.github.io/secsse/reference/plot_idparslist.md)
  : function to visualize the structure of the idparslist

- [`plot_state_exact()`](https://rsetienne.github.io/secsse/reference/plot_state_exact.md)
  : Plot the local probability along a tree

- [`prepare_full_lambdas()`](https://rsetienne.github.io/secsse/reference/prepare_full_lambdas.md)
  : Prepares the entire set of lambda matrices for cla_secsse. It
  provides the set of matrices containing all the speciation rates

- [`q_doubletrans()`](https://rsetienne.github.io/secsse/reference/q_doubletrans.md)
  : Basic Qmatrix Sets a Q matrix where double transitions are not
  allowed

- [`secsse_loglik()`](https://rsetienne.github.io/secsse/reference/secsse_loglik.md)
  : Likelihood for SecSSE model Loglikelihood calculation for the SecSSE
  model given a set of parameters and data

- [`secsse_loglik_eval()`](https://rsetienne.github.io/secsse/reference/secsse_loglik_eval.md)
  : Likelihood for SecSSE model Logikelihood calculation for the SecSSE
  model given a set of parameters and data, returning also the
  likelihoods along the branches

- [`secsse_ml()`](https://rsetienne.github.io/secsse/reference/secsse_ml.md)
  : Maximum likehood estimation for (SecSSE)

- [`secsse_ml_func_def_pars()`](https://rsetienne.github.io/secsse/reference/secsse_ml_func_def_pars.md)
  : Maximum likehood estimation for (SecSSE) complex functions as
  parameter

- [`secsse_sim()`](https://rsetienne.github.io/secsse/reference/secsse_sim.md)
  : Function to simulate a tree, conditional on observing all states.

- [`secsse_single_branch_loglik()`](https://rsetienne.github.io/secsse/reference/secsse_single_branch_loglik.md)
  : Likelihood for SecSSE model Loglikelihood calculation for the SecSSE
  model given a set of parameters and data, calculated for a single
  branch

- [`sortingtraits()`](https://rsetienne.github.io/secsse/reference/sortingtraits.md)
  : Data checking and trait sorting In preparation for likelihood
  calculation, it orders trait data according the tree tips

- [`traits`](https://rsetienne.github.io/secsse/reference/traits.md) : A
  table with trait info to run the vignette
