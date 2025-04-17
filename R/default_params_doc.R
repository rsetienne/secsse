#' Default parameter documentation
#' 
#' This function's purpose is to list all parameter documentation to be 
#' inherited by the relevant functions.
#'
#' @param phy phylogenetic tree of class `phylo`, rooted and with
#'  branch lengths. Alternatively, multiple phylogenetic trees can be provided
#'  as the `multiPhylo` class.
#' @param traits vector with trait states for each tip in the phylogeny. The 
#'  order of the states must be the same as the tree tips. For help, see 
#'  `vignette("starting_secsse", package = "secsse")`. When providing a
#'  `multiPhylo` set of multiple phylognies, traits should be a list where 
#'  each entry in the list corresponds to the matching phylogeny on that
#'  position.
#' @param num_concealed_states number of concealed states, generally equivalent
#'  to the number of examined states in the dataset.
#' @param idparslist overview of parameters and their values.
#' @param idparsopt a numeric vector with the ID of parameters to be estimated.
#' @param idfactorsopt id of the factors that will be optimized. There are not
#'  fixed factors, so use a constant within `functions_defining_params`.
#' @param initfactors the initial guess for a factor (it should be set to `NULL`
#'  when no factors).
#' @param idparsfuncdefpar id of the parameters which will be a function of
#' optimized and/or fixed parameters. The order of id should match
#' `functions_defining_params`.
#' @param functions_defining_params a list of functions. Each element will be a
#'  function which defines a parameter e.g. `id_3 <- (id_1 + id_2) / 2`. See 
#'  example.
#' @param initparsopt a numeric vector with the initial guess of the parameters 
#'  to be estimated.
#' @param idparsfix a numeric vector with the ID of the fixed parameters.
#' @param parsfix a numeric vector with the value of the fixed parameters.
#' @param cond condition on the existence of a node root: `"maddison_cond"`,
#'  `"proper_cond"` (default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:
#'  `"maddison_weights"`, `"proper_weights"` (default) or `"equal_weights"`.
#'  It can also be specified for the root state: the vector `c(1, 0, 0)` 
#'  indicates state 1 was the root state. When
#'  using a `multiPhylo` object, root_state_weight should be list where each
#'  entry in the list corresponds to the root_state_weight for each tree.
#' @param sampling_fraction vector that states the sampling proportion per
#'  trait state. It must have as many elements as there are trait states. When
#'  using a `multiPhylo` object, sampling fraction should be list where each
#'  entry in the list corresponds to the sampling proportion for each tree.
#' @param tol A numeric vector with the maximum tolerance of the optimization 
#'  algorithm. Default is `c(1e-04, 1e-05, 1e-05)`.
#' @param maxiter max number of iterations. Default is
#'  `1000 * round((1.25) ^ length(idparsopt))`.
#' @param num_cycles Number of cycles of the optimization. When set to `Inf`, 
#'  the optimization will be repeated until the result is, within the 
#'  tolerance, equal to the starting values, with a maximum of 10 cycles. 
#' @param is_complete_tree logical specifying whether or not a tree with all its
#'  extinct species is provided. If set to `TRUE`, it also assumes that all 
#'  *all* extinct lineages are present on the tree. Defaults to `FALSE`.
#' @param verbose sets verbose output; default is `TRUE` when `optimmethod` is
#'  `"simplex"`. If `optimmethod` is set to `"simplex"`, then even if set to 
#'  `FALSE`, optimizer output will be shown.
#' @param num_threads number of threads to be used. Default is one thread.
#' @param atol A numeric specifying the absolute tolerance of integration.
#' @param rtol A numeric specifying the relative tolerance of integration.
#' @param method integration method used, available are:
#'  `"odeint::runge_kutta_cash_karp54"`, `"odeint::runge_kutta_fehlberg78"`,
#'  `"odeint::runge_kutta_dopri5"`, `"odeint::bulirsch_stoer"` and
#'  `"odeint::runge_kutta4"`. Default method is: `"odeint::bulirsch_stoer"`.
#' @param parameter list where first vector represents lambdas, the second 
#'  mus and the third transition rates.
#' @param setting_calculation argument used internally to speed up calculation.
#'  It should be left blank (default : `setting_calculation = NULL`).
#' @param loglik_penalty the size of the penalty for all parameters; default is
#'  0 (no penalty).
#' @param num_steps number of substeps to show intermediate likelihoods
#'  along a branch.
#' @param see_ancestral_states Boolean for whether the ancestral states should 
#'  be shown? Defaults to `FALSE`.
#' @param lambdas speciation rates, in the form of a list of matrices.
#' @param mus extinction rates, in the form of a vector.
#' @param qs The Q matrix, for example the result of function q_doubletrans, but
#'  generally in the form of a matrix.
#' @param crown_age crown age of the tree, tree will be simulated conditional
#'  on non-extinction and this crown age.
#' @param pool_init_states pool of initial states at the crown, in case this is
#'  different from all available states, otherwise leave at NULL
#' @param max_spec Maximum number of species in the tree (please note that the
#'  tree is not conditioned on this number, but that this is a safeguard 
#'  against generating extremely large trees).
#' @param min_spec Minimum number of species in the tree.
#' @param max_species_extant Should the maximum number of species be counted in
#' the reconstructed tree (if TRUE) or in the complete tree (if FALSE).
#' @param tree_size_hist if TRUE, returns a vector of all found tree sizes. 
#' @param conditioning can be `"obs_states"`, `"true_states"` or `"none"`, the 
#'  tree is simulated until one is generated that contains all observed states
#'  (`"obs_states"`), all true states (e.g. all combinations of obs and hidden
#'  states), or is always returned (`"none"`). Alternatively, a vector with
#'  the names of required observed states can be provided, e.g. c("S", "N").
#' @param non_extinction boolean stating if the tree should be conditioned on 
#'  non-extinction of the crown lineages. Defaults to `TRUE`.
#' @param max_tries maximum number of simulations to try to obtain a tree.
#' @param drop_extinct boolean stating if extinct species should be dropped from
#'  the tree. Defaults to `TRUE`.
#' @param seed pseudo-random number generator seed.
#' @param parameters list where first vector represents lambdas, the second mus
#' and the third transition rates.
#' @param prob_func a function to calculate the probability of interest, see
#' description.
#' @param masterBlock matrix of transitions among only examined states, `NA` in
#'  the main diagonal, used to build the full transition rates matrix.
#' @param diff.conceal Boolean stating if the concealed states should be 
#'  different. E.g. that the transition rates for the concealed 
#'  states are different from the transition rates for the examined states.
#'  Normally it should be `FALSE` in order to avoid having a huge number of 
#'  parameters.
#' @param trait_info data frame where first column has species ids and the second
#'  one is the trait associated information.
#' @param optimmethod A string with method used for optimization. Default is
#' `"simplex"`. Alternative is `"subplex"`. Both are called from 
#' [DDD::optimizer()], simplex is implemented natively in [DDD], while subplex 
#' is ultimately called from [subplex::subplex()].
#' @param lambd_and_modeSpe a matrix with the 4 models of speciation possible.
#' @param initloglik A numeric with the value of loglikehood obtained prior to
#'  optimisation. Only used internally.
#' @param state_names vector of names of all observed states.
#' @param transition_matrix a matrix containing a description of all speciation
#'  events, where the first column indicates the source state, the second and
#'  third column indicate the two daughter states, and the fourth column gives
#'  the rate indicator used. E.g.: `["SA", "S", "A", 1]` for a trait state 
#'  `"SA"` which upon speciation generates two daughter species with traits 
#'  `"S"` and `"A"`, where the number 1 is used as indicator for optimization 
#'  of the likelihood.
#' @param model used model, choice of `"ETD"` (Examined Traits Diversification),
#'  `"CTD"` (Concealed Traits Diversification) or `"CR"` (Constant Rate).
#' @param shift_matrix matrix of shifts, indicating in order:
#'  1. starting state (typically the column in the transition matrix)
#'  2. ending state (typically the row in the transition matrix)
#'  3. associated rate indicator.
#' @param q_matrix `q_matrix` with only transitions between observed states.
#' @param lambda_list previously generated list of lambda matrices,
#'  used to infer the rate number to start with.
#' @param object lambda matrices, `q_matrix` or mu vector.
#' @param start_at_crown if FALSE, the simulation starts with one species
#' instead of the two assumed by default by secsse (also in ML), and 
#' the resulting crown age will be lower than the set crown age. This allows
#' for direct comparison with BiSSE and facilitates implementing speciation
#' effects at the crown.
#' @param params parameters in order, where each value reflects the value
#'  of the parameter at that position, e.g. `c(0.3, 0.2, 0.1)` will fill out
#'  the value 0.3 for the parameter with rate identifier 1, 0.2 for the
#'  parameter with rate identifier 2 and 0.1 for the parameter with rate
#'  identifier 3.
#' @param param_posit initial parameter structure, consisting of a list with
#'  three entries:
#'    1. lambda matrices
#'    2. mus
#'    3. Q matrix
#'    
#'  In each entry, integers numbers (1-n) indicate the parameter to be 
#'  optimized.
#' @param ml_pars resulting parameter estimates as returned by for instance
#'  [cla_secsse_ml()], having the same structure as `param_post`.
#' @param mu_vector previously defined mus - used to choose indicator number.
#' @param display_warning display a warning if necessary
#' @param take_into_account_root_edge if TRUE, the LL integration is continued
#' along the root edge. This also affects conditioning (as now, conditioning
#' no longer needs to assume a speciation event at the start of the tree)
#' @param use_normalization normalize the density vector during integration,
#' more accurate but slower (default = TRUE)
#' @return Nothing
#' @keywords internal
#' @export
default_params_doc <- function(phy,
                               traits,
                               num_concealed_states,
                               idparslist,
                               initparsopt,
                               idparsfix,
                               idparsopt,
                               idfactorsopt,
                               parsfix,
                               cond,
                               root_state_weight,
                               sampling_fraction,
                               tol,
                               maxiter,
                               optimethod,
                               num_cycles,
                               loglik_penalty,
                               is_complete_tree,
                               verbose,
                               num_threads,
                               atol,
                               rtol,
                               method,
                               parameter,
                               setting_calculation,
                               num_steps,
                               see_ancestral_states,
                               lambdas,
                               mus,
                               qs,
                               crown_age,
                               pool_init_states,
                               maxSpec,
                               conditioning,
                               non_extinction,
                               max_tries,
                               drop_extinct,
                               seed,
                               prob_func,
                               parameters,
                               masterBlock,
                               diff.conceal,
                               trait_info,
                               lambd_and_modeSpe,
                               initloglik,
                               initfactors,
                               idparsfuncdefpar,
                               functions_defining_params,
                               state_names,
                               transition_matrix,
                               model,
                               shift_matrix,
                               q_matrix,
                               lambda_list,
                               object,
                               params,
                               param_posit,
                               ml_pars,
                               mu_vector,
                               max_spec,
                               min_spec,
                               max_species_extant,
                               tree_size_hist,
                               start_at_crown,
                               optimmethod,
                               display_warning,
                               take_into_account_root_edge,
                               use_normalization) {
  # Nothing
}
