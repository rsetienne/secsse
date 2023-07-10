#' Default parameter documentation
#' 
#' This function's purpose is to list all parameter documentation to be 
#' inherited by the relevant functions.
#'
#' @param phy phylogenetic tree of class `phylo`, ultrametric, rooted and with
#'   branch lengths.
#' @param traits vector with trait states for each tip in the phylogeny. The 
#'   order of the states must be the same as the tree tips. For help, see 
#'   `vignette("starting_secsse", package = "secsse")`.
#' @param num_concealed_states number of concealed states, generally equivalent
#'   to the number of examined states in the dataset.
#' @param idparslist overview of parameters and their values.
#' @param initparsopt a numeric vector with the initial guess of the parameters 
#'   to be estimated.
#' @param idparsfix a numeric vector with the ID of the fixed parameters.
#' @param parsfix a numeric vector with the value of the fixed parameters.
#' @param cond condition on the existence of a node root: `"maddison_cond"`,
#'   `"proper_cond"` (default). For details, see vignette.
#' @param sampling_fraction vector that states the sampling proportion per
#'   trait state. It must have as many elements as there are trait states.
#' @param tol A numeric vector with the maximum tolerance of the optimization 
#'   algorithm. Default is `c(1e-04, 1e-05, 1e-05)`.
#' @param maxiter max number of iterations. Default is
#'   `1000 * round((1.25) ^ length(idparsopt))`.
#' @param num_cycles Number of cycles of the optimization. When set to `Inf`, 
#'   the optimization will be repeated until the result is, within the 
#'   tolerance, equal to the starting values, with a maximum of 10 cycles. 
#' @param is_complete_tree logical specifying whether or not a tree with all its
#'   extinct species is provided. If set to `TRUE`, it also assumes that all 
#'   *all* extinct lineages are present on the tree. Defaults to `FALSE`.
#' @param verbose sets verbose output; default is `TRUE` when `optimmethod` is
#'   `"subplex"`.
#' @param num_threads number of threads. Set to -1 to use all available threads.
#'   Default is one thread.
#' @param atol A numeric specifying the absolute tolerance of integration.
#' @param rtol A numeric specifying the relative tolerance of integration.
#' @param method integration method used, available are:
#'   `"odeint::runge_kutta_cash_karp54"`, `"odeint::runge_kutta_fehlberg78"`,
#'   `"odeint::runge_kutta_dopri5"`, `"odeint::bulirsch_stoer"` and
#'   `"odeint::runge_kutta4"`. Default method is: `"odeint::bulirsch_stoer"`.
#' @param parameter list where first vector represents lambdas, the second 
#'   mus and the third transition rates.
#' @param setting_calculation argument used internally to speed up calculation.
#'   It should be left blank (default : `setting_calculation = NULL`).
#' @param loglik_penalty the size of the penalty for all parameters; default is
#'   0 (no penalty).
#' @param num_steps number of substeps to show intermediate likelihoods
#'   along a branch.
#' @param see_ancestral_states Boolean for whether the ancestral states should 
#'   be shown? Defaults to `FALSE`.
#' @param lambdas speciation rates, in the form of a list of matrices.
#' @param mus extinction rates, in the form of a vector.
#' @param qs The Q matrix, for example the result of function q_doubletrans, but
#'   generally in the form of a matrix.
#' @param crown_age crown age of the tree, tree will be simulated conditional
#'   on non-extinction and this crown age.
#' @param pool_init_states pool of initial states at the crown, in case this is
#'   different from all available states, otherwise leave at NULL
#' @param maxSpec Maximum number of species in the tree (please note that the
#'   tree is not conditioned on this number, but that this is a safeguard 
#'   against generating extremely large trees).
#' @param conditioning can be `"obs_states"`, `"true_states"` or `"none"`, the 
#'   tree is simulated until one is generated that contains all observed states
#'   (`"obs_states"`), all true states (e.g. all combinations of obs and hidden
#'   states), or is always returned (`"none"`).
#' @param non_extinction boolean stating if the tree should be conditioned on 
#'   non-extinction of the crown lineages. Defaults to `TRUE`.
#' @param max_tries maximum number of simulations to try to obtain a tree.
#' @param drop_extinct boolean stating if extinct species should be dropped from
#'   the tree. Defaults to `TRUE`.
#' @param seed pseudo-random number generator seed.
#' @param parameters list where first vector represents lambdas, the second mus
#' and the third transition rates.
#' @param prob_func a function to calculate the probability of interest, see
#' description.
#' @param masterBlock matrix of transitions among only examined states, `NA` in
#'   the main diagonal, used to build the full transition rates matrix.
#' @param diff.conceal Boolean stating if the concealed states should be 
#'   different. Normally it should be `FALSE`. E.g. that the transition rates 
#'   for the concealed states are different from the transition rates for the 
#'   examined states.
#' @param trait_info data frame where first column has species ids and the second
#'   one is the trait associated information.
#' @param optimmethod A string with method used for optimization. Default is 
#'   `"subplex"`. Alternative is `"simplex"` and it shouldn't be used in normal 
#'   conditions (only for debugging). Both are called from [DDD:optimizer()], 
#'   simplex is implemented natively in [DDD], while subplex is ultimately
#'   called from [subplex::subplex()].
#' @param lambd_and_modeSpe a matrix with the 4 models of speciation possible.
#'
#' @return Nothing
#' @keywords internal
#' @export
default_params_doc <- function(phy,
                               traits,
                               num_concealed_states,
                               idparslist,
                               initparsopt,
                               idparsfix,
                               parsfix,
                               cond,
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
                               lambd_and_modeSpe) {
  # Nothing
}
