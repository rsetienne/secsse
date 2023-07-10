#' Default parameter documentation
#' 
#' This function's purpose is to list all parameter documentation to be 
#' inherited by the relevant functions.
#'
#' @param phy phylogenetic tree of class `phylo`, ultrametric, rooted and with
#'   branch lengths.
#' @param traits a vector with trait states for each tip in the phylogeny.
#' @param num_concealed_states number of concealed states, generally equivalent
#'   to the number of examined states in the dataset.
#' @param idparslist overview of parameters and their values.
#' @param idparsopt a numeric vector with the ID of parameters to be estimated.
#' @param initparsopt a numeric vector with the initial guess of the parameters 
#'   to be estimated.
#' @param idparsfix a numeric vector with the ID of the fixed parameters.
#' @param parsfix a numeric vector with the value of the fixed parameters.
#' @param cond condition on the existence of a node root: `"maddison_cond"`,
#'   `"proper_cond"` (default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:
#'   `"maddison_weights"`, `"proper_weights"` (default) or `"equal_weights"`.
#'   It can also be specified for the root state: the vector `c(1, 0, 0)` 
#'   indicates state 1 was the root state.
#' @param sampling_fraction vector that states the sampling proportion per
#'   trait state. It must have as many elements as there are trait states.
#' @param tol A numeric vector with the maximum tolerance of the optimization 
#'   algorithm. Default is `c(1e-04, 1e-05, 1e-05)`.
#' @param maxiter max number of iterations. Default is
#'   `1000 * round((1.25) ^ length(idparsopt))`.
#' @param optimmethod method used for optimization. Available are simplex and
#'   subplex, default is `"subplex"`. Simplex should only be used for debugging.
#' @param num_cycles Number of cycles of the optimization. When set to `Inf`, 
#'   the optimization will be repeated until the result is, within the 
#'   tolerance, equal to the starting values, with a maximum of 10 cycles. 
#' @param is_complete_tree logical specifying whether or not a tree with all its
#'   extinct species is provided. If set to `TRUE`, it also assumes that all 
#'   *all* extinct lineages are present on the tree. Defaults to `FALSE`.
#' @param verbose sets verbose output; default is verbose when optimmethod is
#'   `'subplex'`.
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
#'
#' @return Nothing
#' @export
#' @keywords internal
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
                               num_cyles,
                               loglik_penalty,
                               is_complete_tree,
                               verbose,
                               num_threads,
                               atol,
                               rtol,
                               method,
                               parameter,
                               setting_calculation,
                               num_steps) {
  # Nothing
}
