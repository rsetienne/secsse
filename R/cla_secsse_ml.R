#' Maximum likehood estimation under Several examined and concealed
#' States-dependent Speciation and Extinction (SecSSE) with cladogenetic option
#' @title Maximum likehood estimation for (SecSSE)
#' @param phy phylogenetic tree of class phylo, ultrametric, rooted and with
#' branch lengths.
#' @param traits a vector with trait states for each tip in the phylogeny.
#' @param num_concealed_states number of concealed states, generally equivalent
#' to the number of examined states in the dataset.
#' @param idparslist overview of parameters and their values.
#' @param idparsopt id of parameters to be estimated.
#' @param initparsopt initial guess of the parameters to be estimated.
#' @param idparsfix id of the fixed parameters.
#' @param parsfix value of the fixed parameters.
#' @param cond condition on the existence of a node root: 'maddison_cond',
#' 'proper_cond'(default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:
#' 'maddison_weights','proper_weights'(default) or 'equal_weights'.
#' It can also be specified the
#' root state:the vector c(1,0,0) indicates state 1 was the root state.
#' @param sampling_fraction vector that states the sampling proportion per
#' trait state. It must have as many elements as there are trait states.
#' @param tol maximum tolerance. Default is 'c(1e-04, 1e-05, 1e-05)'.
#' @param maxiter max number of iterations. Default is
#' '1000*round((1.25)^length(idparsopt))'.
#' @param optimmethod method used for optimization. Default is 'simplex'.
#' @param num_cycles number of cycles of the optimization (default is 1).
#' @param loglik_penalty the size of the penalty for all parameters; default is
#'  0 (no penalty)
#' @param is_complete_tree whether or not a tree with all its extinct species
#' is provided
#' @param verbose sets verbose output; default is verbose when optimmethod is
#' 'subplex'
#' @param num_threads number of threads. Set to -1 to use all available threads.
#' Default is one thread.
#' @param atol absolute tolerance of integration
#' @param rtol relative tolerance of integration
#' @param method integration method used, available are:
#' "odeint::runge_kutta_cash_karp54", "odeint::runge_kutta_fehlberg78",
#' "odeint::runge_kutta_dopri5", "odeint::bulirsch_stoer" and
#' "odeint::runge_kutta4". Default method is:"odeint::bulirsch_stoer".
#' @return Parameter estimated and maximum likelihood
#' @examples
#'# Example of how to set the arguments for a ML search.
#'library(secsse)
#'library(DDD)
#'set.seed(13)
#'# Check the vignette for a better working exercise.
#'# lambdas for 0A and 1A and 2A are the same but need to be estimated
#'# (CTD model, see Syst Biol paper)
#'# mus are fixed to zero,
#'# the transition rates are constrained to be equal and fixed 0.01
#'phylotree <- ape::rcoal(31, tip.label = 1:31)
#'#get some traits
#'traits <-  sample(c(0,1,2), ape::Ntip(phylotree), replace = TRUE)
#'num_concealed_states <- 3
#'idparslist <- cla_id_paramPos(traits,num_concealed_states)
#'idparslist$lambdas[1,] <- c(1,1,1,2,2,2,3,3,3)
#'idparslist[[2]][] <- 4
#'masterBlock <- matrix(5,ncol = 3,nrow = 3,byrow = TRUE)
#'diag(masterBlock) <- NA
#'diff.conceal <- FALSE
#'idparslist[[3]] <- q_doubletrans(traits,masterBlock,diff.conceal)
#'startingpoint <- bd_ML(brts = ape::branching.times(phylotree))
#'intGuessLamba <- startingpoint$lambda0
#'intGuessMu <- startingpoint$mu0
#'idparsopt <- c(1,2,3)
#'initparsopt <- c(rep(intGuessLamba,3))
#'idparsfix <- c(0,4,5)
#'parsfix <- c(0,0,0.01)
#'tol <- c(1e-04, 1e-05, 1e-07)
#'maxiter <- 1000 * round((1.25) ^ length(idparsopt))
#'optimmethod <- 'simplex'
#'cond <- 'proper_cond'
#'root_state_weight <- 'proper_weights'
#'sampling_fraction <- c(1,1,1)
#' \dontrun{
#'model <- cla_secsse_ml(
#'  phylotree,
#'  traits,
#'  num_concealed_states,
#'  idparslist,
#'  idparsopt,
#'  initparsopt,
#'  idparsfix,
#'  parsfix,
#'  cond,
#'  root_state_weight,
#'  sampling_fraction,
#'  tol,
#'  maxiter,
#'  optimmethod,
#'  num_cycles = 1) }
#' # [1] -90.97626
#' @export
cla_secsse_ml <- function(phy,
                          traits,
                          num_concealed_states,
                          idparslist,
                          idparsopt,
                          initparsopt,
                          idparsfix,
                          parsfix,
                          cond = "proper_cond",
                          root_state_weight = "proper_weights",
                          sampling_fraction,
                          tol = c(1e-04, 1e-05, 1e-07),
                          maxiter = 1000 * round((1.25)^length(idparsopt)),
                          optimmethod = "simplex",
                          num_cycles = 1,
                          loglik_penalty = 0,
                          is_complete_tree = FALSE,
                          verbose = (optimmethod == "subplex"),
                          num_threads = 1,
                          atol = 1e-12,
                          rtol = 1e-12,
                          method = "odeint::bulirsch_stoer") {

    structure_func <- NULL
    if (is.matrix(traits)) {
        cat("you are setting a model where some species have more
            than one trait state \n")
    }

    if (length(initparsopt) != length(idparsopt)) {
        stop("initparsopt must be the same length as idparsopt.
             Number of parameters to optimize does not match the number of
             initial values for the search")
    }

    if (length(idparsfix) != length(parsfix)) {
        stop("idparsfix and parsfix must be the same length.
             Number of fixed elements does not match the fixed figures")
    }

    if (anyDuplicated(c(idparsopt, idparsfix)) != 0) {
        stop("at least one element was asked to be both fixed and estimated ")
    }

    if (identical(as.numeric(sort(c(idparsopt, idparsfix))),
                  as.numeric(sort(unique(unlist(idparslist))))) == FALSE) {
        stop("All elements in idparslist must be included in either
             idparsopt or idparsfix ")
    }

    if (anyDuplicated(c(unique(sort(as.vector(idparslist[[3]]))),
                        idparsfix[which(parsfix == 0)])) != 0) {
        cat("Note: you set some transitions as impossible to happen.", "\n")
    }

    if (is.matrix(idparslist[[1]])) {
        ## it is a tailor case otherwise
        idparslist[[1]] <- prepare_full_lambdas(traits,
                                                num_concealed_states,
                                                idparslist[[1]])
    }

    see_ancestral_states <- FALSE

    cat("Calculating the likelihood for the initial parameters ...", "\n")
    utils::flush.console()
    trparsopt <- initparsopt / (1 + initparsopt)
    trparsopt[which(initparsopt == Inf)] <- 1
    trparsfix <- parsfix / (1 + parsfix)
    trparsfix[which(parsfix == Inf)] <- 1
    mus <- calc_mus(is_complete_tree,
                    idparslist,
                    idparsfix,
                    parsfix,
                    idparsopt,
                    initparsopt)
    optimpars <- c(tol, maxiter)

    setting_calculation <- build_initStates_time(phy,
                                                 traits,
                                                 num_concealed_states,
                                                 sampling_fraction,
                                                 is_complete_tree,
                                                 mus)
    initloglik <- secsse_loglik_choosepar(trparsopt = trparsopt,
                                          trparsfix = trparsfix,
                                          idparsopt = idparsopt,
                                          idparsfix = idparsfix,
                                          idparslist = idparslist,
                                          structure_func = structure_func,
                                          phy = phy,
                                          traits = traits,
                                          num_concealed_states =
                                              num_concealed_states,
                                          cond = cond,
                                          root_state_weight = root_state_weight,
                                          sampling_fraction = sampling_fraction,
                                          setting_calculation =
                                              setting_calculation,
                                          see_ancestral_states =
                                              see_ancestral_states,
                                          loglik_penalty = loglik_penalty,
                                          is_complete_tree = is_complete_tree,
                                          verbose = verbose,
                                          num_threads = num_threads,
                                          atol = atol,
                                          rtol = rtol,
                                          method = method)
    cat("The loglikelihood for the initial parameter values is",
        initloglik, "\n")
    if (initloglik == -Inf) {
        stop("The initial parameter values have a likelihood that is 
             equal to 0 or below machine precision. 
             Try again with different initial values.")
    } else {
        cat("Optimizing the likelihood - this may take a while.", "\n")
        utils::flush.console()
        out <- DDD::optimizer(optimmethod = optimmethod,
                              optimpars = optimpars,
                              fun = secsse_loglik_choosepar,
                              trparsopt = trparsopt,
                              idparsopt = idparsopt,
                              trparsfix = trparsfix,
                              idparsfix = idparsfix,
                              idparslist = idparslist,
                              structure_func = structure_func,
                              phy = phy,
                              traits = traits,
                              num_concealed_states = num_concealed_states,
                              cond = cond,
                              root_state_weight = root_state_weight,
                              sampling_fraction = sampling_fraction,
                              setting_calculation = setting_calculation,
                              see_ancestral_states = see_ancestral_states,
                              num_cycles = num_cycles,
                              loglik_penalty = loglik_penalty,
                              is_complete_tree = is_complete_tree,
                              verbose = verbose,
                              num_threads = num_threads,
                              atol = atol,
                              rtol = rtol,
                              method = method)
        if (out$conv != 0) {
            stop("Optimization has not converged. 
                 Try again with different initial values.\n")
        } else {
            ml_pars1 <- secsse_transform_parameters(as.numeric(unlist(out$par)),
                                                   trparsfix,
                                                   idparsopt,
                                                   idparsfix,
                                                   idparslist,
                                                   structure_func)
            out2 <- list(MLpars = ml_pars1,
                         ML = as.numeric(unlist(out$fvalues)),
                         conv = out$conv)
        }
    }
    return(out2)
}
