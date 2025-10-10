#' Maximum likehood estimation for (SecSSE) complex functions as parameter
#' 
#' Maximum likehood estimation under Several examined and concealed
#' States-dependent Speciation and Extinction (SecSSE) where some parameters
#' are functions of other parameters and/or factors.
#' 
#' @inheritParams default_params_doc
#' 
#' @return Parameter estimates and maximum likelihood
#' @examples
#'# Example of how to set the arguments for an ML search. The ML search is stopped
#'# after 10 iterations to keep run time short.
#'library(secsse)
#'library(DDD)
#'set.seed(16)
#'phylotree <- ape::rbdtree(0.07,0.001,Tmax=50)
#'startingpoint<-bd_ML(brts = ape::branching.times(phylotree))
#'intGuessLamba <- startingpoint$lambda0
#'intGuessMu <- startingpoint$mu0
#'traits <- sample(c(0,1,2), ape::Ntip(phylotree),replace=TRUE) #get some traits
#'num_concealed_states<-3
#'idparslist<-id_paramPos(traits, num_concealed_states)
#'idparslist[[1]][c(1,4,7)] <- 1
#'idparslist[[1]][c(2,5,8)] <- 2
#'idparslist[[1]][c(3,6,9)] <- 3
#'idparslist[[2]][] <- 4
#'masterBlock <- matrix(c(5,6,5,6,5,6,5,6,5),ncol = 3,nrow = 3,byrow = TRUE)
#'diag(masterBlock) <- NA
#'diff.conceal <- FALSE
#'idparslist[[3]] <- q_doubletrans(traits,masterBlock,diff.conceal)
#'idparsfuncdefpar <- c(3,5,6)
#'idparsopt <- c(1,2)
#'idparsfix <- c(0,4)
#'initparsopt <- c(rep(intGuessLamba,2))
#'parsfix <- c(0,0)
#'idfactorsopt <- 1
#'initfactors <- 4
#'# functions_defining_params is a list of functions. Each function has no
#'# arguments and to refer
#'# to parameters ids should be indicated as "par_" i.e. par_3 refers to
#'# parameter 3. When a function is defined, be sure that all the parameters
#'# involved are either estimated, fixed or
#'# defined by previous functions (i.e, a function that defines parameter in
#'# 'functions_defining_params'). The user is responsible for this. In this
#'# exampl3, par_3 (i.e., parameter 3) is needed to calculate par_6. This is
#'# correct because par_3 is defined in
#'# the first function of 'functions_defining_params'. Notice that factor_1
#'# indicates a value that will be estimated to satisfy the equation. The same
#'# factor can be shared to define several parameters.
#'functions_defining_params <- list()
#'functions_defining_params[[1]] <- function(){
#'  par_3 <- par_1 + par_2
#'}
#'functions_defining_params[[2]] <- function(){
#'  par_5 <- par_1 * factor_1
#'}
#'functions_defining_params[[3]] <- function(){
#'  par_6 <- par_3 * factor_1
#'}
#'
#'tol = c(1e-03, 1e-03, 1e-03)
#'optimmethod = "simplex"
#'cond<-"proper_cond"
#'root_state_weight <- "proper_weights"
#'sampling_fraction <- c(1,1,1)
#'maxiter <- 10
#'model <- secsse_ml_func_def_pars(phylotree,
#'traits,
#'num_concealed_states,
#'idparslist,
#'idparsopt,
#'initparsopt,
#'idfactorsopt,
#'initfactors,
#'idparsfix,
#'parsfix,
#'idparsfuncdefpar,
#'functions_defining_params,
#'cond,
#'root_state_weight,
#'sampling_fraction,
#'tol,
#'maxiter,
#'optimmethod,
#'num_cycles = 1)
#'print(model$ML)
#'# ML -136.45265 if run till completion
#' @export
secsse_ml_func_def_pars <- function(phy,
                                    traits,
                                    num_concealed_states,
                                    idparslist,
                                    idparsopt,
                                    initparsopt,
                                    idfactorsopt,
                                    initfactors,
                                    idparsfix,
                                    parsfix,
                                    idparsfuncdefpar,
                                    functions_defining_params = NULL,
                                    cond = "proper_cond",
                                    root_state_weight = "proper_weights",
                                    sampling_fraction,
                                    tol = c(1E-4, 1E-5, 1E-7),
                                    maxiter = 1000 *
                                      round((1.25) ^ length(idparsopt)),
                                    optimmethod = "simplex",
                                    num_cycles = 1,
                                    loglik_penalty = 0,
                                    is_complete_tree = FALSE,
                                    take_into_account_root_edge = FALSE,
                                    num_threads = 1,
                                    atol = 1e-8,
                                    rtol = 1e-6,
                                    method = "odeint::runge_kutta_cash_karp54",
                                    return_root_state = FALSE) {
  return(master_ml(phy = phy,
                   traits = traits,
                   num_concealed_states = num_concealed_states,
                   idparslist = idparslist,
                   idparsopt = idparsopt,
                   initparsopt = initparsopt,
                   idparsfix = idparsfix,
                   parsfix = parsfix,
                   idfactorsopt = idfactorsopt,
                   initfactors = initfactors,
                   idparsfuncdefpar = idparsfuncdefpar,
                   functions_defining_params = functions_defining_params,
                   cond = cond,
                   root_state_weight = root_state_weight,
                   sampling_fraction = sampling_fraction,
                   tol = tol,
                   maxiter = maxiter,
                   optimmethod = optimmethod,
                   num_cycles = num_cycles,
                   loglik_penalty = loglik_penalty,
                   is_complete_tree = is_complete_tree,
                   take_into_account_root_edge = 
                     take_into_account_root_edge,
                   num_threads = num_threads,
                   atol = atol,
                   rtol = rtol,
                   method = method,
                   return_root_state = return_root_state))
}

#' Maximum likehood estimation for (SecSSE) with parameter as complex
#' functions. Cladogenetic version
#' 
#' Maximum likehood estimation under cla Several examined and concealed
#' States-dependent Speciation and Extinction (SecSSE) where some paramaters are
#' functions of other parameters and/or factors. Offers the option of
#' cladogenesis
#' 
#' @inheritParams default_params_doc
#' 
#' @return Parameter estimated and maximum likelihood
#' @examples
#'# Example of how to set the arguments for an ML search. The ML search is
#'# stopped after 10 iterations to keep the run time short.
#'library(secsse)
#'library(DDD)
#'set.seed(16)
#'phylotree <- ape::rbdtree(0.07,0.001,Tmax=50)
#'startingpoint <- bd_ML(brts = ape::branching.times(phylotree))
#'intGuessLamba <- startingpoint$lambda0
#'intGuessMu <- startingpoint$mu0
#'traits <-  sample(c(0,1,2),
#'                  ape::Ntip(phylotree), replace = TRUE) # get some traits
#'num_concealed_states <- 3
#'idparslist <- cla_id_paramPos(traits, num_concealed_states)
#'idparslist$lambdas[1,] <- c(1,2,3,1,2,3,1,2,3)
#'idparslist[[2]][] <- 4
#'masterBlock <- matrix(c(5,6,5,6,5,6,5,6,5),ncol = 3, nrow=3, byrow = TRUE)
#'diag(masterBlock) <- NA
#'diff.conceal <- FALSE
#'idparslist[[3]] <- q_doubletrans(traits,masterBlock,diff.conceal)
#'idparsfuncdefpar <- c(3,5,6)
#'idparsopt <- c(1,2)
#'idparsfix <- c(0,4)
#'initparsopt <- c(rep(intGuessLamba,2))
#'parsfix <- c(0,0)
#'idfactorsopt <- 1
#'initfactors <- 4
#'# functions_defining_params is a list of functions. Each function has no
#'# arguments and to refer
#'# to parameters ids should be indicated as 'par_' i.e. par_3 refers to
#'# parameter 3. When a
#'# function is defined, be sure that all the parameters involved are either
#'# estimated, fixed or
#'# defined by previous functions (i.e, a function that defines parameter in
#'# 'functions_defining_params'). The user is responsible for this. In this
#'# example, par_3
#'# (i.e., parameter 3) is needed to calculate par_6. This is correct because
#'# par_3 is defined
#'# in the first function of 'functions_defining_params'. Notice that factor_1
#'# indicates a value
#'# that will be estimated to satisfy the equation. The same factor can be
#'# shared to define several parameters.
#'functions_defining_params <- list()
#'functions_defining_params[[1]] <- function() {
#'  par_3 <- par_1 + par_2
#'}
#'functions_defining_params[[2]] <- function() {
#'  par_5 <- par_1 * factor_1
#'}
#'functions_defining_params[[3]] <- function() {
#'  par_6 <- par_3 * factor_1
#'}
#'
#'tol = c(1e-02, 1e-03, 1e-04)
#'optimmethod = 'subplex'
#'cond <- 'proper_cond'
#'root_state_weight <- 'proper_weights'
#'sampling_fraction <- c(1,1,1)
#'maxiter <- 10
#'model <- cla_secsse_ml_func_def_pars(phylotree,
#'traits,
#'num_concealed_states,
#'idparslist,
#'idparsopt,
#'initparsopt,
#'idfactorsopt,
#'initfactors,
#'idparsfix,
#'parsfix,
#'idparsfuncdefpar,
#'functions_defining_params,
#'cond,
#'root_state_weight,
#'sampling_fraction,
#'tol,
#'maxiter,
#'optimmethod,
#'num_cycles = 1)
#'# ML -136.5796 if run till completion
#' @export
cla_secsse_ml_func_def_pars <- function(phy,
                                        traits,
                                        num_concealed_states,
                                        idparslist,
                                        idparsopt,
                                        initparsopt,
                                        idfactorsopt,
                                        initfactors,
                                        idparsfix,
                                        parsfix,
                                        idparsfuncdefpar,
                                        functions_defining_params,
                                        cond = "proper_cond",
                                        root_state_weight = "proper_weights",
                                        sampling_fraction,
                                        tol = c(1e-04, 1e-05, 1e-07),
                                        maxiter = 1000 *
                                          round((1.25) ^ length(idparsopt)),
                                        optimmethod = "simplex",
                                        num_cycles = 1,
                                        loglik_penalty = 0,
                                        is_complete_tree = FALSE,
                                        take_into_account_root_edge = FALSE,
                                        verbose = TRUE,
                                        num_threads = 1,
                                        atol = 1e-12,
                                        rtol = 1e-12,
                                        method = "odeint::runge_kutta_cash_karp54",
                                        use_normalization = TRUE,
                                        return_root_state = FALSE) {
  return(master_ml(phy = phy,
                   traits = traits,
                   num_concealed_states = num_concealed_states,
                   idparslist = idparslist,
                   idparsopt = idparsopt,
                   initparsopt = initparsopt,
                   idparsfix = idparsfix,
                   parsfix = parsfix,
                   idfactorsopt = idfactorsopt,
                   initfactors = initfactors,
                   idparsfuncdefpar = idparsfuncdefpar,
                   functions_defining_params = functions_defining_params,
                   cond = cond,
                   root_state_weight = root_state_weight,
                   sampling_fraction = sampling_fraction,
                   tol = tol,
                   maxiter = maxiter,
                   optimmethod = optimmethod,
                   num_cycles = num_cycles,
                   loglik_penalty = loglik_penalty,
                   is_complete_tree = is_complete_tree,
                   take_into_account_root_edge = take_into_account_root_edge,
                   verbose = verbose,
                   num_threads = num_threads,
                   atol = atol,
                   rtol = rtol,
                   method = method,
                   use_normalization = use_normalization,
                   return_root_state = return_root_state))
}
