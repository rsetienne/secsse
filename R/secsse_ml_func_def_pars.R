#' Maximum likehood estimation under Several examined and concealed
#' States-dependent Speciation and Extinction (SecSSE) where some paramaters
#' are functions of other parameters and/or factors.
#' @title Maximum likehood estimation for (SecSSE) with parameter as complex
#' functions.
#' @param phy phylogenetic tree of class phylo, ultrametric, rooted and with
#' branch lengths.
#' @param traits a vector with trait states for each tip in the phylogeny.
#' @param num_concealed_states number of concealed states, generally equivalent
#' to the number of examined states in the dataset.
#' @param idparslist overview of parameters and their values.
#' @param idparsopt id of parameters to be estimated.
#' @param initparsopt initial guess of the parameters to be estimated.
#' @param idfactorsopt id of the factors that will be optimized. There are not
#' fixed factors, so use a constant within 'functions_defining_params'.
#' @param initfactors the initial guess for a factor (it should be set to NULL
#' when no factors).
#' @param idparsfix id of the fixed parameters (it should be set to NULL when
#' there are no factors).
#' @param parsfix value of the fixed parameters.
#' @param idparsfuncdefpar id of the parameters which will be a function of
#' optimized and/or fixed parameters. The order of id should match
#' functions_defining_params
#' @param functions_defining_params a list of functions. Each element will be a
#' function which defines a parameter e.g. id_3 <- (id_1+id_2)/2. See example
#' and vigenette
#' @param cond condition on the existence of a node root:
#' "maddison_cond","proper_cond"(default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:
#' "maddison_weights","proper_weights"(default) or "equal_weights". It can also
#' be specified the root state:the vector c(1, 0, 0) indicates state
#' 1 was the root state.
#' @param sampling_fraction vector that states the sampling proportion per trait
#' state. It must have as many elements as there are trait states.
#' @param tol maximum tolerance. Default is "c(1e-04, 1e-05, 1e-05)".
#' @param maxiter max number of iterations. Default is
#' "1000 *round((1.25)^length(idparsopt))".
#' @param optimmethod method used for optimization. Default is "subplex".
#' @param num_cycles number of cycles of the optimization (default is 1).
#' @param loglik_penalty the size of the penalty for all parameters;
#' default is 0 (no penalty)
#' @param is_complete_tree whether or not a tree with all its extinct species
#' is provided
#' @param num_threads number of threads. Set to -1 to use all available threads.
#' Default is one thread.
#' @param atol absolute tolerance of integration
#' @param rtol relative tolerance of integration
#' @param method integration method used, available are:
#' "odeint::runge_kutta_cash_karp54", "odeint::runge_kutta_fehlberg78",
#' "odeint::runge_kutta_dopri5", "odeint::bulirsch_stoer" and
#' "odeint::runge_kutta4". Default method is:"odeint::bulirsch_stoer".
#' @return Parameter estimated and maximum likelihood
#' @return Parameter estimated and maximum likelihood
#' @examples
#'# Example of how to set the arguments for a ML search.
#'rm(list=ls(all=TRUE))
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
#'tol = c(1e-02, 1e-03, 1e-04)
#'maxiter = 1000 * round((1.25)^length(idparsopt))
#'optimmethod = "subplex"
#'cond<-"proper_cond"
#'root_state_weight <- "proper_weights"
#'sampling_fraction <- c(1,1,1)
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
#'# ML -136.5796
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
                                    optimmethod = "subplex",
                                    num_cycles = 1,
                                    loglik_penalty = 0,
                                    is_complete_tree = FALSE,
                                    num_threads = 1,
                                    atol = 1e-8,
                                    rtol = 1e-6,
                                    method = "odeint::bulirsch_stoer") {
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
                   num_threads = num_threads,
                   atol = atol,
                   rtol = rtol,
                   method = method))
} 


#' Maximum likehood estimation under cla Several examined and concealed
#' States-dependent Speciation and Extinction (SecSSE) where some paramaters are
#' functions of other parameters and/or factors. Offers the option of
#' cladogenesis
#' @title Maximum likehood estimation for (SecSSE) with parameter as complex
#' functions. Cladogenetic version
#' @param phy phylogenetic tree of class phylo, ultrametric, rooted and with
#' branch lengths.
#' @param traits a vector with trait states for each tip in the phylogeny.
#' @param num_concealed_states number of concealed states, generally equivalent
#' to the number of examined states in the dataset.
#' @param idparslist overview of parameters and their values.
#' @param idparsopt id of parameters to be estimated.
#' @param initparsopt initial guess of the parameters to be estimated.
#' @param idfactorsopt id of the factors that will be optimized. There are not
#' fixed factors, so use a constant within 'functions_defining_params'.
#' @param initfactors the initial guess for a factor (it should be set to NULL
#' when no factors).
#' @param idparsfix id of the fixed parameters (it should be set to NULL when
#' no factors).
#' @param parsfix value of the fixed parameters.
#' @param idparsfuncdefpar id of the parameters which will be a function of
#' optimized and/or fixed parameters. The order of id should match
#' functions_defining_params
#' @param functions_defining_params a list of functions. Each element will be a
#' function which defines a parameter e.g. id_3 <- (id_1+id_2)/2. See example
#' and vigenette
#' @param cond condition on the existence of a node root: 'maddison_cond',
#' 'proper_cond'(default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:'maddison_weights',
#' 'proper_weights'(default) or 'equal_weights'. It can also be specified the
#' root
#' state:the vector c(1,0,0) indicates state 1 was the root state.
#' @param sampling_fraction vector that states the sampling proportion per trait
#' state. It must have as many elements as there are trait states.
#' @param tol maximum tolerance. Default is 'c(1e-04, 1e-05, 1e-05)'.
#' @param maxiter max number of iterations. Default is
#' '1000*round((1.25)^length(idparsopt))'.
#' @param optimmethod method used for optimization. Default is 'simplex'.
#' @param num_cycles number of cycles of the optimization (default is 1).
#' @param loglik_penalty the size of the penalty for all parameters; default
#' is 0 (no penalty)
#' @param is_complete_tree whether or not a tree with all its extinct species
#' is provided
#' @param verbose sets verbose output; default is verbose when optimmethod is
#' 'subplex'
#' @param num_threads number of threads. Set to -1 to use all available
#' threads. Default is one thread.
#' @param atol absolute tolerance of integration
#' @param rtol relative tolerance of integration
#' @param method integration method used, available are:
#' "odeint::runge_kutta_cash_karp54", "odeint::runge_kutta_fehlberg78",
#' "odeint::runge_kutta_dopri5", "odeint::bulirsch_stoer" and
#' "odeint::runge_kutta4". Default method is:"odeint::bulirsch_stoer".
#' @return Parameter estimated and maximum likelihood
#' @return Parameter estimated and maximum likelihood
#' @examples
#'# Example of how to set the arguments for a ML search.
#'rm(list=ls(all=TRUE))
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
#'maxiter = 1000 * round((1.25)^length(idparsopt))
#'optimmethod = 'subplex'
#'cond <- 'proper_cond'
#'root_state_weight <- 'proper_weights'
#'sampling_fraction <- c(1,1,1)
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
#'# ML -136.5796
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
                                        verbose = (optimmethod == "subplex"),
                                        num_threads = 1,
                                        atol = 1e-12,
                                        rtol = 1e-12,
                                        method = "odeint::bulirsch_stoer") {
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
                   num_threads = num_threads,
                   atol = atol,
                   rtol = rtol,
                   method = method))
} 
