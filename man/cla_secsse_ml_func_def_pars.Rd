% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cla_secsse_ml_func_def_pars.R
\name{cla_secsse_ml_func_def_pars}
\alias{cla_secsse_ml_func_def_pars}
\title{Maximum likehood estimation for (SecSSE) with parameter as complex functions. Cladogenetic version}
\usage{
cla_secsse_ml_func_def_pars(phy, traits, num_concealed_states, idparslist,
  idparsopt, initparsopt, idfactosopt, initfactos, idparsfix, parsfix,
  idparsfuncdefpar, functions_defining_params, cond = "proper_cond",
  root_state_weight = "proper_weights", sampling_fraction,
  tol = c(1e-04, 1e-05, 1e-07), maxiter = 1000 *
  round((1.25)^length(idparsopt)), use_fortran = TRUE,
  methode = "ode45", optimmethod = "simplex", run_parallel = FALSE)
}
\arguments{
\item{phy}{phylogenetic tree of class phylo, ultrametric, rooted and with branch lengths.}

\item{traits}{a vector with trait states for each tip in the phylogeny.}

\item{num_concealed_states}{number of concealed states, generally equivalent to the number of examined states in the dataset.}

\item{idparslist}{overview of parameters and their values.}

\item{idparsopt}{id of parameters to be estimated.}

\item{initparsopt}{initial guess of the parameters to be estimated.}

\item{idfactosopt}{id of the factors that will be optimized. There are not fixed factors, so use a constant within 'functions_defining_params'.}

\item{initfactos}{the initial guess for a factor (it should be set to NULL when no factors).}

\item{idparsfix}{id of the fixed parameters (it should be set to NULL when no factors).}

\item{parsfix}{value of the fixed parameters.}

\item{idparsfuncdefpar}{id of the parameters which will be a function of optimized and/or fixed parameters. The order of id should match functions_defining_params}

\item{functions_defining_params}{a list of functions. Each element will be a function which defines a parameter e.g. id_3 <- (id_1+id_2)/2. See example and vigenette}

\item{cond}{condition on the existence of a node root: "maddison_cond","proper_cond"(default). For details, see vignette.}

\item{root_state_weight}{the method to weigh the states:"maddison_weights","proper_weights"(default) or "equal_weights". It can also be specified the root state:the vector c(1,0,0) indicates state 1 was the root state.}

\item{sampling_fraction}{vector that states the sampling proportion per trait state. It must have as many elements as there are trait states.}

\item{tol}{maximum tolerance. Default is "c(1e-04, 1e-05, 1e-05)".}

\item{maxiter}{max number of iterations. Default is "1000 *round((1.25)^length(idparsopt))".}

\item{use_fortran}{Should the Fortran code for numerical integration be called? Default is TRUE.}

\item{methode}{method used for integration calculation. Default is "ode45".}

\item{optimmethod}{method used for optimization. Default is "simplex".}

\item{run_parallel}{should the routine to run in parallel be called? Read note below}
}
\value{
Parameter estimated and maximum likelihood
}
\description{
Maximum likehood estimation under cla Several examined and concealed States-dependent Speciation and Extinction (SecSSE) where some paramaters are functions of other parameters and/or factors. Offers the option of cladogenesis
}
\note{
To run in parallel it is needed to load the following libraries when windows: apTreeshape, doparallel and foreach. When unix, it requires: apTreeshape, doparallel, foreach and doMC
}
\examples{
# Example of how to set the arguments for a ML search.
rm(list=ls(all=TRUE))
library(secsse)
library(DDD)
set.seed(16)
phylotree <- ape::rbdtree(0.07,0.001,Tmax=50)
startingpoint<-bd_ML(brts = ape::branching.times(phylotree))
intGuessLamba <- startingpoint$lambda0
intGuessMu <- startingpoint$mu0
traits <-  sample(c(0,1,2), ape::Ntip(phylotree),replace=TRUE) #get some traits
num_concealed_states<-3
idparslist<-cla_id_paramPos(traits,num_concealed_states)
idparslist$lambdas[1,]<-c(1,2,3,1,2,3,1,2,3)
idparslist[[2]][]<-4
masterBlock<-matrix(c(5,6),ncol=3,nrow=3,byrow=TRUE)
diag(masterBlock)<-NA
diff.conceal <- FALSE
idparslist[[3]]<-q_doubletrans(traits,masterBlock,diff.conceal)
idparsfuncdefpar<-c(3,5,6)
idparsopt<-c(1,2)
idparsfix<-c(0,4)
initparsopt<-c(rep(intGuessLamba,2))
parsfix<-c(0,0)
idfactosopt<-1
initfactos<- 4
# functions_defining_params is a list of functions. Each function has no arguments and to refer
# to parameters ids should be indicated as "par_" i.e. par_3 refers to parameter 3. When a
# function is defined, be sure that all the parameters involved are either estimated, fixed or
# defined by previous functions (i.e, a function that defines parameter in
# 'functions_defining_params'). The user is responsible for this. In this example, par_3
# (i.e., parameter 3) is needed to calculate par_6. This is correct because par_3 is defined
# in the first function of 'functions_defining_params'. Notice that factor_1 indicates a value
# that will be estimated to satisfy the equation. The same factor can be shared to define
# several parameters.
functions_defining_params <- list()
functions_defining_params[[1]] <- function(){
 par_3 <- par_1 + par_2
}
functions_defining_params[[2]] <- function(){
 par_5 <- par_1 * factor_1
}
functions_defining_params[[3]] <- function(){
 par_6 <- par_3 * factor_1
}

tol = c(1e-04, 1e-05, 1e-07)
maxiter = 1000 * round((1.25)^length(idparsopt))
use_fortran = TRUE
methode = "ode45"
optimmethod = "simplex"
run_parallel = FALSE
cond<-"proper_cond"
root_state_weight <- "proper_weights"
sampling_fraction<-c(1,1,1)
#model<-cla_secsse_ml_func_def_pars(phylotree, traits, num_concealed_states, idparslist,
#                                   idparsopt, initparsopt, idfactosopt, initfactos,
#                                   idparsfix, parsfix, idparsfuncdefpar,
#                                   functions_defining_params, cond, root_state_weight,
#                                   sampling_fraction, tol, maxiter, use_fortran,
#                                   methode, optimmethod, run_parallel)

# ML -136.5796
}