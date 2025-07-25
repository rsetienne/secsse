% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/secsse_ml_func_def_pars.R
\name{cla_secsse_ml_func_def_pars}
\alias{cla_secsse_ml_func_def_pars}
\title{Maximum likehood estimation for (SecSSE) with parameter as complex
functions. Cladogenetic version}
\usage{
cla_secsse_ml_func_def_pars(
  phy,
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
  maxiter = 1000 * round((1.25)^length(idparsopt)),
  optimmethod = "simplex",
  num_cycles = 1,
  loglik_penalty = 0,
  is_complete_tree = FALSE,
  take_into_account_root_edge = FALSE,
  verbose = TRUE,
  num_threads = 1,
  atol = 1e-12,
  rtol = 1e-12,
  method = "odeint::bulirsch_stoer",
  use_normalization = TRUE
)
}
\arguments{
\item{phy}{phylogenetic tree of class \code{phylo}, rooted and with
branch lengths. Alternatively, multiple phylogenetic trees can be provided
as the \code{multiPhylo} class.}

\item{traits}{vector with trait states for each tip in the phylogeny. The
order of the states must be the same as the tree tips. For help, see
\code{vignette("starting_secsse", package = "secsse")}. When providing a
\code{multiPhylo} set of multiple phylognies, traits should be a list where
each entry in the list corresponds to the matching phylogeny on that
position.}

\item{num_concealed_states}{number of concealed states, generally equivalent
to the number of examined states in the dataset.}

\item{idparslist}{overview of parameters and their values.}

\item{idparsopt}{a numeric vector with the ID of parameters to be estimated.}

\item{initparsopt}{a numeric vector with the initial guess of the parameters
to be estimated.}

\item{idfactorsopt}{id of the factors that will be optimized. There are not
fixed factors, so use a constant within \code{functions_defining_params}.}

\item{initfactors}{the initial guess for a factor (it should be set to \code{NULL}
when no factors).}

\item{idparsfix}{a numeric vector with the ID of the fixed parameters.}

\item{parsfix}{a numeric vector with the value of the fixed parameters.}

\item{idparsfuncdefpar}{id of the parameters which will be a function of
optimized and/or fixed parameters. The order of id should match
\code{functions_defining_params}.}

\item{functions_defining_params}{a list of functions. Each element will be a
function which defines a parameter e.g. \code{id_3 <- (id_1 + id_2) / 2}. See
example.}

\item{cond}{condition on the existence of a node root: \code{"maddison_cond"},
\code{"proper_cond"} (default). For details, see vignette.}

\item{root_state_weight}{the method to weigh the states:
\code{"maddison_weights"}, \code{"proper_weights"} (default) or \code{"equal_weights"}.
It can also be specified for the root state: the vector \code{c(1, 0, 0)}
indicates state 1 was the root state. When
using a \code{multiPhylo} object, root_state_weight should be list where each
entry in the list corresponds to the root_state_weight for each tree.}

\item{sampling_fraction}{vector that states the sampling proportion per
trait state. It must have as many elements as there are trait states. When
using a \code{multiPhylo} object, sampling fraction should be list where each
entry in the list corresponds to the sampling proportion for each tree.}

\item{tol}{A numeric vector with the maximum tolerance of the optimization
algorithm. Default is \code{c(1e-04, 1e-05, 1e-05)}.}

\item{maxiter}{max number of iterations. Default is
\code{1000 * round((1.25) ^ length(idparsopt))}.}

\item{optimmethod}{A string with method used for optimization. Default is
\code{"simplex"}. Alternative is \code{"subplex"}. Both are called from
\code{\link[DDD:optimizer]{DDD::optimizer()}}, simplex is implemented natively in \pkg{DDD}, while subplex
is ultimately called from \code{\link[subplex:subplex]{subplex::subplex()}}.}

\item{num_cycles}{Number of cycles of the optimization. When set to \code{Inf},
the optimization will be repeated until the result is, within the
tolerance, equal to the starting values, with a maximum of 10 cycles.}

\item{loglik_penalty}{the size of the penalty for all parameters; default is
0 (no penalty).}

\item{is_complete_tree}{logical specifying whether or not a tree with all its
extinct species is provided. If set to \code{TRUE}, it also assumes that all
\emph{all} extinct lineages are present on the tree. Defaults to \code{FALSE}.}

\item{take_into_account_root_edge}{if TRUE, the LL integration is continued
along the root edge. This also affects conditioning (as now, conditioning
no longer needs to assume a speciation event at the start of the tree)}

\item{verbose}{sets verbose output; default is \code{TRUE} when \code{optimmethod} is
\code{"simplex"}. If \code{optimmethod} is set to \code{"simplex"}, then even if set to
\code{FALSE}, optimizer output will be shown.}

\item{num_threads}{number of threads to be used. Default is one thread.}

\item{atol}{A numeric specifying the absolute tolerance of integration.}

\item{rtol}{A numeric specifying the relative tolerance of integration.}

\item{method}{integration method used, available are:
\code{"odeint::runge_kutta_cash_karp54"}, \code{"odeint::runge_kutta_fehlberg78"},
\code{"odeint::runge_kutta_dopri5"}, \code{"odeint::bulirsch_stoer"} and
\code{"odeint::runge_kutta4"}. Default method is: \code{"odeint::bulirsch_stoer"}.}

\item{use_normalization}{normalize the density vector during integration,
more accurate but slower (default = TRUE)}
}
\value{
Parameter estimated and maximum likelihood
}
\description{
Maximum likehood estimation under cla Several examined and concealed
States-dependent Speciation and Extinction (SecSSE) where some paramaters are
functions of other parameters and/or factors. Offers the option of
cladogenesis
}
\examples{
# Example of how to set the arguments for an ML search. The ML search is
# stopped after 10 iterations to keep the run time short.
library(secsse)
library(DDD)
set.seed(16)
phylotree <- ape::rbdtree(0.07,0.001,Tmax=50)
startingpoint <- bd_ML(brts = ape::branching.times(phylotree))
intGuessLamba <- startingpoint$lambda0
intGuessMu <- startingpoint$mu0
traits <-  sample(c(0,1,2),
                 ape::Ntip(phylotree), replace = TRUE) # get some traits
num_concealed_states <- 3
idparslist <- cla_id_paramPos(traits, num_concealed_states)
idparslist$lambdas[1,] <- c(1,2,3,1,2,3,1,2,3)
idparslist[[2]][] <- 4
masterBlock <- matrix(c(5,6,5,6,5,6,5,6,5),ncol = 3, nrow=3, byrow = TRUE)
diag(masterBlock) <- NA
diff.conceal <- FALSE
idparslist[[3]] <- q_doubletrans(traits,masterBlock,diff.conceal)
idparsfuncdefpar <- c(3,5,6)
idparsopt <- c(1,2)
idparsfix <- c(0,4)
initparsopt <- c(rep(intGuessLamba,2))
parsfix <- c(0,0)
idfactorsopt <- 1
initfactors <- 4
# functions_defining_params is a list of functions. Each function has no
# arguments and to refer
# to parameters ids should be indicated as 'par_' i.e. par_3 refers to
# parameter 3. When a
# function is defined, be sure that all the parameters involved are either
# estimated, fixed or
# defined by previous functions (i.e, a function that defines parameter in
# 'functions_defining_params'). The user is responsible for this. In this
# example, par_3
# (i.e., parameter 3) is needed to calculate par_6. This is correct because
# par_3 is defined
# in the first function of 'functions_defining_params'. Notice that factor_1
# indicates a value
# that will be estimated to satisfy the equation. The same factor can be
# shared to define several parameters.
functions_defining_params <- list()
functions_defining_params[[1]] <- function() {
 par_3 <- par_1 + par_2
}
functions_defining_params[[2]] <- function() {
 par_5 <- par_1 * factor_1
}
functions_defining_params[[3]] <- function() {
 par_6 <- par_3 * factor_1
}

tol = c(1e-02, 1e-03, 1e-04)
optimmethod = 'subplex'
cond <- 'proper_cond'
root_state_weight <- 'proper_weights'
sampling_fraction <- c(1,1,1)
maxiter <- 10
model <- cla_secsse_ml_func_def_pars(phylotree,
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
cond,
root_state_weight,
sampling_fraction,
tol,
maxiter,
optimmethod,
num_cycles = 1)
# ML -136.5796 if run till completion
}
