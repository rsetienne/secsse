% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/secsse_ml.R
\name{cla_secsse_ml}
\alias{cla_secsse_ml}
\title{Maximum likehood estimation for (SecSSE)}
\usage{
cla_secsse_ml(
  phy,
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
  take_into_account_root_edge = FALSE,
  verbose = FALSE,
  num_threads = 1,
  atol = 1e-08,
  rtol = 1e-07,
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

\item{idparsfix}{a numeric vector with the ID of the fixed parameters.}

\item{parsfix}{a numeric vector with the value of the fixed parameters.}

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
Maximum likehood estimation under Several examined and concealed
States-dependent Speciation and Extinction (SecSSE) with cladogenetic option
}
\examples{
# Example of how to set the arguments for a ML search.
library(secsse)
library(DDD)
set.seed(13)
# Check the vignette for a better working exercise.
# lambdas for 0A and 1A and 2A are the same but need to be estimated
# (CTD model, see Syst Biol paper)
# mus are fixed to zero,
# the transition rates are constrained to be equal and fixed 0.01
phylotree <- ape::rcoal(31, tip.label = 1:31)
#get some traits
traits <-  sample(c(0,1,2), ape::Ntip(phylotree), replace = TRUE)
num_concealed_states <- 3
idparslist <- cla_id_paramPos(traits,num_concealed_states)
idparslist$lambdas[1,] <- c(1,1,1,2,2,2,3,3,3)
idparslist[[2]][] <- 4
masterBlock <- matrix(5,ncol = 3,nrow = 3,byrow = TRUE)
diag(masterBlock) <- NA
diff.conceal <- FALSE
idparslist[[3]] <- q_doubletrans(traits,masterBlock,diff.conceal)
startingpoint <- bd_ML(brts = ape::branching.times(phylotree))
intGuessLamba <- startingpoint$lambda0
intGuessMu <- startingpoint$mu0
idparsopt <- c(1,2,3)
initparsopt <- c(rep(intGuessLamba,3))
idparsfix <- c(0,4,5)
parsfix <- c(0,0,0.01)
tol <- c(1e-03, 1e-03, 1e-03)
maxiter <- 1000 * round((1.25) ^ length(idparsopt))
optimmethod <- 'simplex'
cond <- 'proper_cond'
root_state_weight <- 'proper_weights'
sampling_fraction <- c(1,1,1)
model <- cla_secsse_ml(
 phylotree,
 traits,
 num_concealed_states,
 idparslist,
 idparsopt,
 initparsopt,
 idparsfix,
 parsfix,
 cond,
 root_state_weight,
 sampling_fraction,
 tol,
 maxiter,
 optimmethod,
 num_cycles = 1,
 num_threads = 1,
 verbose = FALSE)
# [1] -90.9763
}
