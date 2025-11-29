# Maximum likehood estimation for (SecSSE) complex functions as parameter

Maximum likehood estimation under Several examined and concealed
States-dependent Speciation and Extinction (SecSSE) where some
parameters are functions of other parameters and/or factors.

## Usage

``` r
secsse_ml_func_def_pars(
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
  functions_defining_params = NULL,
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
  num_threads = 1,
  atol = 1e-08,
  rtol = 1e-06,
  method = "odeint::runge_kutta_cash_karp54",
  return_root_state = FALSE
)
```

## Arguments

- phy:

  phylogenetic tree of class `phylo`, rooted and with branch lengths.
  Alternatively, multiple phylogenetic trees can be provided as the
  `multiPhylo` class.

- traits:

  vector with trait states for each tip in the phylogeny. The order of
  the states must be the same as the tree tips. For help, see
  [`vignette("starting_secsse", package = "secsse")`](https://rsetienne.github.io/secsse/articles/starting_secsse.md).
  When providing a `multiPhylo` set of multiple phylognies, traits
  should be a list where each entry in the list corresponds to the
  matching phylogeny on that position.

- num_concealed_states:

  number of concealed states, generally equivalent to the number of
  examined states in the dataset.

- idparslist:

  overview of parameters and their values.

- idparsopt:

  a numeric vector with the ID of parameters to be estimated.

- initparsopt:

  a numeric vector with the initial guess of the parameters to be
  estimated.

- idfactorsopt:

  id of the factors that will be optimized. There are not fixed factors,
  so use a constant within `functions_defining_params`.

- initfactors:

  the initial guess for a factor (it should be set to `NULL` when no
  factors).

- idparsfix:

  a numeric vector with the ID of the fixed parameters.

- parsfix:

  a numeric vector with the value of the fixed parameters.

- idparsfuncdefpar:

  id of the parameters which will be a function of optimized and/or
  fixed parameters. The order of id should match
  `functions_defining_params`.

- functions_defining_params:

  a list of functions. Each element will be a function which defines a
  parameter e.g. `id_3 <- (id_1 + id_2) / 2`. See example.

- cond:

  condition on the existence of a node root: `"maddison_cond"`,
  `"proper_cond"` (default). For details, see vignette.

- root_state_weight:

  the method to weigh the states: `"maddison_weights"`,
  `"proper_weights"` (default) or `"equal_weights"`. It can also be
  specified for the root state: the vector `c(1, 0, 0)` indicates state
  1 was the root state. When using a `multiPhylo` object,
  root_state_weight should be list where each entry in the list
  corresponds to the root_state_weight for each tree.

- sampling_fraction:

  vector that states the sampling proportion per trait state. It must
  have as many elements as there are trait states. When using a
  `multiPhylo` object, sampling fraction should be list where each entry
  in the list corresponds to the sampling proportion for each tree.

- tol:

  A numeric vector with the maximum tolerance of the optimization
  algorithm. Default is `c(1e-04, 1e-05, 1e-05)`.

- maxiter:

  max number of iterations. Default is
  `1000 * round((1.25) ^ length(idparsopt))`.

- optimmethod:

  A string with method used for optimization. Default is `"simplex"`.
  Alternative is `"subplex"`. Both are called from
  [`DDD::optimizer()`](https://rsetienne.github.io/DDD/reference/optimizer.html),
  simplex is implemented natively in DDD, while subplex is ultimately
  called from
  [`subplex::subplex()`](https://rdrr.io/pkg/subplex/man/subplex.html).

- num_cycles:

  Number of cycles of the optimization. When set to `Inf`, the
  optimization will be repeated until the result is, within the
  tolerance, equal to the starting values, with a maximum of 10 cycles.

- loglik_penalty:

  the size of the penalty for all parameters; default is 0 (no penalty).

- is_complete_tree:

  logical specifying whether or not a tree with all its extinct species
  is provided. If set to `TRUE`, it also assumes that all *all* extinct
  lineages are present on the tree. Defaults to `FALSE`.

- take_into_account_root_edge:

  if TRUE, the LL integration is continued along the root edge. This
  also affects conditioning (as now, conditioning no longer needs to
  assume a speciation event at the start of the tree)

- num_threads:

  number of threads to be used. Default is one thread.

- atol:

  A numeric specifying the absolute tolerance of integration.

- rtol:

  A numeric specifying the relative tolerance of integration.

- method:

  ODE integration method. Choose from:
  `"odeint::runge_kutta_cash_karp54"`,
  `"odeint::runge_kutta_fehlberg78"`, `"odeint::runge_kutta_dopri5"`,
  `"odeint::bulirsch_stoer"` and `"odeint::runge_kutta4"`. Default
  method is: `"odeint::runge_kutta_cash_karp54"`.

- return_root_state:

  if TRUE, returns the state of the system at the root, this can be
  useful to use as the starting point of a simulation. When used in ML,
  after finishing the ML optimization, the found optimum is evaluated
  one more time to retrieve the root state (to avoid having to store the
  root state every ML evaluation).

## Value

Parameter estimates and maximum likelihood

## Examples

``` r
# Example of how to set the arguments for an ML search. The ML search is stopped
# after 10 iterations to keep run time short.
library(secsse)
library(DDD)
set.seed(16)
phylotree <- ape::rbdtree(0.07,0.001,Tmax=50)
startingpoint<-bd_ML(brts = ape::branching.times(phylotree))
#> You are optimizing lambda0 mu0 
#> You are fixing lambda1 mu1 
#> Optimizing the likelihood - this may take a while. 
#> The loglikelihood for the initial parameter values is -106.425581129189.
#> 1 0.1 0.05 -106.425581129189 initial 
#> 2 0.090134529147982 0.0537544696066746 -105.907164659399 expand 
#> 3 0.090134529147982 0.0537544696066746 -105.907164659399 reflect 
#> 4 0.090134529147982 0.0537544696066746 -105.907164659399 contract outside 
#> 5 0.0828507795100222 0.0525 -105.824147025887 expand 
#> 6 0.0876957494407158 0.05 -105.783359852731 reflect 
#> 7 0.0756637168141592 0.0462721893491124 -105.702009387679 expand 
#> 8 0.0792452830188678 0.0395061728395062 -105.427240145262 expand 
#> 9 0.0792452830188678 0.0395061728395062 -105.427240145262 reflect 
#> 10 0.0685714285714284 0.0209006928406466 -104.954532462709 expand 
#> 11 0.0685714285714284 0.0209006928406466 -104.954532462709 contract outside 
#> 12 0.0665204277488344 0.00910258239406473 -104.839232262949 reflect 
#> 13 0.0581066376496188 0.00294951786727165 -104.62976274662 reflect 
#> 14 0.0581066376496188 0.00294951786727165 -104.62976274662 contract inside 
#> 15 0.0581066376496188 0.00294951786727165 -104.62976274662 reflect 
#> 16 0.0581066376496188 0.00294951786727165 -104.62976274662 contract inside 
#> 17 0.0581066376496188 0.00294951786727165 -104.62976274662 reflect 
#> 18 0.0581066376496188 0.00294951786727165 -104.62976274662 contract inside 
#> 19 0.0581066376496188 0.00294951786727165 -104.62976274662 reflect 
#> 20 0.0581066376496188 0.00294951786727165 -104.62976274662 reflect 
#> 21 0.0581066376496188 0.00294951786727165 -104.62976274662 contract inside 
#> 22 0.0581066376496188 0.00294951786727165 -104.62976274662 contract inside 
#> 23 0.0579203041153689 0.000953362994805256 -104.608294411085 expand 
#> 24 0.0579203041153689 0.000953362994805256 -104.608294411085 contract inside 
#> 25 0.0579203041153689 0.000953362994805256 -104.608294411085 contract inside 
#> 26 0.0569185700009541 0.000733154024100243 -104.601565921372 reflect 
#> 27 0.0569185700009541 0.000733154024100243 -104.601565921372 contract inside 
#> 28 0.0569185700009541 0.000733154024100243 -104.601565921372 reflect 
#> 29 0.0567732085054733 5.15034386826397e-05 -104.593729444779 reflect 
#> 30 0.0567732085054733 5.15034386826397e-05 -104.593729444779 contract outside 
#> 31 0.0567732085054733 5.15034386826397e-05 -104.593729444779 contract inside 
#> 32 0.0564072277699549 1.16314076479027e-05 -104.592993248356 reflect 
#> 33 0.0564072277699549 1.16314076479027e-05 -104.592993248356 contract inside 
#> 34 0.0564072277699549 1.16314076479027e-05 -104.592993248356 contract inside 
#> 35 0.0564072277699549 1.16314076479027e-05 -104.592993248356 contract inside 
#> 36 0.0564072277699549 1.16314076479027e-05 -104.592993248356 contract inside 
#> 37 0.0564072277699549 1.16314076479027e-05 -104.592993248356 reflect 
#> 38 0.0564072277699549 1.16314076479027e-05 -104.592993248356 contract inside 
#> 39 0.0564072277699549 1.16314076479027e-05 -104.592993248356 reflect 
#> 40 0.0566031943033848 5.26275596158352e-06 -104.592948906985 reflect 
#> 41 0.0566031943033848 5.26275596158352e-06 -104.592948906985 contract inside 
#> 42 0.0566031943033848 5.26275596158352e-06 -104.592948906985 contract inside 
#> 43 0.0564717106877172 3.06528647061127e-06 -104.592865779747 reflect 
#> 44 0.0564717106877172 3.06528647061127e-06 -104.592865779747 contract inside 
#> 45 0.0564717106877172 3.06528647061127e-06 -104.592865779747 contract inside 
#> 46 0.0565369808069373 5.39732631051912e-07 -104.592847596263 reflect 
#> 47 0.0565369808069373 5.39732631051912e-07 -104.592847596263 contract inside 
#> 48 0.0564877157586484 1.77106711409237e-08 -104.592829553347 reflect 
#> 49 0.0564877157586484 1.77106711409237e-08 -104.592829553347 contract inside 
#> 50 0.0564877157586484 1.77106711409237e-08 -104.592829553347 contract inside 
#> 51 0.0564877157586484 1.77106711409237e-08 -104.592829553347 reflect 
#> 52 0.0564877157586484 1.77106711409237e-08 -104.592829553347 contract inside 
#> 53 0.0564877157586484 1.77106711409237e-08 -104.592829553347 contract outside 
#> 54 0.0564877157586484 1.77106711409237e-08 -104.592829553347 contract inside 
#> 55 0.0564877157586484 1.77106711409237e-08 -104.592829553347 reflect 
#> 56 0.0564877157586484 1.77106711409237e-08 -104.592829553347 contract inside 
#> 57 0.0564877157586484 1.77106711409237e-08 -104.592829553347 contract inside 
#> 58 0.0564877157586484 1.77106711409237e-08 -104.592829553347 reflect 
#> 59 0.0564877157586484 1.77106711409237e-08 -104.592829553347 contract inside 
#> 60 0.0564877157586484 1.77106711409237e-08 -104.592829553347 contract inside 
#> 61 0.0564877157586484 1.77106711409237e-08 -104.592829553347 reflect 
#> 62 0.0564787828875198 9.44968080415927e-09 -104.592829468775 reflect 
#> 63 0.0564787828875198 9.44968080415927e-09 -104.592829468775 contract inside 
#> 64 0.0564787828875198 9.44968080415927e-09 -104.592829468775 contract inside 
#> 65 0.0564842332801909 1.66558900374773e-08 -104.592829465403 contract inside 
#> 66 0.056480297387231 4.35303382301016e-09 -104.592829359687 reflect 
#> 67 0.056480297387231 4.35303382301016e-09 -104.592829359687 reflect 
#> 68 0.056480297387231 4.35303382301016e-09 -104.592829359687 contract inside 
#> 69 0.056480297387231 4.35303382301016e-09 -104.592829359687 contract inside 
#> 70 0.0564805246767246 1.99140315115161e-09 -104.592829326339 reflect 
#> 71 0.0564805246767246 1.99140315115161e-09 -104.592829326339 contract inside 
#> 72 0.0564833919317324 4.11848852619872e-09 -104.592829316166 expand 
#> 73 0.0564833919317324 4.11848852619872e-09 -104.592829316166 contract inside 
#> 74 0.0564818708911459 1.30326827848874e-09 -104.592829293202 reflect 
#> 75 0.0564818708911459 1.30326827848874e-09 -104.592829293202 reflect 
#> 76 0.056483217108998 6.15133406772924e-10 -104.592829275411 reflect 
#> 77 0.056483217108998 6.15133406772924e-10 -104.592829275411 contract inside 
#> 78 0.056483217108998 6.15133406772924e-10 -104.592829275411 contract inside 
#> 79 0.056483217108998 6.15133406772924e-10 -104.592829275411 reflect 
#> 80 0.056483217108998 6.15133406772924e-10 -104.592829275411 contract inside 
#> 81 0.056483217108998 6.15133406772924e-10 -104.592829275411 contract inside 
#> 82 0.0564836455885223 3.39975282354587e-10 -104.592829272312 reflect 
#> 83 0.0564829739829893 6.57753981152692e-11 -104.592829269624 reflect 
#> 84 0.0564829739829893 6.57753981152692e-11 -104.592829269624 contract inside 
#> 85 0.0564829739829893 6.57753981152692e-11 -104.592829269624 contract inside 
#> 86 0.0564829739829893 6.57753981152692e-11 -104.592829269624 reflect 
#> 87 0.0564829739829893 6.57753981152692e-11 -104.592829269624 contract inside 
#> 88 0.0564829739829893 6.57753981152692e-11 -104.592829269624 reflect 
#> 89 0.0564834791538976 6.1041849912291e-11 -104.592829268828 expand 
#> 90 0.0564834791538976 6.1041849912291e-11 -104.592829268828 contract inside 
#> 91 0.0564834791538976 6.1041849912291e-11 -104.592829268828 contract inside 
#> 92 0.0564834791538976 6.1041849912291e-11 -104.592829268828 contract inside 
#> 93 0.0564834404619696 6.30926800524269e-13 -104.592829268103 expand 
#> 94 0.0564834404619696 6.30926800524269e-13 -104.592829268103 contract inside 
#> 95 0.0564834404619696 6.30926800524269e-13 -104.592829268103 contract inside 
#> 96 0.0564834404619696 6.30926800524269e-13 -104.592829268103 contract inside 
#> 97 0.0564834404619696 6.30926800524269e-13 -104.592829268103 contract inside 
#> 98 0.0564834404619696 6.30926800524269e-13 -104.592829268103 contract inside 
#> 99 0.0564834404619696 6.30926800524269e-13 -104.592829268103 contract inside 
#> 100 0.0564834404619696 6.30926800524269e-13 -104.592829268103 contract inside 
#> 101 0.0564834404619696 6.30926800524269e-13 -104.592829268103 contract inside 
#> 102 0.0564834404619696 6.30926800524269e-13 -104.592829268103 contract inside 
#> 103 0.0564834404619696 6.30926800524269e-13 -104.592829268103 contract inside 
#> 104 0.0564834404619696 6.30926800524269e-13 -104.592829268103 contract inside 
#> 105 0.0564834404619696 6.30926800524269e-13 -104.592829268103 contract inside 
#> 106 0.0564834404619696 6.30926800524269e-13 -104.592829268103 contract inside 
#> 107 0.0564834404619696 6.30926800524269e-13 -104.592829268103 contract inside 
#> 108 0.0564834404619696 6.30926800524269e-13 -104.592829268103 contract inside 
#> 109 0.0564834404619696 6.30926800524269e-13 -104.592829268103 contract inside 
#> 110 0.0564834404619696 6.30926800524269e-13 -104.592829268103 contract inside 
#> 111 0.0564834411366911 1.0121875039335e-14 -104.592829268096 reflect 
#> 112 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 113 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 114 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 115 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 116 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 117 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 118 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 119 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 120 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 121 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 122 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 123 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 124 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 125 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 126 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 127 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 128 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 129 0.0564834411366911 1.0121875039335e-14 -104.592829268096 contract inside 
#> 130 0.0564834411466793 1.05558453530427e-15 -104.592829268096 reflect 
#> 131 0.0564834411466793 1.05558453530427e-15 -104.592829268096 contract inside 
#> 132 0.0564834411466793 1.05558453530427e-15 -104.592829268096 contract inside 
#> 133 0.0564834411466793 1.05558453530427e-15 -104.592829268096 contract inside 
#> 134 0.0564834411466793 1.05558453530427e-15 -104.592829268096 contract inside 
#> 135 0.0564834411466793 1.05558453530427e-15 -104.592829268096 reflect 
#> 136 0.0564834411466793 1.05558453530427e-15 -104.592829268096 contract inside 
#> 137 0.0564834411466793 1.05558453530427e-15 -104.592829268096 shrink 
#> 138 0.0564834411466793 1.05558453530427e-15 -104.592829268096 contract inside 
#> 139 0.0564834411466793 1.05558453530427e-15 -104.592829268096 shrink 
#> 140 0.0564834411466793 1.05558453530427e-15 -104.592829268096 reflect 
#> 141 0.0564834411466793 1.05558453530427e-15 -104.592829268096 contract outside 
#> 142 0.0564834411466793 1.05558453530427e-15 -104.592829268096 shrink 
#> 143 0.0564834411466793 1.05558453530427e-15 -104.592829268096 reflect 
#> 144 0.0564834411465597 1.1641772011865e-15 -104.592829268096 shrink 
#> 145 0.0564834411465597 1.1641772011865e-15 -104.592829268096 shrink 
#> 146 0.0564834411465597 1.1641772011865e-15 -104.592829268096 shrink 
#> 147 0.0564834411465597 1.1641772011865e-15 -104.592829268096 shrink 
#> 148 0.0564834411465597 1.1641772011865e-15 -104.592829268096 contract outside 
#> 149 0.0564834411465597 1.1641772011865e-15 -104.592829268096 shrink 
#> 150 0.0564834411465597 1.1641772011865e-15 -104.592829268096 contract inside 
#> 151 0.0564834411465597 1.1641772011865e-15 -104.592829268096 shrink 
#> 152 0.0564834411465597 1.1641772011865e-15 -104.592829268096 shrink 
#> 153 0.0564834411465597 1.1641772011865e-15 -104.592829268096 reflect 
#> 154 0.0564834411465597 1.1641772011865e-15 -104.592829268096 contract outside 
#> Optimization has terminated successfully. 
#> 
#> Maximum likelihood parameter estimates: lambda0: 0.056483, mu0: 0.000000, lambda1: 0.000000, mu1: 0.000000: 
#> Maximum loglikelihood: -104.592829
intGuessLamba <- startingpoint$lambda0
intGuessMu <- startingpoint$mu0
traits <- sample(c(0,1,2), ape::Ntip(phylotree),replace=TRUE) #get some traits
num_concealed_states<-3
idparslist<-id_paramPos(traits, num_concealed_states)
idparslist[[1]][c(1,4,7)] <- 1
idparslist[[1]][c(2,5,8)] <- 2
idparslist[[1]][c(3,6,9)] <- 3
idparslist[[2]][] <- 4
masterBlock <- matrix(c(5,6,5,6,5,6,5,6,5),ncol = 3,nrow = 3,byrow = TRUE)
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
# to parameters ids should be indicated as "par_" i.e. par_3 refers to
# parameter 3. When a function is defined, be sure that all the parameters
# involved are either estimated, fixed or
# defined by previous functions (i.e, a function that defines parameter in
# 'functions_defining_params'). The user is responsible for this. In this
# exampl3, par_3 (i.e., parameter 3) is needed to calculate par_6. This is
# correct because par_3 is defined in
# the first function of 'functions_defining_params'. Notice that factor_1
# indicates a value that will be estimated to satisfy the equation. The same
# factor can be shared to define several parameters.
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

tol = c(1e-03, 1e-03, 1e-03)
optimmethod = "simplex"
cond<-"proper_cond"
root_state_weight <- "proper_weights"
sampling_fraction <- c(1,1,1)
maxiter <- 10
model <- secsse_ml_func_def_pars(phylotree,
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
#> Note: you set some transitions as impossible to happen.
#> Warning: Optimization has not converged. Try again with different initial values or increase the number of iterations.
print(model$ML)
#> NULL
# ML -136.45265 if run till completion
```
