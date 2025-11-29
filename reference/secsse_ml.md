# Maximum likehood estimation for (SecSSE)

Maximum likehood estimation under Several examined and concealed trait
States dependent Speciation and Extinction (SecSSE)

## Usage

``` r
secsse_ml(
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
  method = "odeint::runge_kutta_cash_karp54",
  use_normalization = TRUE,
  return_root_state = FALSE,
  see_ancestral_states = FALSE
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

- idparsfix:

  a numeric vector with the ID of the fixed parameters.

- parsfix:

  a numeric vector with the value of the fixed parameters.

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

- verbose:

  sets verbose output; default is `TRUE`.

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

- use_normalization:

  normalize the density vector during integration, more accurate but
  slower (default = TRUE)

- return_root_state:

  if TRUE, returns the state of the system at the root, this can be
  useful to use as the starting point of a simulation. When used in ML,
  after finishing the ML optimization, the found optimum is evaluated
  one more time to retrieve the root state (to avoid having to store the
  root state every ML evaluation).

- see_ancestral_states:

  Boolean for whether the ancestral states for each of the internal
  nodes should be output. Defaults to `FALSE`.

## Value

A list with the following elements \$MLpars: the maximum likelihood
parameter estimates \$ML: the maximum likelihood of the data
(phylogeny + tip states) given the parameters (speciation, extinction,
transition rates). \$conv: whether the optimization converged or not If
see_ancestral_states = TRUE, then there will be two additional elements:
\$ancestral_states: a matrix with the probabilities of each state at the
internal nodes \$states: a matrix with the probabilities E, D
(normalized) and S that are used in the calculations. The
ancestral_states matrix is a submatrix of this matrix. This matrix is
mostly used for package developers. If return_root_state = TRUE, then
there will be one additional element: \$root_state: vector with
probabilities of each state at the root. This vector is the same as the
top row of \$ancestral_states We have used the shorthand description of
"probabilities of each state", but technically, the probabilities are
the normalized probabilities D of the data given each state at the
internal nodes.

## Examples

``` r
# Example of how to set the arguments for a ML search.
library(secsse)
library(DDD)
set.seed(13)
# lambdas for 0A and 1A and 2A are the same but need to be estimated
# mus are fixed to
# the transition rates are constrained to be equal and fixed 0.01
phylotree <- ape::rcoal(31, tip.label = 1:31)
traits <-  sample(c(0,1,2), ape::Ntip(phylotree),replace=TRUE)#get some traits
num_concealed_states<-3
idparslist <- id_paramPos(traits, num_concealed_states)
idparslist[[1]][c(1,4,7)] <- 1
idparslist[[1]][c(2,5,8)] <- 2
idparslist[[1]][c(3,6,9)] <- 3
idparslist[[2]][]<-4
masterBlock <- matrix(5,ncol = 3,nrow = 3,byrow = TRUE)
diag(masterBlock) <- NA
diff.conceal <- FALSE
idparslist[[3]] <- q_doubletrans(traits,masterBlock,diff.conceal)
startingpoint <- DDD::bd_ML(brts = ape::branching.times(phylotree))
#> You are optimizing lambda0 mu0 
#> You are fixing lambda1 mu1 
#> Optimizing the likelihood - this may take a while. 
#> The loglikelihood for the initial parameter values is -67.5669039297515.
#> 1 0.1 0.05 -67.5669039297515 initial 
#> 2 0.10751708428246 0.0450354609929078 -65.4917085799141 expand 
#> 3 0.118987341772152 0.0425707547169811 -62.6116799074596 expand 
#> 4 0.130040674026729 0.0316219369894982 -60.0623957590543 expand 
#> 5 0.16003578884581 0.0214904679376082 -54.1880637307688 expand 
#> 6 0.171919252786984 0.0109777015437391 -52.1459129951388 reflect 
#> 7 0.20421052631579 0.00124575311438257 -47.3167705826442 reflect 
#> 8 0.20421052631579 0.00124575311438257 -47.3167705826442 contract inside 
#> 9 0.224058848983124 0.000466793974113919 -44.7449478078427 expand 
#> 10 0.224058848983124 0.000466793974113919 -44.7449478078427 contract inside 
#> 11 0.224058848983124 0.000466793974113919 -44.7449478078427 contract inside 
#> 12 0.224058848983124 0.000466793974113919 -44.7449478078427 contract inside 
#> 13 0.241585109022692 0.00187733366385477 -42.6722136437239 expand 
#> 14 0.241585109022692 0.00187733366385477 -42.6722136437239 contract inside 
#> 15 0.245063027607108 0.000534793931928738 -42.2744067243519 reflect 
#> 16 0.283725124695804 0.002686531706008 -38.268404787894 expand 
#> 17 0.311666925320418 0.00107428358172823 -35.7092732077242 expand 
#> 18 0.417006504462819 0.00458056624685354 -27.9166533594376 expand 
#> 19 0.552345063361616 0.00310013051989255 -20.5841775642041 expand 
#> 20 0.999710978419122 0.00941693450855305 -6.11082263685136 expand 
#> 21 2.27909404225444 0.00960139832193025 10.0151263830053 expand 
#> 22 5.21705445573082 0.0160006094561096 16.1785096523468 reflect 
#> 23 5.21705445573082 0.0160006094561096 16.1785096523468 contract inside 
#> 24 5.01998823368004 0.0136379132713058 16.2560080596617 contract outside 
#> 25 5.01998823368004 0.0136379132713058 16.2560080596617 contract inside 
#> 26 5.01998823368004 0.0136379132713058 16.2560080596617 contract inside 
#> 27 5.01998823368004 0.0136379132713058 16.2560080596617 contract inside 
#> 28 4.55446419871381 0.0128534405447574 16.2631263946886 contract outside 
#> 29 4.64588795177035 0.0137039945459958 16.2848709659965 contract inside 
#> 30 4.80229839456717 0.0134581943044441 16.2937332389957 contract inside 
#> 31 4.81120051567174 0.0139452912078601 16.2940812796361 contract outside 
#> 32 4.81120051567174 0.0139452912078601 16.2940812796361 contract inside 
#> 33 4.76567827800176 0.0137022618170484 16.2947894537823 contract inside 
#> 34 4.76065312690264 0.0145556889964371 16.2960740367271 expand 
#> 35 4.76065312690264 0.0145556889964371 16.2960740367271 contract inside 
#> 36 4.79022620079006 0.0154864548902051 16.2973233277974 expand 
#> 37 4.75218055478964 0.0169942434321319 16.2997041566054 expand 
#> 38 4.79223074981811 0.0196248205710623 16.3037561005318 expand 
#> 39 4.73629313942515 0.0239978850776008 16.3099408464347 expand 
#> 40 4.78816706953557 0.0315694753860139 16.3225147746448 expand 
#> 41 4.70280849451185 0.0444562799202822 16.3386546709549 expand 
#> 42 4.76300865274624 0.0670993241803681 16.3772660505063 expand 
#> 43 4.62504104143953 0.107370095721328 16.4187629975376 expand 
#> 44 4.67404733819262 0.182915943276343 16.5392504556424 expand 
#> 45 4.43521596124657 0.336230605439469 16.6473252817553 expand 
#> 46 4.41164804858929 0.710768233889407 17.0706425877938 expand 
#> 47 3.98316180283265 2.24020350664308 17.2171748932065 expand 
#> 48 3.98316180283265 2.24020350664308 17.2171748932065 contract inside 
#> 49 4.07346035699536 2.38203955173376 17.4662453522485 reflect 
#> 50 4.07346035699536 2.38203955173376 17.4662453522485 contract inside 
#> 51 4.31161638009175 1.32342613753295 17.5226302331618 reflect 
#> 52 4.16709185895139 2.53686144433066 17.7154198744796 reflect 
#> 53 4.41433359817449 1.3954628618119 17.7185736950648 reflect 
#> 54 4.26424430847961 2.70653812576365 17.9642758896243 reflect 
#> 55 4.26424430847961 2.70653812576365 17.9642758896243 reflect 
#> 56 4.36512011550272 2.89331525059336 18.2120603942574 reflect 
#> 57 4.36512011550272 2.89331525059336 18.2120603942574 contract inside 
#> 58 4.66136773871767 1.88546289068739 18.3961601459291 expand 
#> 59 4.70648903260947 3.21428499862903 18.9854345168409 expand 
#> 60 4.70648903260947 3.21428499862903 18.9854345168409 reflect 
#> 61 5.09425066821065 3.5929320228858 19.7460720492775 reflect 
#> 62 5.09425066821065 3.5929320228858 19.7460720492775 contract inside 
#> 63 5.80189742917594 2.73281194988049 19.891495086665 expand 
#> 64 6.6039561842407 4.74973441119815 21.7148252295986 expand 
#> 65 6.6039561842407 4.74973441119815 21.7148252295986 contract outside 
#> 66 7.98114599108672 6.78056785143906 23.1538234475791 reflect 
#> 67 7.98114599108672 6.78056785143906 23.1538234475791 contract inside 
#> 68 7.98114599108672 6.78056785143906 23.1538234475791 reflect 
#> 69 7.98114599108672 6.78056785143906 23.1538234475791 contract inside 
#> 70 7.98114599108672 6.78056785143906 23.1538234475791 reflect 
#> 71 10.1395538920206 8.45171163553439 23.9139163249958 expand 
#> 72 10.1395538920206 8.45171163553439 23.9139163249958 contract inside 
#> 73 13.8199382485835 15.2298157008651 24.8174728317274 expand 
#> 74 13.8199382485835 15.2298157008651 24.8174728317274 contract inside 
#> 75 13.8199382485835 15.2298157008651 24.8174728317274 contract inside 
#> 76 15.4285564045627 16.0830170757239 25.2469358976012 reflect 
#> 77 15.4285564045627 16.0830170757239 25.2469358976012 contract inside 
#> 78 15.4285564045627 16.0830170757239 25.2469358976012 contract outside 
#> 79 15.4285564045627 16.0830170757239 25.2469358976012 reflect 
#> 80 14.9726890772342 14.9063686208314 25.2528175353109 contract inside 
#> 81 16.2514277486149 16.6500060792838 25.2706584263598 contract inside 
#> 82 15.5076209515346 15.9061391050424 25.2791276917661 contract inside 
#> 83 15.5076209515346 15.9061391050424 25.2791276917661 contract inside 
#> 84 15.8456501667037 16.1784190497506 25.2798613824538 contract inside 
#> 85 15.5412562982367 15.7972482423416 25.2820875172611 contract inside 
#> 86 15.5412562982367 15.7972482423416 25.2820875172611 contract inside 
#> 87 15.5412562982367 15.7972482423416 25.2820875172611 reflect 
#> 88 15.5099956829711 15.8145425305294 25.2822760245013 contract inside 
#> 89 15.5099956829711 15.8145425305294 25.2822760245013 contract inside 
#> 90 15.5014371780236 15.7744980008058 25.2823446966981 contract inside 
#> 91 15.459681324352 15.7418549288367 25.2823501534742 contract inside 
#> 92 15.4952510865558 15.7863042056356 25.2823633180108 contract inside 
#> 93 15.4894333680145 15.7692724247296 25.2823729393679 contract inside 
#> 94 15.4894333680145 15.7692724247296 25.2823729393679 contract outside 
#> 95 15.497160986466 15.784416356866 25.2823738601726 contract inside 
#> 96 15.5010048128956 15.7863183379813 25.2823758721389 contract inside 
#> 97 15.4942566106739 15.7773159990925 25.2823765858442 contract inside 
#> 98 15.4942566106739 15.7773159990925 25.2823765858442 contract outside 
#> 99 15.4942566106739 15.7773159990925 25.2823765858442 contract inside 
#> 100 15.4956590365794 15.7796905716934 25.2823767093344 contract outside 
#> 101 15.4956590365794 15.7796905716934 25.2823767093344 contract inside 
#> 102 15.4956590365794 15.7796905716934 25.2823767093344 contract inside 
#> 103 15.4957155413331 15.779422678147 25.2823767239752 contract inside 
#> 104 15.4957155413331 15.779422678147 25.2823767239752 contract inside 
#> 105 15.4957155413331 15.779422678147 25.2823767239752 contract inside 
#> 106 15.4959900024837 15.7798135534044 25.2823767246439 contract inside 
#> 107 15.4958325900067 15.7796667883301 25.2823767251057 contract inside 
#> Optimization has terminated successfully. 
#> 
#> Maximum likelihood parameter estimates: lambda0: 15.495833, mu0: 15.779667, lambda1: 0.000000, mu1: 0.000000: 
#> Maximum loglikelihood: 25.282377
intGuessLamba <- startingpoint$lambda0
intGuessMu <- startingpoint$mu0
idparsopt <- c(1,2,3,5)
initparsopt <- c(rep(intGuessLamba,3),rep((intGuessLamba/5),1))
idparsfix <- c(0,4)
parsfix <- c(0,0)
tol <- c(1e-02, 1e-03, 1e-03)
maxiter <- 1000 * round((1.25)^length(idparsopt))
optimmethod <- 'simplex'
cond <- 'proper_cond'
root_state_weight <- 'proper_weights'
sampling_fraction <- c(1,1,1)
model<-secsse_ml(
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
verbose = FALSE)
#> Note: you set some transitions as impossible to happen.
# model$ML
# [1] -16.47099
```
