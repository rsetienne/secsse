# Likelihood for SecSSE model Logikelihood calculation for the SecSSE model given a set of parameters and data, returning also the likelihoods along the branches

Likelihood for SecSSE model Logikelihood calculation for the SecSSE
model given a set of parameters and data, returning also the likelihoods
along the branches

## Usage

``` r
secsse_loglik_eval(
  parameter,
  phy,
  traits,
  num_concealed_states,
  cond = "proper_cond",
  root_state_weight = "proper_weights",
  sampling_fraction,
  setting_calculation = NULL,
  loglik_penalty = 0,
  is_complete_tree = FALSE,
  num_threads = 1,
  atol = 1e-08,
  rtol = 1e-07,
  method = "odeint::runge_kutta_cash_karp54",
  num_steps = 100
)
```

## Arguments

- parameter:

  list where first vector represents lambdas, the second mus and the
  third transition rates.

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

- setting_calculation:

  argument used internally to speed up calculation. This should be left
  blank (default : `setting_calculation = NULL`).

- loglik_penalty:

  the size of the penalty for all parameters; default is 0 (no penalty).

- is_complete_tree:

  logical specifying whether or not a tree with all its extinct species
  is provided. If set to `TRUE`, it also assumes that all *all* extinct
  lineages are present on the tree. Defaults to `FALSE`.

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

- num_steps:

  number of substeps to show intermediate likelihoods along a branch.

## Value

A list containing: "output", observed states along evaluated time points
along all branches, used for plotting. "states" all ancestral states on
the nodes and "duration", indicating the time taken for the total
evaluation

## Examples

``` r
set.seed(5)
phy <- ape::rphylo(n = 4, birth = 1, death = 0)
traits <- c(0, 1, 1, 0)
params <- secsse::id_paramPos(c(0, 1), 2)
params[[1]][] <- c(0.2, 0.2, 0.1, 0.1)
params[[2]][] <- 0.0
params[[3]][, ] <- 0.1
diag(params[[3]]) <- NA

secsse_loglik_eval(parameter = params,
                   phy = phy,
                   traits = traits,
                   num_concealed_states = 2,
                   sampling_fraction = c(1, 1),
                   num_steps = 10)
#> Warning: Deduced names and order of used states to be: 0, 1
#> if this is incorrect, consider passing states as matching numeric 
#>   ordering, e.g. 1 for the first state, 2 for the second etc.
#> $output
#>       [,1] [,2]        [,3] [,4] [,5] [,6] [,7]       [,8]       [,9]
#>  [1,]    7    1 0.000000000    0    0    0    0 1.00000000 0.00000000
#>  [2,]    7    1 0.049700250    0    0    0    0 0.98037680 0.00975672
#>  [3,]    7    1 0.099400499    0    0    0    0 0.96125773 0.01915438
#>  [4,]    7    1 0.149100749    0    0    0    0 0.94262856 0.02820413
#>  [5,]    7    1 0.198800998    0    0    0    0 0.92447549 0.03691677
#>  [6,]    7    1 0.248501248    0    0    0    0 0.90678511 0.04530280
#>  [7,]    7    1 0.298201497    0    0    0    0 0.88954442 0.05337238
#>  [8,]    7    1 0.347901747    0    0    0    0 0.87274075 0.06113538
#>  [9,]    7    1 0.397601996    0    0    0    0 0.85636185 0.06860137
#> [10,]    7    1 0.447302246    0    0    0    0 0.84039578 0.07577964
#> [11,]    7    1 0.497002495    0    0    0    0 0.82483097 0.08267919
#> [12,]    7    3 0.000000000    0    0    0    0 0.00000000 1.00000000
#> [13,]    7    3 0.049700250    0    0    0    0 0.00975672 0.98037680
#> [14,]    7    3 0.099400499    0    0    0    0 0.01915438 0.96125773
#> [15,]    7    3 0.149100749    0    0    0    0 0.02820413 0.94262856
#> [16,]    7    3 0.198800998    0    0    0    0 0.03691677 0.92447549
#> [17,]    7    3 0.248501248    0    0    0    0 0.04530280 0.90678511
#> [18,]    7    3 0.298201497    0    0    0    0 0.05337238 0.88954442
#> [19,]    7    3 0.347901747    0    0    0    0 0.06113538 0.87274075
#> [20,]    7    3 0.397601996    0    0    0    0 0.06860137 0.85636185
#> [21,]    7    3 0.447302246    0    0    0    0 0.07577964 0.84039578
#> [22,]    7    3 0.497002495    0    0    0    0 0.08267919 0.82483097
#> [23,]    6    2 0.000000000    0    0    0    0 0.00000000 1.00000000
#> [24,]    6    2 0.077391750    0    0    0    0 0.01503639 0.96966295
#> [25,]    6    2 0.154783500    0    0    0    0 0.02921721 0.94052908
#> [26,]    6    2 0.232175250    0    0    0    0 0.04258347 0.91254599
#> [27,]    6    2 0.309566999    0    0    0    0 0.05517431 0.88566364
#> [28,]    6    2 0.386958749    0    0    0    0 0.06702706 0.85983423
#> [29,]    6    2 0.464350499    0    0    0    0 0.07817729 0.83501208
#> [30,]    6    2 0.541742249    0    0    0    0 0.08865896 0.81115355
#> [31,]    6    2 0.619133999    0    0    0    0 0.09850443 0.78821695
#> [32,]    6    2 0.696525749    0    0    0    0 0.10774457 0.76616246
#> [33,]    6    2 0.773917499    0    0    0    0 0.11640881 0.74495202
#> [34,]    6    7 0.000000000    0    0    0    0 0.32533265 0.32533265
#> [35,]    6    7 0.027691500    0    0    0    0 0.32271192 0.32271192
#> [36,]    6    7 0.055383001    0    0    0    0 0.32012191 0.32012191
#> [37,]    6    7 0.083074501    0    0    0    0 0.31756216 0.31756216
#> [38,]    6    7 0.110766001    0    0    0    0 0.31503227 0.31503227
#> [39,]    6    7 0.138457502    0    0    0    0 0.31253180 0.31253180
#> [40,]    6    7 0.166149002    0    0    0    0 0.31006034 0.31006034
#> [41,]    6    7 0.193840502    0    0    0    0 0.30761750 0.30761750
#> [42,]    6    7 0.221532003    0    0    0    0 0.30520286 0.30520286
#> [43,]    6    7 0.249223503    0    0    0    0 0.30281603 0.30281603
#> [44,]    6    7 0.276915003    0    0    0    0 0.30045663 0.30045663
#> [45,]    5    4 0.000000000    0    0    0    0 1.00000000 0.00000000
#> [46,]    5    4 0.080187748    0    0    0    0 0.96858978 0.01556336
#> [47,]    5    4 0.160375497    0    0    0    0 0.93846918 0.03020982
#> [48,]    5    4 0.240563245    0    0    0    0 0.90958007 0.04398490
#> [49,]    5    4 0.320750993    0    0    0    0 0.88186701 0.05693194
#> [50,]    5    4 0.400938741    0    0    0    0 0.85527713 0.06909221
#> [51,]    5    4 0.481126490    0    0    0    0 0.82975999 0.08050499
#> [52,]    5    4 0.561314238    0    0    0    0 0.80526751 0.09120769
#> [53,]    5    4 0.641501986    0    0    0    0 0.78175378 0.10123591
#> [54,]    5    4 0.721689735    0    0    0    0 0.75917506 0.11062355
#> [55,]    5    4 0.801877483    0    0    0    0 0.73748961 0.11940288
#> [56,]    5    6 0.000000000    0    0    0    0 0.10273725 0.65746159
#> [57,]    5    6 0.002795998    0    0    0    0 0.10284434 0.65663886
#> [58,]    5    6 0.005591997    0    0    0    0 0.10295111 0.65581736
#> [59,]    5    6 0.008387995    0    0    0    0 0.10305757 0.65499711
#> [60,]    5    6 0.011183994    0    0    0    0 0.10316371 0.65417809
#> [61,]    5    6 0.013979992    0    0    0    0 0.10326953 0.65336031
#> [62,]    5    6 0.016775991    0    0    0    0 0.10337505 0.65254376
#> [63,]    5    6 0.019571989    0    0    0    0 0.10348024 0.65172845
#> [64,]    5    6 0.022367987    0    0    0    0 0.10358513 0.65091436
#> [65,]    5    6 0.025163986    0    0    0    0 0.10368970 0.65010151
#> [66,]    5    6 0.027959984    0    0    0    0 0.10379396 0.64928987
#>             [,10]       [,11]
#>  [1,] 1.000000000 0.000000000
#>  [2,] 0.985237013 0.009780916
#>  [3,] 0.970763781 0.019249183
#>  [4,] 0.956573789 0.028413071
#>  [5,] 0.942660674 0.037280642
#>  [6,] 0.929018225 0.045859752
#>  [7,] 0.915640376 0.054158058
#>  [8,] 0.902521204 0.062183024
#>  [9,] 0.889654925 0.069941924
#> [10,] 0.877035893 0.077441847
#> [11,] 0.864658592 0.084689704
#> [12,] 0.000000000 1.000000000
#> [13,] 0.009780916 0.985237013
#> [14,] 0.019249183 0.970763781
#> [15,] 0.028413071 0.956573789
#> [16,] 0.037280642 0.942660674
#> [17,] 0.045859752 0.929018225
#> [18,] 0.054158058 0.915640376
#> [19,] 0.062183024 0.902521204
#> [20,] 0.069941924 0.889654925
#> [21,] 0.077441847 0.877035893
#> [22,] 0.084689704 0.864658592
#> [23,] 0.000000000 1.000000000
#> [24,] 0.015094388 0.977137634
#> [25,] 0.029441848 0.954969055
#> [26,] 0.043072947 0.933470146
#> [27,] 0.056017056 0.912617668
#> [28,] 0.068302395 0.892389228
#> [29,] 0.079956077 0.872763249
#> [30,] 0.091004153 0.853718938
#> [31,] 0.101471650 0.835236258
#> [32,] 0.111382616 0.817295901
#> [33,] 0.120760151 0.799879260
#> [34,] 0.174667355 0.174667355
#> [35,] 0.175009400 0.175009400
#> [36,] 0.175334246 0.175334246
#> [37,] 0.175642204 0.175642204
#> [38,] 0.175933579 0.175933579
#> [39,] 0.176208671 0.176208671
#> [40,] 0.176467777 0.176467777
#> [41,] 0.176711187 0.176711187
#> [42,] 0.176939189 0.176939189
#> [43,] 0.177152064 0.177152064
#> [44,] 0.177350089 0.177350089
#> [45,] 1.000000000 0.000000000
#> [46,] 0.976324799 0.015625549
#> [47,] 0.953393474 0.030450425
#> [48,] 0.931179252 0.044508557
#> [49,] 0.909656370 0.057832498
#> [50,] 0.888800036 0.070453481
#> [51,] 0.868586391 0.082401469
#> [52,] 0.848992477 0.093705212
#> [53,] 0.829996198 0.104392289
#> [54,] 0.811576291 0.114489160
#> [55,] 0.793712291 0.124021207
#> [56,] 0.031454687 0.208346470
#> [57,] 0.031690079 0.208334741
#> [58,] 0.031925006 0.208322891
#> [59,] 0.032159467 0.208310921
#> [60,] 0.032393463 0.208298829
#> [61,] 0.032626996 0.208286618
#> [62,] 0.032860065 0.208274286
#> [63,] 0.033092671 0.208261835
#> [64,] 0.033324816 0.208249264
#> [65,] 0.033556499 0.208236575
#> [66,] 0.033787722 0.208223766
#> 
#> $states
#>      [,1] [,2] [,3] [,4]      [,5]      [,6]       [,7]       [,8] [,9] [,10]
#> [1,]    0    0    0    0 1.0000000 0.0000000 1.00000000 0.00000000    1     1
#> [2,]    0    0    0    0 0.0000000 1.0000000 0.00000000 1.00000000    1     1
#> [3,]    0    0    0    0 0.0000000 1.0000000 0.00000000 1.00000000    1     1
#> [4,]    0    0    0    0 1.0000000 0.0000000 1.00000000 0.00000000    1     1
#> [5,]    0    0    0    0 0.4243298 0.4297629 0.07433059 0.07157672    0     0
#> [6,]    0    0    0    0 0.1027373 0.6574616 0.03145469 0.20834647    0     0
#> [7,]    0    0    0    0 0.3253326 0.3253326 0.17466735 0.17466735    0     0
#>      [,11] [,12]
#> [1,]     1     1
#> [2,]     1     1
#> [3,]     1     1
#> [4,]     1     1
#> [5,]     0     0
#> [6,]     0     0
#> [7,]     0     0
#> 
#> $duration
#> [1] 0.00013379
#> 
```
