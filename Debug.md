# secsse_hanno_dev

```R
> source("secsse_cla.R")
this tree has:  4126  tips and  4125  internal nodes
Unit: milliseconds
        expr      min       lq     mean   median       uq      max neval  cld
 single thr. 51.45336 51.79697 52.17586 52.14852 52.63214 52.88532    10 a   
   2 threads 28.44547 28.48321 28.92067 28.66617 29.18035 29.95662    10  b  
   4 threads 16.77417 16.87427 18.00470 17.02503 18.24348 22.82410    10   c 
   8 threads 11.23132 11.58550 13.41421 13.44186 14.31095 16.13349    10    d
```

## Breaking changes 

* The `num_steps = NULL` option is gone. This value defaults to `100` in `secsse_loglik_eval` now.
* `eval_cpp` returns a `List` [[output]],[[states]],[[duration]].
* Some superfluous wrapper (`master_xyz`) might still lingering in the code.

## Remaining issues

### misleading comment(s)

This *might* have been an issue for some reasons:

```
#' @note Multithreading might lead to a slightly reduced accuracy
#' (in the order of 1e-10) and is therefore not enabled by default.
#' Please use at your own discretion.
```

Multithreading leads to slightly *different* results due to reordering but
this has nothing to do with accuracy. In fact, the integration itself is
not affected at all.

* `eval_cpp` must return full states because of `collect_node_bars` (info available in the stored matrix).
* Too much `state` data (i.e. `ances` NA states) transfered to C++.
* Inefficient column major matrix memory layout.
* Some superfluous wrapper (`master_xyz`) might still lingering in the code.
* `num_threads` is not passed through all the layers leading to `eval_cpp`. Should be a global(ish) setting anyhow.
* The data layout for the `stored` result is cumbersome on the C++ side. Please double-check in `secsse_eval.cpp`.
* `secsse_sim.h\cpp` could need some tinkering.

## Benching `store`

### `hanno_dev`

```R
> source("secsse_store.R")
this tree has: 569 tips and 568 internal nodes, num_steps = 10 
Unit: milliseconds
        expr       min        lq     mean    median        uq      max neval
 single thr. 51.039688 51.493067 51.71075 51.587962 51.708617 53.00567    10
   2 threads 26.333897 26.602829 26.78888 26.791953 26.890987 27.65554    10
   4 threads 14.000106 14.289136 14.62157 14.485783 15.145808 15.55330    10
   8 threads  8.605093  8.821089 13.82669  9.095943  9.291371 55.87885    10

> source("secsse_store.R")
this tree has: 569 tips and 568 internal nodes, num_steps = 100 
Unit: milliseconds
        expr       min        lq      mean    median        uq      max neval
 single thr. 434.63324 436.94616 442.20108 437.84397 439.50956 484.0255    10
   2 threads 217.26720 218.61650 223.89733 219.31238 219.44311 268.6833    10
   4 threads 112.11049 114.06196 114.82934 115.24166 115.79207 116.1735    10
   8 threads  60.84545  63.67806  69.89058  64.82028  66.40901 115.8395    10

> source("secsse_store.R")
this tree has: 957 tips and 956 internal nodes, num_steps = 10 
Unit: milliseconds
        expr      min       lq     mean   median       uq      max neval
 single thr. 83.63308 84.24488 84.60058 84.52680 84.93523 85.68493    10   
   2 threads 42.82260 42.94673 43.11096 43.10289 43.32173 43.43227    10  
   4 threads 22.73277 23.02302 23.59615 23.50217 24.20230 24.50711    10 
   8 threads 12.93618 13.02415 14.55027 14.21917 15.99571 16.84265    10

> source("secsse_store.R")
this tree has: 957 tips and 956 internal nodes, num_steps = 100 
Unit: milliseconds
        expr      min       lq     mean   median       uq      max neval
 single thr. 724.5631 732.3597 733.8956 734.4381 736.2300 738.4202    10   
   2 threads 365.2519 366.5167 373.6500 367.9019 371.7572 417.8656    10  
   4 threads 191.7030 191.8802 192.7503 192.5686 192.9304 195.2887    10 
   8 threads 102.3770 102.8724 104.4372 103.5731 104.8228 109.2572    10

# a big one?
> source("secsse_store.R")
this tree has: 8531 tips and 8530 internal nodes, num_steps = 100 
Unit: milliseconds
        expr       min        lq      mean    median        uq       max neval
 single thr. 6569.2933 6577.9312 6600.2636 6583.8110 6625.1909 6662.5209    10
   2 threads 3288.8340 3293.3418 3310.7211 3311.7254 3321.3449 3341.6627    10
   4 threads 1707.3803 1710.9579 1717.7718 1717.9355 1723.1349 1728.6820    10
   8 threads  897.1503  905.0342  919.0528  909.4819  916.3939  985.9016    10
```

### `develop` (note the units)

```R
source("secsse_store.R")
this tree has: 569 tips and 568 internal nodes, num_steps = 10 
Unit: milliseconds
        expr      min       lq     mean   median       uq      max neval
 single thr. 230.2166 238.1931 240.7193 239.9897 246.0346 248.2763    10

> source("secsse_store.R")
this tree has: 569 tips and 568 internal nodes, num_steps = 100 
Unit: seconds
        expr      min       lq     mean   median       uq      max neval
 single thr. 2.068572 2.072582 2.093156 2.082578 2.110746 2.140897    10

> source("secsse_store.R")
this tree has: 957 tips and 956 internal nodes, num_steps = 10 
Unit: milliseconds
        expr      min       lq     mean   median       uq      max neval
 single thr. 391.0567 391.5107 394.3454 394.1759 395.6691 399.1920    10

> source("secsse_store.R")
this tree has: 957 tips and 956 internal nodes, num_steps = 100 
Unit: seconds
        expr      min       lq     mean   median       uq      max neval
 single thr. 3.447375 3.484134 3.502839 3.511399 3.527202 3.534494    10

# a big one?
source("secsse_store.R")
this tree has: 8531 tips and 8530 internal nodes, num_steps = 100 
Unit: seconds
        expr      min       lq     mean   median       uq      max neval
 single thr. 30.77723 30.86942 31.02351 30.94320 31.18533 31.50008    10
```

