# Breaking changes 

* The `num_threads = NULL` is gone. This value defaults to `100` in `secsse_loglik_eval` now.
* `eval_cpp` returns a `List` [[output]],[[states]],[[duration]].
* Some superfluous wrapper (`master_xyz`) might still lingering in the code.

# Remaining issues

* Some superfluous wrapper (`master_xyz`) might still lingering in the code.
* `num_threads` is not passed through all the layers leading to `eval_cpp`. Should be a global(ish) setting anyhow.
* `secsse_sim.h\cpp` could need some tinkering.

# 'store' bench `hanno_dev`

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
        expr      min       lq     mean   median       uq      max neval  cld
 single thr. 83.63308 84.24488 84.60058 84.52680 84.93523 85.68493    10 a   
   2 threads 42.82260 42.94673 43.11096 43.10289 43.32173 43.43227    10  b  
   4 threads 22.73277 23.02302 23.59615 23.50217 24.20230 24.50711    10   c 
   8 threads 12.93618 13.02415 14.55027 14.21917 15.99571 16.84265    10    d

> source("secsse_store.R")
this tree has: 957 tips and 956 internal nodes, num_steps = 100 
Unit: milliseconds
        expr      min       lq     mean   median       uq      max neval  cld
 single thr. 724.5631 732.3597 733.8956 734.4381 736.2300 738.4202    10 a   
   2 threads 365.2519 366.5167 373.6500 367.9019 371.7572 417.8656    10  b  
   4 threads 191.7030 191.8802 192.7503 192.5686 192.9304 195.2887    10   c 
   8 threads 102.3770 102.8724 104.4372 103.5731 104.8228 109.2572    10    d
> 
```

# 'store' bench `develop` (note the units)

```R
source("secsse_store.R")
this tree has: 569 tips and 568 internal nodes, num_steps = 10 
Unit: milliseconds
        expr      min       lq     mean   median       uq      max neval cld
 single thr. 230.2166 238.1931 240.7193 239.9897 246.0346 248.2763    10   a
   2 threads 235.5266 236.4017 239.4577 238.2896 243.1415 246.1469    10   a
   4 threads 233.1688 233.6139 239.7728 239.3966 246.1307 248.3981    10   a
   8 threads 232.2317 233.5301 237.3019 235.8108 240.7786 244.1390    10   a
There were 30 warnings (use warnings() to see them)

> source("secsse_store.R")
this tree has: 569 tips and 568 internal nodes, num_steps = 100 
Unit: seconds
        expr      min       lq     mean   median       uq      max neval cld
 single thr. 2.068572 2.072582 2.093156 2.082578 2.110746 2.140897    10   a
   2 threads 2.055411 2.067366 2.081735 2.078422 2.084541 2.119112    10   a
   4 threads 2.066754 2.070876 2.075975 2.075220 2.082172 2.085049    10   a
   8 threads 2.064539 2.066838 2.078086 2.078250 2.086002 2.097862    10   a
There were 30 warnings (use warnings() to see them)

> source("secsse_store.R")
this tree has: 957 tips and 956 internal nodes, num_steps = 10 
Unit: milliseconds
        expr      min       lq     mean   median       uq      max neval cld
 single thr. 391.0567 391.5107 394.3454 394.1759 395.6691 399.1920    10   a
   2 threads 389.8366 391.5455 394.2471 394.1692 396.4653 399.5887    10   a
   4 threads 390.8509 391.5043 392.8117 392.7376 394.4531 395.1885    10   a
   8 threads 393.2036 393.5428 395.1897 394.2929 397.2716 399.3219    10   a
There were 30 warnings (use warnings() to see them)

> source("secsse_store.R")
this tree has: 957 tips and 956 internal nodes, num_steps = 100 
Unit: seconds
        expr      min       lq     mean   median       uq      max neval cld
 single thr. 3.447375 3.484134 3.502839 3.511399 3.527202 3.534494    10   a
   2 threads 3.473739 3.489593 3.513760 3.505985 3.520598 3.596552    10   a
   4 threads 3.462059 3.501752 3.517745 3.518999 3.538560 3.567041    10   a
   8 threads 3.475357 3.481082 3.515387 3.507610 3.532091 3.581478    10   a
There were 30 warnings (use warnings() to see them)
```

