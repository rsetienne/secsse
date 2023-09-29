# secsse 3.0.1

Version 3.0.1 patches some inaccuracies in simulation functions, and
deprecates expand_q_matrix, as this was making some incorrect
assumptions.

## Breaking changes

-   The function expand_q_matrix is now deprecated, please use
    q_doubletrans

## Minor changes

-   Added conditioning of simulation of complete trees on complete tree
    size
-   Fixed some issues with setting `pool_init_states`
-   Added option to return a histogram of simulated tree sizes, across
    all successful and failed trees.

# 3.0.0

Version 3.0.0 extends the C++ code base used for the standard likelihood
to the "cla\_" likelihood, harnessing the same computation improvement.

## Breaking changes

-   Function name changes:
    -   `create_lambda_matrices()` is now called `create_lambda_list()`
    -   `create_transition_matrix()` is now called `create_q_matrix()`
    -   `create_mus()` is now called `create_mu_vector()`
    -   `create_default_q_list()` is now called
        `create_default_shift_matrix()`
    -   `create_default_lambda_list ()` is now called
        `create_default_lambda_transition_matrix()`
    -   `create_default_q_list()` is now called
        `create_default_shift_matrix()`
-   Package data files renamed:
    -   `phylo_Vign` is now called `phylo_vignette`
    -   `traitinfo` is now called `traits`
    -   `phy` is now called `example_phy_GeoSSE`
-   `plot_state_exact()` argument `steps` renamed to `num_steps` and
    argument `focal_tree` renamed to `phy` for consistency with other
    functions.

## Major changes

-   Vastly improve the computational speed of "cla\_" likelihood
    calculation.
-   Optimization of parallelization resulting in better scaling with
    more threads and faster run time for standard secsse and cla_secsse
    likelihood calculations.

## Minor changes

-   Added a `NEWS.md` file to track changes to the package.
-   Documentation reworked into `default_params_doc()`.
-   Several documentation formatting improvements and linking.
    Documentation now follows and allows for roxygen2 markdown.
-   A new vignette:
    -   *Using secsse with complete phylogenies (with extinction)*
        `vignette("complete_tree", package = "secsse")`
-   A new [pkgdown
    website](https://rsetienne.github.io/secsse/index.html)!
    -   It contains all the documentation and vignettes of the package,
        along with additional interesting information like the *Secsse
        versions* article with details on performance and the
        development history of secsse.
-   Revise, combine and simplify the *Using SecSSE ML search* and
    *Setting up a secsse analysis* into the *Starting secsse* vignette
    `vignette("starting_secsse", package = "secsse")`.
-   `secsse_sim()` argument `conditioning` now defaults to
    `"obs_states"` from `"none"`.
-   No longer Import package 'stringr' and Suggest package 'testit'.
-   New organisation of code in .R, .cpp and .h files. (Developer side).
-   Start archiving in Zenodo, with new .zenodo.json metadata file.

## Bug fixes

-   `secsse_sim()` fix bug causing error when simulating trees with
    extinct species.

# 2.6.0

## Major changes

-   C++ code base for the standard likelihood, making smarter use of
    parallelization, this marks another 10-fold increase in speed.

## Minor changes

-   Add a number of helper functions: `fill_in()`,
    `create_default_q_list()`, `create_default_transition_list()`,
    `create_mus()`
-   Implemented necessary changes to comply with CRAN clang16 build and
    solve issue with the boost odeint library uninitialized variable
    (see <https://github.com/boostorg/odeint/issues/59> and more details
    at <https://github.com/rsetienne/DAISIE/pull/158>)
-   Updated Copyright license to the Boost Software License, Version 1.0
    for included C++ code (R code remains GPL\>=3).

## Bug fixes

-   Fix memory leaks

# 2.5.0

Version 2.5.0 appeared in 2021 on GitHub and was published in May 2023
on CRAN. Version 2.5.0 marks the first version using C++ to perform the
integration, and it used tbb (from the RcppParallel package) to perform
multithreading. This marks a ten fold increase in speed over previous
versions. Secondly, 2.5.0 introduces the function `secsse_sim()` to
simulate a diversification process using the (cla) secsse framework.
Lastly, in version 2.5.0 functions were added to allow visualisation of
inferred rates of speciation across the tree (e.g. `plot_state_exact()`
and `secsse_loglik_eval()`).

# 2.0.0

Version 2.0.0 appeared in June of 2019 on CRAN and extended the package
with the cla framework, e.g. including state shifts during speciation /
asymmetric inheritance during speciation.

# 1.0.0

The first version of secsse appeared in January of 2019 on CRAN. It used
the package deSolve to solve all integrations, and could switch between
either using a fully R based evaluation, or use FORTRAN to speed up
calculations. Furthermore, using the foreach package, within-R
parallelization was implemented. However, parallelization only
situationally improved computation times, and generally, computation was
relatively slow.
