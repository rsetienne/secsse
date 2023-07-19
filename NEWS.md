# secsse 3.0.0

Version 3.0.0 is expected to arrive to CRAN in the second half of 2023. It 
extends the C++ code base used for the standard likelihood to the cla
likelihood, harnessing the same computation improvement. 

## Breaking changes



## Major changes

## Minor changes
* Added a `NEWS.md` file to track changes to the package.

# 2.6.0

* `fill`Version 2.6.0 appeared on CRAN in July 2023, and introduced many functions 
suited to prepare the parameter structure for secsse. It also introduced a new
C++ code base for the standard likelihood, making smarter use of
parallelization, this marks another 10-fold increase in speed.

### 2.5.0
Version 2.5.0 appeared in 2021 on GitHub and was published in May 2023 on CRAN.
Version 2.5.0 marks the first version using C++ to perform the integration,
and it used tbb (from the RcppParallel package) to perform multithreading. This
marks a ten fold increase in speed over previous versions.

# 2.0.0
Version 2.0.0 appeared in June of 2019 on CRAN and extended the package with the
cla framework, e.g. including state shifts during speciation / asymmetric 
inheritance during speciation. 

# 1.0.0
The first version of secsse appeared in January of 2019 on CRAN. It used the
package deSolve to solve all integrations, and could switch between either using
a fully R based evaluation, or use FORTRAN to speed up calculations.
Furthermore, using the foreach package, within-R parallelization was
implemented. However, parallelization only situationally improved computation
times, and generally, computation was relatively slow.