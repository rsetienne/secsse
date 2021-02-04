
# SecSSE: Several Examined and Concealed States-Dependent Speciation and Extinction

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/secsse)](https://CRAN.R-project.org/package=secsse)
[![](http://cranlogs.r-pkg.org/badges/grand-total/secsse)]( https://CRAN.R-project.org/package=secsse)
[![](http://cranlogs.r-pkg.org/badges/secsse)](https://CRAN.R-project.org/package=secsse)
<!-- badges: end -->

Branch|[![GitHub Actions logo](pics/github_actions_logo.png)](https://github.com/features/actions)|[![Codecov logo](pics/Codecov.png)](https://www.codecov.io)
--------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------
`master`|[![R build status](https://github.com/rsetienne/secsse/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/rsetienne/secsse/actions)|[![codecov.io](https://codecov.io/github/rsetienne/secsse/coverage.svg?branch=master)](https://codecov.io/github/rsetienne/secsse/branch/master)
`develop`|[![R build status](https://github.com/rsetienne/secsse/workflows/R-CMD-check/badge.svg?branch=develop)](https://github.com/rsetienne/secsse/actions)|[![codecov.io](https://codecov.io/github/rsetienne/secsse/coverage.svg?branch=develop)](https://codecov.io/github/rsetienne/secsse/branch/develop)

## What is SecSSE?
SecSSE is an R package designed for multistate data sets under a concealed state and speciation (`hisse`) framework. In this sense, it is parallel to the 'MuSSE' functionality implemented in `diversitree`, but it accounts for finding possible spurious relationships between traits and diversification rates ("false positives", Rabosky & Goldberg 2015) by testing against a "hidden trait" (Beaulieu et al. 2013), which is responsible for more variation in diversification rates than the trait being investigated. 

## Installing secsse

The `secsse` package has a stable version on [CRAN](https://CRAN.R-project.org/package=secsse) and a development version on GitHub.

### From CRAN

From within R, do:

``` r
install.packages("secsse")
```

### From GitHub

Install secsse from this GitHub repository by running:

``` r
install.packages("remotes")
remotes::install_github("rsetienne/secsse")
```

## Using secsse as a package dependency

### From CRAN

To your DESCRIPTION file, add `secsse` as any normal package.

If your package directly uses `secsse`:

```
Imports:
  secsse
```

If your package uses `secsse` in its peripherals (e.g. vignettes and tests):

```
Suggests:
  secsse
```

### From GitHub

```
Remotes:
  rsetienne/secsse
```

## Cite
If you use `secsse` in your publications, please cite:

* Herrera-Alsina, Leonel et al. “Detecting the Dependence of Diversification on Multiple Traits from Phylogenetic Trees and Trait Data.” Systematic biology vol. 68,2 (2019): 317-328. doi:10.1093/sysbio/syy057

## References
* Beaulieu, Jeremy M et al. “Identifying hidden rate changes in the evolution of a binary morphological character: the evolution of plant habit in campanulid angiosperms.” Systematic biology vol. 62,5 (2013): 725-37. doi:10.1093/sysbio/syt034

* Beaulieu, Jeremy M, and Brian C O'Meara. “Detecting Hidden Diversification Shifts in Models of Trait-Dependent Speciation and Extinction.” Systematic biology vol. 65,4 (2016): 583-601. doi:10.1093/sysbio/syw022

* Rabosky, Daniel L., and Emma E. Goldberg. "Model inadequacy and mistaken inferences of trait-dependent speciation." Systematic biology 64.2 (2015): 340-355. doi:10.1093/sysbio/syu131