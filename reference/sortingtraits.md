# Data checking and trait sorting In preparation for likelihood calculation, it orders trait data according the tree tips

Data checking and trait sorting In preparation for likelihood
calculation, it orders trait data according the tree tips

## Usage

``` r
sortingtraits(trait_info, phy)
```

## Arguments

- trait_info:

  data frame where first column has species ids and the second one is
  the trait associated information.

- phy:

  phylogenetic tree of class `phylo`, rooted and with branch lengths.
  Alternatively, multiple phylogenetic trees can be provided as the
  `multiPhylo` class.

## Value

Vector of traits

## Examples

``` r
# Some data we have prepared
data(traits)
data('phylo_vignette')
traits <- sortingtraits(traits, phylo_vignette)
```
