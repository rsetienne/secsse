# Event times of a (possibly non-ultrametric) phylogenetic tree

Times at which speciation or extinction occurs

## Usage

``` r
event_times(phy)
```

## Arguments

- phy:

  phylogenetic tree of class phylo, without polytomies, rooted and with
  branch lengths. Need not be ultrametric.

## Value

times at which speciation or extinction happens.

## Note

This script has been modified from BAMMtools' internal function
NU.branching.times
