# Basic Qmatrix Sets a Q matrix where double transitions are not allowed

This function expands the Q_matrix, but it does so assuming that the
number of concealed traits is equal to the number of examined traits, if
you have a different number, you should consider looking at the function
[`expand_q_matrix()`](https://rsetienne.github.io/secsse/reference/expand_q_matrix.md).

## Usage

``` r
q_doubletrans(traits, masterBlock, diff.conceal)
```

## Arguments

- traits:

  vector with trait states for each tip in the phylogeny. The order of
  the states must be the same as the tree tips. For help, see
  [`vignette("starting_secsse", package = "secsse")`](https://rsetienne.github.io/secsse/articles/starting_secsse.md).
  When providing a `multiPhylo` set of multiple phylognies, traits
  should be a list where each entry in the list corresponds to the
  matching phylogeny on that position.

- masterBlock:

  matrix of transitions among only examined states, `NA` in the main
  diagonal, used to build the full transition rates matrix.

- diff.conceal:

  Boolean stating if the concealed states should be different. E.g. that
  the transition rates for the concealed states are different from the
  transition rates for the examined states. Normally it should be
  `FALSE` in order to avoid having a huge number of parameters.

## Value

Q matrix that includes both examined and concealed states, it should be
declared as the third element of idparslist.

## Examples

``` r
traits <- sample(c(0,1,2), 45,replace = TRUE) #get some traits
# For a three-state trait
masterBlock <- matrix(99,ncol = 3,nrow = 3,byrow = TRUE)
diag(masterBlock) <- NA
masterBlock[1,2] <- 6
masterBlock[1,3] <- 7
masterBlock[2,1] <- 8
masterBlock[2,3] <- 9
masterBlock[3,1] <- 10
masterBlock[3,2] <- 11
myQ <- q_doubletrans(traits,masterBlock,diff.conceal = FALSE)
# now, it can replace the Q matrix from id_paramPos
num_concealed_states <- 3
param_posit <- id_paramPos(traits,num_concealed_states)
param_posit[[3]] <- myQ
```
