# function to visualize the structure of the idparslist

function to visualize the structure of the idparslist

## Usage

``` r
plot_idparslist(idparslist, state_names)
```

## Arguments

- idparslist:

  idparslist setup list

- state_names:

  names of all states (including concealed states)

## Value

list with two ggplot objects: plot_qmat which visualizes the q_matrix,
and plot_lambda which visualizes the lambda matrices.

## Examples

``` r
idparslist <- list()
focal_matrix <-
secsse::create_default_lambda_transition_matrix(state_names = c("1", "2"),
                                                model = "CR")
idparslist[[1]] <- 
  secsse::create_lambda_list(state_names = c("1", "2"),
                             num_concealed_states = 2,
                             transition_matrix = focal_matrix,
                             model = "CR")
idparslist[[2]] <- secsse::create_mu_vector(state_names = c("1", "2"),
                                            num_concealed_states = 2,
                                            model = "CR",
                                            lambda_list = idparslist[[1]])
shift_mat <- secsse::create_default_shift_matrix(state_names = c("1", "2"),
                                                 num_concealed_states = 2,
                                                 mu_vector = idparslist[[2]])
idparslist[[3]] <- secsse::create_q_matrix(state_names = c("1", "2"),
                                           num_concealed_states = 2,
                                           shift_matrix = shift_mat,
                                           diff.conceal = FALSE)
p <- plot_idparslist(idparslist, state_names = names(idparslist[[1]]))
p$plot_lambda

p$plot_qmat
```
