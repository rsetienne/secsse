
state_names <- c("H", "W")

lambda_trans <- secsse::create_default_lambda_transition_matrix(model = "ETD",
                                                                state_names = state_names)
lambdas <- secsse::create_lambda_list(transition_matrix = lambda_trans,
                                      state_names = state_names,
                                      model = "ETD")
mus <- secsse::create_mu_vector(state_names = state_names,
                                num_concealed_states = 2,
                                lambda_list = lambdas)
shiftmat <- secsse::create_default_shift_matrix(mu_vector = mus,
                                                state_names = state_names)
qmat <- secsse::create_q_matrix(state_names = state_names,
                                num_concealed_states = 2,
                                shift_matrix = shiftmat)
idparslist <- list()
idparslist[[1]] <- lambdas
idparslist[[2]] <- mus
idparslist[[3]] <- qmat

all_trees <- readRDS("/Users/thijsjanzen/Library/CloudStorage/Dropbox/projects/multi_trees/multi_trees_woody.rds")

focal_tree <- all_trees[[155]]
focal_tree$root.edge <- NULL
focal_tree$edge.length <- 1

set.seed(42)
#initpars <- c(0.3, 0.0, 0, 0)
initpars <- runif(4)
cat(initpars, "\n")
filled_in_pars <- list()
filled_in_pars[[1]] <- secsse::fill_in(idparslist[[1]], initpars)
filled_in_pars[[2]] <- secsse::fill_in(idparslist[[2]], initpars)
filled_in_pars[[3]] <- secsse::fill_in(idparslist[[3]], initpars)

found <- c()
for (bl in seq(1, 50, by = 0.5)) {
  focal_tree$edge.length <- bl
  for (bup in c(0, 1, 2)) {
    l2 <- secsse::secsse_single_branch_loglik(parameter = filled_in_pars,
                                              phy = focal_tree,
                                              traits = "H",
                                              num_concealed_states = 2,
                                              sampling_fraction = c(1, 1),
                                              take_into_account_root_edge = FALSE,
                                              display_warning = FALSE,
                                              break_up = bup)
    to_add <- c(bl, l2$loglik, bup)
    #cat(to_add, "\n")
    found <- rbind(found, to_add)
  }
}
require(tidyverse)

colnames(found) <- c("branch_length", "ll", "using_breaking")
found <- as_tibble(found)
found %>%
  ggplot(aes(x = branch_length, y = ll, col = as.factor(using_breaking))) +
    geom_line()

run_estim <- function(bup) {
  secsse::secsse_single_branch_loglik(parameter = filled_in_pars,
                                      phy = focal_tree,
                                      traits = "H",
                                      num_concealed_states = 2,
                                      sampling_fraction = c(1, 1),
                                      take_into_account_root_edge = FALSE,
                                      display_warning = FALSE,
                                      break_up = bup)
}

run_estim(0)
run_estim(1)
run_estim(2)

vx <- microbenchmark::microbenchmark(run_estim(0),
                               run_estim(1),
                               run_estim(2))
autoplot(vx)
