
state_names <- c("H", "W")

lambda_trans <- secsse::create_default_lambda_transition_matrix(model = "CR",
                                                                state_names = state_names)
lambdas <- secsse::create_lambda_list(transition_matrix = lambda_trans,
                                      state_names = state_names,
                                      model = "CR")
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

set.seed(48)
#initpars <- c(0.3, 0.0, 0, 0)
initpars <- runif(4)
cat(initpars, "\n")
filled_in_pars <- list()
filled_in_pars[[1]] <- secsse::fill_in(idparslist[[1]], initpars)
filled_in_pars[[2]] <- secsse::fill_in(idparslist[[2]], initpars)
filled_in_pars[[3]] <- secsse::fill_in(idparslist[[3]], initpars)

l2 <- secsse::secsse_single_branch_loglik(parameter = filled_in_pars,
                                    phy = focal_tree,
                                    traits = "H",
                                    num_concealed_states = 2,
                                    sampling_fraction = c(1, 1),
                                    take_into_account_root_edge = FALSE,
                                    display_warning = FALSE,
                                    use_log_transform = TRUE)

l1 <- secsse::secsse_single_branch_loglik(parameter = filled_in_pars,
                                    phy = focal_tree,
                                    traits = "H",
                                    num_concealed_states = 2,
                                    sampling_fraction = c(1, 1),
                                    take_into_account_root_edge = FALSE,
                                    display_warning = FALSE)



all.equal(l1 ,l2)
l1
l2


focal_tree <- all_trees[[194]]
