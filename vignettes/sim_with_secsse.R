## ----setup_params-------------------------------------------------------------
spec_matrix <- c(0, 0, 0, 1)
spec_matrix <- rbind(spec_matrix, c(1, 1, 1, 1))
lambda_list <- secsse::create_lambda_list(state_names = c(0, 1),
                                          num_concealed_states = 2,
                                          transition_matrix = spec_matrix,
                                          model = "CR")

mu_vector <- secsse::create_mu_vector(state_names = c(0, 1),
                                   num_concealed_states = 2,
                                   model = "CR",
                                   lambda_list = lambda_list)

shift_matrix <- c(0, 1, 3)
shift_matrix <- rbind(shift_matrix, c(1, 0, 4))

q_matrix <- secsse::create_q_matrix(state_names = c(0, 1),
                                    num_concealed_states = 2,
                                    shift_matrix = shift_matrix,
                                    diff.conceal = FALSE)

## ----enter parameters---------------------------------------------------------
speciation_rate <- 0.5
extinction_rate <- 0.05
q_ab <- 0.1
q_ba <- 0.1
used_params <- c(speciation_rate, extinction_rate, q_ab, q_ba)

sim_lambda_list <- secsse::fill_in(lambda_list, used_params)
sim_mu_vector   <- secsse::fill_in(mu_vector, used_params)
sim_q_matrix    <- secsse::fill_in(q_matrix, used_params)

## ----simulate_tree------------------------------------------------------------
sim_tree <- secsse::secsse_sim(lambdas = sim_lambda_list,
                               mus = sim_mu_vector,
                               qs = sim_q_matrix,
                               crown_age = 5,
                               num_concealed_states = 2,
                               seed = 5)

if (requireNamespace("diversitree")) {
  traits_for_plot <- data.frame(trait = as.numeric(sim_tree$obs_traits),
                                row.names = sim_tree$phy$tip.label)
  diversitree::trait.plot(tree = sim_tree$phy,
                          dat = traits_for_plot,
                          cols = list("trait" = c("blue", "red")),
                          type = "p")
} else {
  plot(sim_tree$phy)
}


## ----conditioning-------------------------------------------------------------
sim_tree <- secsse::secsse_sim(lambdas = sim_lambda_list,
                               mus = sim_mu_vector,
                               qs = sim_q_matrix,
                               crown_age = 5,
                               num_concealed_states = 2,
                               conditioning = "obs_states",
                               seed = 6)
sim_tree$obs_traits
sim_tree$true_traits

sim_tree <- secsse::secsse_sim(lambdas = sim_lambda_list,
                               mus = sim_mu_vector,
                               qs = sim_q_matrix,
                               crown_age = 5,
                               num_concealed_states = 2,
                               conditioning = "true_states",
                               seed = 6)
sim_tree$obs_traits
sim_tree$true_traits

