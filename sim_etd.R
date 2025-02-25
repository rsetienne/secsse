all_states <- c(0, 1)

focal_model <- "CTD"

idparslist <- list()
focal_matrix <-
  secsse::create_default_lambda_transition_matrix(state_names = all_states,
                                                  model = focal_model)
idparslist[[1]] <- 
  secsse::create_lambda_list(state_names = all_states,
                             num_concealed_states = 2,
                             transition_matrix = focal_matrix,
                             model = focal_model)
idparslist[[2]] <- secsse::create_mu_vector(state_names = all_states,
                                            num_concealed_states = 2,
                                            model = focal_model,
                                            lambda_list = idparslist[[1]])
shift_mat <- secsse::create_default_shift_matrix(state_names = all_states,
                                                 num_concealed_states = 2,
                                                 mu_vector = idparslist[[2]])
idparslist[[3]] <- secsse::create_q_matrix(state_names = all_states,
                                           num_concealed_states = 2,
                                           shift_matrix = shift_mat,
                                           diff.conceal = FALSE)

all_trees <- readRDS("/Users/thijsjanzen/Library/CloudStorage/Dropbox/projects/multi_trees/multi_trees3.rds")
all_traits <- readRDS("/Users/thijsjanzen/Library/CloudStorage/Dropbox/projects/multi_trees/multi_traits3.rds")
all_sf <- readRDS("/Users/thijsjanzen/Library/CloudStorage/Dropbox/projects/multi_trees/multi_sf3.rds")

all_sizes <- c()
for (i in seq_along(all_trees)) {
  all_sizes[i] <- treestats::number_of_lineages(all_trees[[i]])
}
quantile(all_sizes, c(0.025, 0.5, 0.975))

params <- readxl::read_xlsx("/Users/thijsjanzen/Library/CloudStorage/Dropbox/projects/multi_trees/res_rootedge.xlsx")
focal_params <- subset(params, params$Model == focal_model)

used_params <- as.numeric(as.vector(unlist(focal_params[1, 2:7])))

params <- list()
params[[1]] <- secsse::fill_in(idparslist[[1]], used_params)
params[[2]] <- secsse::fill_in(idparslist[[2]], used_params)
params[[3]] <- secsse::fill_in(idparslist[[3]], used_params)

secsse::cla_secsse_loglik(parameter = params,
                          phy = all_trees,
                          traits = all_traits,
                          num_concealed_states = 2,
                          sampling_fraction = all_sf,
                          take_into_account_root_edge = TRUE,
                          display_warning = FALSE)

# simulate data
test_trees <- list()
test_sizes <- c()
for (i in 1:length(all_trees)) {
  ca <- 0 
  if (length(all_trees[[i]]$tip.label) == 1) {
    ca <- all_trees[[i]]$root.edge + all_trees[[i]]$edge.length[[1]]
  } else {
    ca <- treestats::tree_height(all_trees[[i]])
  }
  
  
  
  #treestats::crown_age(all_trees[[i]])
  local_tree <- secsse::secsse_sim(lambdas = params[[1]],
                                   mus = params[[2]],
                                   qs = params[[3]],
                                   crown_age = ca,
                                   num_concealed_states = 2,
                                   min_spec = 0,
                                   start_at_crown = FALSE,
                                   conditioning = "none")
  test_trees[[i]] <- local_tree$phy
  test_sizes[i] <- length(local_tree$phy$tip.label) 
  cat(length(all_trees[[i]]$tip.label), ca, length(local_tree$phy$tip.label), "\n")
}

ref <- cbind(all_sizes, "empirical")
simm <- cbind(test_sizes, "simulated")

to_plot <- rbind(ref, simm)
colnames(to_plot) <- c("size", "type")

require(tidyverse)

to_plot <- as_tibble(to_plot)
to_plot$size <- as.numeric(to_plot$size)

#ggplot(to_plot, aes(x = type, y = size)) +
#  geom_boxplot()
ggplot(to_plot, aes(x = size, fill = type)) +
  geom_density(alpha = 0.7)  +
   scale_x_log10()


quantile(all_sizes, c(0.025, 0.5, 0.975))
quantile(test_sizes, c(0.025, 0.5, 0.975))











