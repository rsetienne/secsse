#' function to extract local likelihoods along all branches
#' @param parameters used parameters for the likelihood caluclation
#' @param focal_tree used phylogeny
#' @param traits used traits
#' @param num_concealed_states number of concealed states
#' @param states the states matrix returned by cla_secsse_loglik
#' @param dt chosen time step
#' @param res resolution of the color scale used.
#' @return nothing
#' @description this function will evaluate the log likelihood locally along
#' all branches and plot the result, using the phytools package.
#' @export
plot_state_exact <- function(parameters,
                             focal_tree,
                             traits,
                             num_concealed_states,
                             states,
                             dt,
                             res) {
  
  min_branch_len <- min(focal_tree$edge.length)
  if (dt > min_branch_len) {
    warning("dt chosen too large, changed to be 1/100th of the shortest branch length")
    dt <- min_branch_len / 100
  }
  
  eval_res <- secsse::cla_secsse_eval(parameter = parameters,
                                      phy = focal_tree,
                                      traits = traits,
                                      num_concealed_states = num_concealed_states,
                                      ancestral_states = states,
                                      dt = dt,
                                      sampling_fraction = c(1, 1, 1),
                                      is_complete_tree = TRUE)
  
  res <- 50
  
  local_maps <- list()
  for (i in 1:length(focal_tree$edge[, 1])) {
    anc <- focal_tree$edge[i, 1]
    desc <- focal_tree$edge[i, 2]
    
    # probs are binned in 1/1000 bins
    # we need stretches of same bins
    focal_data <- subset(eval_res, 
                         eval_res[, 1] == anc - 1 & 
                         eval_res[, 2] == desc - 1)
    to_add <- diff(focal_data[, 3])
    used_bins <- floor(focal_data[1:length(to_add), 4] * res)
    names(to_add) <- used_bins
    local_maps[[i]] <- to_add
  }
  
  tree_to_plot <- focal_tree
  tree_to_plot$maps <- local_maps
  class(tree_to_plot) <- c("simmap", setdiff(class(focal_tree), "simmap"))
  used_cols <- viridis::plasma(n = res, begin = 0, end = 1) 
  names(used_cols) <- 1:res
  plot(tree_to_plot, colors = used_cols, lwd = 3)
  N <- length(focal_tree$tip.label)
  phytools::add.color.bar(leg = 0.3, 
                          title = "prob A\n", 
                          cols = used_cols, 
                          x = 0.1, 
                          y = 1 + 0.08*(N - 1),
                          prompt = FALSE)
}