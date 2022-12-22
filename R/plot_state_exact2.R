#' function to extract local likelihoods along all branches
#' @param parameters used parameters for the likelihood caluclation
#' @param focal_tree used phylogeny
#' @param traits used traits
#' @param num_concealed_states number of concealed states
#' @param states the states matrix returned by cla_secsse_loglik
#' @param steps number of substeps evalualed per branch
#' @return nothing
#' @description this function will evaluate the log likelihood locally along
#' all branches and plot the result, using the phytools package.
#' @export
plot_state_exact2 <- function(parameters,
                             focal_tree,
                             traits,
                             num_concealed_states,
                             states,
                             steps) {
  
  message("collecting branch likelihoods\n")
  eval_res <- secsse::cla_secsse_eval(parameter = parameters,
                                      phy = focal_tree,
                                      traits = traits,
                                      num_concealed_states = num_concealed_states,
                                      ancestral_states = states,
                                      num_steps = steps,
                                      sampling_fraction = c(1, 1, 1),
                                      is_complete_tree = FALSE)
  
  message("converting collected likelihoods:\n")
  
  xs <- ape::node.depth.edgelength(focal_tree)
  ys <- ape::node.height(focal_tree)
  num_tips <- length(focal_tree$tip.label)
  num_nodes <- (1 + num_tips):length(ys)
  
  nodes <- data.frame(x = xs, y = ys, n = c(1:num_tips, num_nodes))
  
  to_plot <- eval_res
  to_plot[, c(1, 2)] <- to_plot[, c(1, 2)] + 1
  
  
  for_plot <- c()
  for (parent in unique(to_plot[, 1])) {
    for (daughter in unique(to_plot[, 2])) {
      indices <- which(to_plot[, 1] == parent & to_plot[, 2] == daughter)
      if (length(indices) > 0) {
        # we have a branch
        focal_branch <- to_plot[indices, ]
        start_x <- nodes$x[which(nodes$n == parent)]
        end_x <- nodes$x[which(nodes$n == daughter)]
        y <- nodes$y[which(nodes$n == daughter)]
        
        bl <- end_x - start_x
        
        probs <- focal_branch[, 4 + 6:11]
        probA <- (probs[, 1] + probs[,2] + probs[,3] ) / rowSums(probs)
        
        for (s in 1:(length(focal_branch[, 1]) - 1)) {
          x0 <- start_x + bl - focal_branch[s, 3]
          x1 <- start_x + bl - focal_branch[s + 1, 3]
          ps <- focal_branch[s, 4]
          for_plot <- rbind(for_plot, c(x0, x1, y, ps, parent, daughter))
        }
      }
    }  
  }
  
  node_bars <- c()
  for (parent in unique(to_plot[, 1])) {
    focal_data <- subset(to_plot, to_plot[, 1] == parent)
    daughters <- unique(focal_data[, 2])
    start_x <- nodes$x[which(nodes$n == parent)]
    y <- c()
    for (i in 1:length(daughters)) {
      y <- c(y, nodes$y[nodes$n == daughters[i]])
    }
    y <- sort(y)
    
    probs <- states[parent, 1 + 6:11]
    rel_prob <- sum(probs[1:3]) / sum(probs)
    
    node_bars <- rbind(node_bars, c(start_x, y, rel_prob))
  }
 colnames(for_plot) <- c("x0", "x1", "y", "prob", "p", "d")
  for_plot <- tibble::as_tibble(for_plot)
  colnames(node_bars) <- c("x", "y0", "y1", "prob")
  node_bars <- tibble::as_tibble(node_bars)
  
  px <- ggplot2::ggplot(for_plot, 
                        ggplot2::aes(x = x0, y = y, xend = x1, yend = y, col = prob)) +
    ggplot2::geom_segment() +
    ggplot2::geom_segment(data = node_bars, 
                          ggplot2::aes(x = x, y = y0, yend = y1, xend = x,
                                       col = prob)) +
    ggplot2::theme_classic() +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.line.y = ggplot2::element_blank())
  
  return(px)
}