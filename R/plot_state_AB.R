#' function to plot the hidden-states across a tree
#' @param focal_tree the used phylogeny
#' @param ancestral_states ancestral states matrix, as returned by 
#' cla_secsse_loglik
#' @return ggplot object
#' @description This plots the probability of observing hidden state A along the
#' tree. Please note that this takes the probabilities at the nodes and 
#' interpolates the probabilities along the branches. Furthermore, it assumes
#' the tips have 50/50 probability of either state. 
#' @export
plot_state_AB <- function(focal_tree, 
                          ancestral_states) {

  if (length(ancestral_states[1, ]) != 6) {
    stop("this only works for two hidden states")
  }
  prob_A <- ancestral_states[, 1] + ancestral_states[, 2] + ancestral_states[, 3]
  prob_B <- ancestral_states[, 4] + ancestral_states[, 5] + ancestral_states[, 6]
  if (sum((prob_A + prob_B) >= 1-1e-9) != length(prob_A)) {
    stop("probabilities don't add up to 1")
  }
  
  require(ggtree)
  
  c <- data.frame(node = rownames(ancestral_states), signal = prob_A)
  a <- data.frame(node = 1:length(focal_tree$tip.label),
                  signal = rep(0.5, length(focal_tree$tip.label)))
  d.1 <- rbind(a, c)
  d.1$node <- as.numeric(d.1$node)
  d.1$signal <- as.numeric(d.1$signal)
  
  colour.vector <- rep("blue", length(d.1$node))
  d.2 <- cbind(d.1, colour.vector)
  
  tree.2 <- dplyr::full_join(focal_tree, d.2, by = 'node')
  
  t1 <- ggtree::ggtree(tree.2, 
                       aes(color = signal),  
                       ladderize = TRUE, continuous = "color") +
    ggplot2::scale_color_gradientn(colours = c('red', 'blue'), 
                                   limits = c(0, 1)) +
    ggtree::geom_tiplab(size = 1) + 
    ggtree::theme(legend.position = c(.05, .85))
  
  return(t1) 
}

#' function to plot the observed states across a tree
#' @param focal_tree the used phylogeny
#' @param ancestral_states ancestral states matrix, as returned by 
#' cla_secsse_loglik
#' @return list with as many plots as traits
#' @description This plots the probability of observing each tip state across
#' the tree. Coloring is interpolated between known probabilities at the nodes.
#' Shown are separate plots for each trait.
#' @export
plot_state_traits <- function(focal_tree, traits, ancestral_states) {
  if (length(ancestral_states[1, ]) != 6) {
    stop("this only works for two hidden states")
  }
  require(ggtree)
  
  plots <- list()
  max_col <- c("red", "blue", "purple")
  
  namez <- c("Mainland", "Island", "Cosmopolitan")
  
  cnt <- 1
  for (tt in unique(traits)) {
    focal_prob <- ancestral_states[, 1 + tt] + ancestral_states[, 4 + tt]
    
    c <- data.frame(node = rownames(ancestral_states), signal = focal_prob)
    a <- data.frame(node = 1:length(focal_tree$tip.label),
                    signal = as.numeric(traits == tt))
    d.1 <- rbind(a, c)
    d.1$node <- as.numeric(d.1$node)
    d.1$signal <- as.numeric(d.1$signal)
    
    colour.vector <- rep("blue", length(d.1$node))
    d.2 <- cbind(d.1, colour.vector)
    
    tree.2 <- dplyr::full_join(focal_tree, d.2, by = 'node')
    
    t1 <- ggtree::ggtree(tree.2, 
                         aes(color = signal),  
                         ladderize = TRUE, continuous = "color") +
      ggplot2::scale_color_gradientn(colours = c('grey', max_col[cnt]), 
                                     limits = c(0, 1)) +
      ggtree::geom_tiplab(size = 1) + 
      ggtree::theme(legend.position = c(.05, .85)) +
      ggplot2::ggtitle(namez[cnt])
    plots[[cnt]] <- t1
    cnt <- cnt + 1
  }
  return(plots) 
}