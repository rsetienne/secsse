#' function to extract local likelihoods along all branches
#' @param parameters used parameters for the likelihood calculation
#' @param focal_tree used phylogeny
#' @param traits used traits
#' @param num_concealed_states number of concealed states
#' @param sampling_fraction sampling fraction
#' @param cond condition on the existence of a node root: 'maddison_cond',
#' 'proper_cond'(default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:'maddison_weigh
#' ,'proper_weights'(default) or 'equal_weights'. It can also be specified the
#' root state:the vector c(1,0,0) indicates state 1 was the root state.
#' @param method integration method used, available are: 
#' "odeint::runge_kutta_cash_karp54", "odeint::runge_kutta_fehlberg78", 
#' "odeint::runge_kutta_dopri5", "odeint::bulirsch_stoer" and 
#' "odeint::runge_kutta4". Default method is:"odeint::bulirsch_stoer".
#' @param atol absolute tolerance of integration
#' @param rtol relative tolerance of integration
#' @param steps number of substeps evaluated per branch
#' @param prob_func a function to calculate the probability of interest, see 
#' description
#' @param is_complete_tree whether or not a tree with all its extinct species is
#' provided
#' @return ggplot2 object
#' @description this function will evaluate the log likelihood locally along
#' all branches and plot the result.
#' the function used for prob_func will be highly dependent on your system.
#' for instance, for a 3 observed, 2 hidden states model, the probability
#' of state A is prob[1] + prob[2] + prob[3], normalized by the row sum.
#' probfunc will be applied to each row of the 'states' matrix (you can thus
#' test your function on the states matrix returned when 
#' 'see_ancestral_states = TRUE'). A typical probfunc function will look like:
#' my_prob_func <- function(x) {
#'  return(sum(x[5:8]) / sum(x))
#' }
#' @export
plot_state_exact <- function(parameters,
                             focal_tree,
                             traits,
                             num_concealed_states,
                             sampling_fraction,
                             cond = "proper_cond",
                             root_state_weight = "proper_weights",
                             is_complete_tree = FALSE,
                             method = "odeint::bulirsch_stoer",
                             atol = 1e-16,
                             rtol = 1e-16,
                             steps = NULL,
                             prob_func = NULL) {
  
  if (is.null(prob_func)) {
    stop("need to set a probability function, check description to how")
  }
  
  message("collecting all states on nodes")
  ll1 <- secsse::secsse_loglik(parameter = parameters,
                               phy = focal_tree,
                               traits = traits,
                               num_concealed_states = num_concealed_states,
                               cond = cond,
                               root_state_weight = root_state_weight,
                               sampling_fraction = sampling_fraction,
                               see_ancestral_states = TRUE,
                               loglik_penalty = 0,
                               is_complete_tree = is_complete_tree,
                               num_threads = 1,
                               atol = atol,
                               rtol = rtol,
                               method = method)
  
  message("collecting branch likelihoods\n")
  eval_res <- secsse::secsse_loglik_eval(parameter = parameters,
                                         phy = focal_tree,
                                         traits = traits,
                                         num_concealed_states = num_concealed_states,
                                         ancestral_states = ll1$states,
                                         cond = cond,
                                         root_state_weight = root_state_weight,
                                         sampling_fraction = sampling_fraction,
                                         is_complete_tree = is_complete_tree,
                                         atol = atol,
                                         rtol = rtol,
                                         method = method,
                                         num_steps = steps)
  
  message("\nconverting collected likelihoods to graph positions:\n")
  
  xs <- ape::node.depth.edgelength(focal_tree)
  ys <- ape::node.height(focal_tree)
  num_tips <- length(focal_tree$tip.label)
  num_nodes <- (1 + num_tips):length(ys)
  
  nodes <- data.frame(x = xs, y = ys, n = c(1:num_tips, num_nodes))
  
  to_plot <- eval_res
  # to_plot[, c(1, 2)] <- to_plot[, c(1, 2)] + 1
  
  
  num_rows <- length(to_plot[, 1])   #length(unique(to_plot[, 1])) * 2 * steps
  
  for_plot <- matrix(nrow = num_rows, ncol = 6)
  for_plot_cnt <- 1
  pb <- txtProgressBar(max = length(unique(to_plot[, 1])), style = 3)
  cnt <- 1
  for (parent in unique(to_plot[, 1])) {
    setTxtProgressBar(pb, cnt)
    cnt <- cnt + 1
    
    to_plot2 <- subset(to_plot, to_plot[, 1] == parent)
    for (daughter in unique(to_plot2[, 2])) {
      indices <- which(to_plot2[, 2] == daughter)
      if (length(indices) > 0) {
        # we have a branch
        focal_branch <- to_plot2[indices, ]
        start_x <- nodes$x[which(nodes$n == parent)]
        end_x <- nodes$x[which(nodes$n == daughter)]
        y <- nodes$y[which(nodes$n == daughter)]
        
        bl <- end_x - start_x
        
        probs <- apply(focal_branch[, 4:length(focal_branch[1, ])], 1, prob_func)
        
        for (s in 1:(length(focal_branch[, 1]) - 1)) {
          x0 <- start_x + bl - focal_branch[s, 3]
          x1 <- start_x + bl - focal_branch[s + 1, 3]
          ps <- probs[s]
          for_plot[for_plot_cnt, ] <- c(x0, x1, y, ps, parent, daughter)
          for_plot_cnt <- for_plot_cnt + 1
        }
      }
    }  
  }
  
  node_bars <- matrix(nrow = length(unique(to_plot[, 1])), ncol = 4)
  node_bars_cnt <- 1
  for (parent in unique(to_plot[, 1])) {
    focal_data <- subset(to_plot, to_plot[, 1] == parent)
    daughters <- unique(focal_data[, 2])
    start_x <- nodes$x[which(nodes$n == parent)]
    y <- c()
    for (i in 1:length(daughters)) {
      y <- c(y, nodes$y[nodes$n == daughters[i]])
    }
    y <- sort(y)
    
    num_states <- length(ll1$states[1, ]) / 2
    max_row_size <- length(ll1$states[1, ])
    
    probs <- ll1$states[parent, ]
    rel_prob <- prob_func(probs)
    # node_bars <- rbind(node_bars, c(start_x, y, rel_prob))
    node_bars[node_bars_cnt, ] <- c(start_x, y, rel_prob)
    node_bars_cnt <- node_bars_cnt + 1
  }
  
  colnames(for_plot) <- c("x0", "x1", "y", "prob", "p", "d")
  for_plot <- tibble::as_tibble(for_plot)
  colnames(node_bars) <- c("x", "y0", "y1", "prob")
  node_bars <- tibble::as_tibble(node_bars)
  
  message("\ngenerating ggplot object\n")
  
  px <- ggplot2::ggplot(for_plot, 
                        ggplot2::aes(x = x0, y = y, xend = x1, yend = y, col = prob)) +
    ggplot2::geom_segment(linewidth = 2) +
    ggplot2::geom_segment(data = node_bars, 
                          ggplot2::aes(x = x, y = y0, yend = y1, xend = x,
                                       col = prob),
                          linewidth = 1.5
    ) +
    ggplot2::theme_classic() +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank())
  
  return(px)
}