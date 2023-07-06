#' Evaluation of probabilities of observing states along branches.
#' @title Likelihood for SecSSE model, using Rcpp
#' @param parameter list where the first is a table where lambdas across
#' different modes of speciation are shown, the second mus and the third
#'  transition rates.
#' @param phy phylogenetic tree of class phylo, ultrametric, fully-resolved,
#' rooted and with branch lengths.
#' @param traits vector with trait states, order of states must be the same as
#'  tree tips, for help, see vignette.
#' @param num_concealed_states number of concealed states, generally equivalent
#' to number of examined states.
#' @param ancestral_states ancestral states matrix provided by
#' cla_secsse_loglik, this is used as starting points for manual integration
#' @param num_steps number of steps to integrate along a branch
#' @param cond condition on the existence of a node root: 'maddison_cond',
#' 'proper_cond'(default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:'maddison_weigh
#' ,'proper_weights'(default) or 'equal_weights'. It can also be specified the
#' root state:the vector c(1,0,0) indicates state 1 was the root state.
#' @param sampling_fraction vector that states the sampling proportion per trait
#' state. It must have as many elements as trait states.
#' @param setting_calculation argument used internally to speed up calculation.
#' It should be leave blank (default : setting_calculation = NULL)
#' @param loglik_penalty the size of the penalty for all parameters; default is
#' 0 (no penalty)
#' @param is_complete_tree whether or not a tree with all its extinct species is
#' provided
#' @param method integration method used, available are:
#' "odeint::runge_kutta_cash_karp54", "odeint::runge_kutta_fehlberg78",
#' "odeint::runge_kutta_dopri5", "odeint::bulirsch_stoer" and
#' "odeint::runge_kutta4". Default method is:"odeint::bulirsch_stoer".
#' @param atol absolute tolerance of integration
#' @param rtol relative tolerance of integration
#' @param verbose provide intermediate verbose output if TRUE
#' @return The loglikelihood of the data given the parameters
#' @description Using see_ancestral_states = TRUE in the function
#' cla_secsse_loglik will provide posterior probabilities of the states of the
#' model on the nodes of the tree, but will not give the values on the branches.
#' This function evaluates these probabilities at fixed time intervals dt.
#' Because dt is fixed, this may lead to some inaccuracies, and dt is best
#' chosen as small as possible.
#' @export
cla_secsse_eval <- function(parameter,
                            phy,
                            traits,
                            num_concealed_states,
                            ancestral_states,
                            num_steps = NULL,
                            cond = "proper_cond",
                            root_state_weight = "proper_weights",
                            sampling_fraction,
                            setting_calculation = NULL,
                            loglik_penalty = 0,
                            is_complete_tree = FALSE,
                            method = "odeint::bulirsch_stoer",
                            atol = 1e-8,
                            rtol = 1e-7,
                            verbose = FALSE) {
  master_eval(parameter = parameter,
              phy = phy,
              traits = traits,
              num_concealed_states = num_concealed_states,
              ancestral_states = ancestral_states,
              cond = cond,
              root_state_weight = root_state_weight,
              sampling_fraction = sampling_fraction,
              setting_calculation = setting_calculation,
              loglik_penalty = loglik_penalty,
              is_complete_tree = is_complete_tree,
              atol = atol,
              rtol = rtol,
              method = method,
              num_steps = num_steps,
              verbose = verbose)
}

#' Logikelihood calculation for the SecSSE model given a set of parameters and
#' data, returning also the likelihoods along the branches
#' @title Likelihood for SecSSE model
#' @param parameter list where first vector represents lambdas, the second mus
#' and the third transition rates.
#' @param phy phylogenetic tree of class phylo, ultrametric, fully-resolved,
#' rooted and with branch lengths.
#' @param traits vector with trait states, order of states must be the same as
#' tree tips, for help, see vignette.
#' @param num_concealed_states number of concealed states, generally equivalent
#' to number of examined states.
#' @param ancestral_states ancestral states matrix provided by
#' secsse_loglik, this is used as starting points for the branch integration
#' @param cond condition on the existence of a node root: "maddison_cond",
#' "proper_cond"(default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:"maddison_weights",
#' "proper_weights"(default) or "equal_weights". It can also be specified the
#' root state:the vector c(1,0,0) indicates state 1 was the root state.
#' @param sampling_fraction vector that states the sampling proportion per
#' trait state. It must have as many elements as trait states.
#' @param setting_calculation argument used internally to speed up calculation.
#' It should be left blank (default : setting_calculation = NULL)
#' @param loglik_penalty the size of the penalty for all parameters; default is
#' 0 (no penalty)
#' @param is_complete_tree whether or not a tree with all its extinct species
#' is provided
#' @param atol absolute tolerance of integration
#' @param rtol relative tolerance of integration
#' @param method integration method used, available are:
#' "odeint::runge_kutta_cash_karp54", "odeint::runge_kutta_fehlberg78",
#' "odeint::runge_kutta_dopri5", "odeint::bulirsch_stoer" and
#' "odeint::runge_kutta4". Default method is:"odeint::bulirsch_stoer".
#' @param num_steps number of substeps to show intermediate likelihoods
#' along a branch, if left to NULL, the intermediate likelihoods at every
#' integration evaluation are stored, which is more exact, but can lead to
#' huge datasets / memory usage.
#' @param verbose provides intermediate output if TRUE
#' @return The loglikelihood of the data given the parameters
#' @examples
#' #' set.seed(5)
#' focal_tree <- ape::rphylo(n = 4, birth = 1, death = 0)
#' traits <- c(0, 1, 1, 0)
#' params <- secsse::id_paramPos(c(0, 1), 2)
#' params[[1]][] <- c(0.2, 0.2, 0.1, 0.1)
#' params[[2]][] <- 0.0
#' params[[3]][, ] <- 0.1
#' diag(params[[3]]) <- NA
#' #  Thus, we have for both, rates
#' # 0A, 1A, 0B and 1B. If we are interested in the posterior probability of
#' # trait 0 we have to provide a helper function that sums the probabilities of
#' # 0A and 0B, e.g.:
#' helper_function <- function(x) {
#'   return(sum(x[c(5, 7)]) / sum(x)) # normalized by total sum, just in case.
#' }
#' ll <- secsse::secsse_loglik(parameter = params,
#'                             phy = focal_tree,
#'                             traits = traits,
#'                             num_concealed_states = 2,
#'                             sampling_fraction = c(1, 1),
#'                             see_ancestral_states = TRUE)
#'
#' secsse_loglik_eval(parameter = params,
#'                    phy = focal_tree,
#'                    traits = traits,
#'                    ancestral_states = ll$states,
#'                    num_concealed_states = 2,
#'                    sampling_fraction = c(1, 1),
#'                    num_steps = 10)
#' @export
secsse_loglik_eval <- function(parameter,
                               phy,
                               traits,
                               num_concealed_states,
                               ancestral_states,
                               cond = "proper_cond",
                               root_state_weight = "proper_weights",
                               sampling_fraction,
                               setting_calculation = NULL,
                               loglik_penalty = 0,
                               is_complete_tree = FALSE,
                               atol = 1e-8,
                               rtol = 1e-7,
                               method = "odeint::bulirsch_stoer",
                               num_steps = NULL,
                               verbose = FALSE) {
  master_eval(parameter = parameter,
              phy = phy,
              traits = traits,
              num_concealed_states = num_concealed_states,
              ancestral_states = ancestral_states,
              cond = cond,
              root_state_weight = root_state_weight,
              sampling_fraction = sampling_fraction,
              setting_calculation = setting_calculation,
              loglik_penalty = loglik_penalty,
              is_complete_tree = is_complete_tree,
              atol = atol,
              rtol = rtol,
              method = method,
              num_steps = num_steps,
              verbose = verbose)
}

#' @keywords internal
master_eval <- function(parameter,
                        phy,
                        traits,
                        num_concealed_states,
                        ancestral_states,
                        cond = "proper_cond",
                        root_state_weight = "proper_weights",
                        sampling_fraction,
                        setting_calculation = NULL,
                        loglik_penalty = 0,
                        is_complete_tree = FALSE,
                        atol = 1e-8,
                        rtol = 1e-7,
                        method = "odeint::bulirsch_stoer",
                        num_steps = NULL,
                        verbose = FALSE) {
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  q_matrix <- parameter[[3]]

  check_input(traits,
              phy,
              sampling_fraction,
              root_state_weight,
              is_complete_tree)
  setting_calculation <- build_initStates_time(phy,
                                               traits,
                                               num_concealed_states,
                                               sampling_fraction,
                                               is_complete_tree,
                                               mus)

  for_time <- setting_calculation$forTime
  ances <- setting_calculation$ances

  if (is.list(lambdas)) {
    calcul <- c()
    ancescpp <- ances - 1
    forTimecpp <- for_time # nolint
    forTimecpp[, c(1, 2)] <- forTimecpp[, c(1, 2)] - 1 # nolint
    calcul <- cla_calThruNodes_store_cpp(ancescpp,
                                         ancestral_states,
                                         forTimecpp,
                                         lambdas,
                                         mus,
                                         q_matrix,
                                         method,
                                         atol,
                                         rtol,
                                         is_complete_tree,
                                         ifelse(is.null(num_steps),
                                                0,
                                                num_steps),
                                         verbose)
  } else {
    calcul <- calThruNodes_store_cpp(ances,
                                     ancestral_states,
                                     for_time,
                                     lambdas,
                                     mus,
                                     q_matrix,
                                     1,
                                     atol,
                                     rtol,
                                     method,
                                     is_complete_tree,
                                     ifelse(is.null(num_steps), 0, num_steps),
                                     verbose)
  }
  # if the number of steps == NULL, pass a 0.
  return(calcul)
}

#' function to plot the local probability along the tree, including the branches
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
#' @param steps number of substeps evaluated per branch, see description.
#' @param prob_func a function to calculate the probability of interest, see
#' description
#' @param is_complete_tree whether or not a tree with all its extinct species is
#' provided
#' @param verbose provides intermediate output (progressbars etc) when TRUE.
#' @return ggplot2 object
#' @description this function will evaluate the log likelihood locally along
#' all branches and plot the result. When steps is left to NULL, all likelihood
#' evaluations during integration are used for plotting. This may work for not
#' too large trees, but may become very memory heavy for larger trees. Instead,
#' the user can indicate a number of steps, which causes the probabilities to be
#' evaluated at a distinct amount of steps along each branch (and the
#' probabilities to be properly integrated in between these steps). This
#' provides an approximation, but generally results look very similar to using
#' the full evaluation.
#' The function used for prob_func will be highly dependent on your system.
#' for instance, for a 3 observed, 2 hidden states model, the probability
#' of state A is prob[1] + prob[2] + prob[3], normalized by the row sum.
#' prob_func will be applied to each row of the 'states' matrix (you can thus
#' test your function on the states matrix returned when
#' 'see_ancestral_states = TRUE'). Please note that the first N columns of the
#' states matrix are the extinction rates, and the (N+1):2N columns belong to
#' the speciation rates, where N = num_obs_states * num_concealed_states.
#'  A typical probfunc function will look like:
#' my_prob_func <- function(x) {
#'  return(sum(x[5:8]) / sum(x))
#' }
#' @examples
#' set.seed(5)
#' focal_tree <- ape::rphylo(n = 4, birth = 1, death = 0)
#' traits <- c(0, 1, 1, 0)
#' params <- secsse::id_paramPos(c(0, 1), 2)
#' params[[1]][] <- c(0.2, 0.2, 0.1, 0.1)
#' params[[2]][] <- 0.0
#' params[[3]][, ] <- 0.1
#' diag(params[[3]]) <- NA
#' #  Thus, we have for both, rates
#' # 0A, 1A, 0B and 1B. If we are interested in the posterior probability of
#' # trait 0,we have to provide a helper function that sums the probabilities of
#' # 0A and 0B, e.g.:
#' helper_function <- function(x) {
#'   return(sum(x[c(5, 7)]) / sum(x)) # normalized by total sum, just in case.
#' }
#'
#' out_plot <- plot_state_exact(parameters = params,
#'                              focal_tree = focal_tree,
#'                              traits = traits,
#'                              num_concealed_states = 2,
#'                              sampling_fraction = c(1, 1),
#'                              steps = 10,
#'                              prob_func = helper_function)
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
                             prob_func = NULL,
                             verbose = FALSE) {
  master_plot(parameters = parameters,
              focal_tree = focal_tree,
              traits = traits,
              num_concealed_states = num_concealed_states,
              sampling_fraction = sampling_fraction,
              cond = cond,
              root_state_weight = root_state_weight,
              is_complete_tree = is_complete_tree,
              method = method,
              atol = atol,
              rtol = rtol,
              steps = steps,
              prob_func = prob_func,
              verbose = verbose)
}

#' function to plot the local probability along the tree,
#' including the branches, for the CLA model.
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
#' @param steps number of substeps evaluated per branch, see description.
#' @param prob_func a function to calculate the probability of interest, see
#' description
#' @param is_complete_tree whether or not a tree with all its extinct species is
#' provided
#' @param verbose return verbose output / progress bars when true.
#' @return ggplot2 object
#' @description this function will evaluate the log likelihood locally along
#' all branches and plot the result. When steps is left to NULL, all likelihood
#' evaluations during integration are used for plotting. This may work for not
#' too large trees, but may become very memory heavy for larger trees. Instead,
#' the user can indicate a number of steps, which causes the probabilities to be
#' evaluated at a distinct amount of steps along each branch (and the
#' probabilities to be properly integrated in between these steps). This
#' provides an approximation, but generally results look very similar to using
#' the full evaluation.
#' The function used for prob_func will be highly dependent on your system.
#' for instance, for a 3 observed, 2 hidden states model, the probability
#' of state A is prob[1] + prob[2] + prob[3], normalized by the row sum.
#' prob_func will be applied to each row of the 'states' matrix (you can thus
#' test your function on the states matrix returned when
#' 'see_ancestral_states = TRUE'). Please note that the first N columns of the
#' states matrix are the extinction rates, and the (N+1):2N columns belong to
#' the speciation rates, where N = num_obs_states * num_concealed_states.
#' A typical probfunc function will look like:
#' my_prob_func <- function(x) {
#'  return(sum(x[5:8]) / sum(x))
#' }
#'
#' @examples
#' set.seed(13)
#'phylotree <- ape::rcoal(12, tip.label = 1:12)
#'traits <- sample(c(0, 1, 2), ape::Ntip(phylotree), replace = TRUE)
#'num_concealed_states <- 3
#'sampling_fraction <- c(1,1,1)
#'phy <- phylotree
#'# the idparlist for a ETD model (dual state inheritance model of evolution)
#'# would be set like this:
#'idparlist <- secsse::cla_id_paramPos(traits,num_concealed_states)
#'lambd_and_modeSpe <- idparlist$lambdas
#'lambd_and_modeSpe[1,] <- c(1,1,1,2,2,2,3,3,3)
#'idparlist[[1]] <- lambd_and_modeSpe
#'idparlist[[2]][] <- 0
#'masterBlock <- matrix(4,ncol = 3, nrow = 3, byrow = TRUE)
#'diag(masterBlock) <- NA
#'idparlist[[3]] <- q_doubletrans(traits, masterBlock, diff.conceal = FALSE)
#'# Now, internally, clasecsse sorts the lambda matrices, so they look like
#'#  a list with 9 matrices, corresponding to the 9 states
#'# (0A,1A,2A,0B, etc)

#'parameter <- idparlist
#'lambda_and_modeSpe <- parameter$lambdas
#'lambda_and_modeSpe[1,] <- c(0.2,0.2,0.2,0.4,0.4,0.4,0.01,0.01,0.01)
#'parameter[[1]] <- prepare_full_lambdas(traits,num_concealed_states,
#'                                       lambda_and_modeSpe)
#'parameter[[2]] <- rep(0,9)
#'masterBlock <- matrix(0.07, ncol = 3, nrow = 3, byrow = TRUE)
#'diag(masterBlock) <- NA
#'parameter[[3]] <- q_doubletrans(traits, masterBlock, diff.conceal = FALSE)
#'helper_function <- function(x) {
#'  return(sum(x[c(10, 13, 16)]) / sum(x))
#'}
#'out_plot <- plot_state_exact_cla(parameters = parameter,
#'                                 focal_tree = phy,
#'                                 traits = traits,
#'                                 num_concealed_states = 3,
#'                                 sampling_fraction = sampling_fraction,
#'                                 cond = 'maddison_cond',
#'                                 root_state_weight = 'maddison_weights',
#'                                 is_complete_tree = FALSE,
#'                                 prob_func = helper_function,
#'                                 steps = 10)
#' @export
plot_state_exact_cla <- function(parameters,
                                 focal_tree,
                                 traits,
                                 num_concealed_states,
                                 sampling_fraction,
                                 cond = "proper_cond",
                                 root_state_weight = "proper_weights",
                                 is_complete_tree = FALSE,
                                 method = "odeint::bulirsch_stoer",
                                 atol = 1e-8,
                                 rtol = 1e-7,
                                 steps = 10,
                                 prob_func = NULL,
                                 verbose = FALSE) {

  master_plot(parameters = parameters,
                     focal_tree = focal_tree,
                     traits = traits,
                     num_concealed_states = num_concealed_states,
                     sampling_fraction = sampling_fraction,
                     cond = cond,
                     root_state_weight = root_state_weight,
                     is_complete_tree = is_complete_tree,
                     method = method,
                     atol = atol,
                     rtol = rtol,
                     steps = steps,
                     prob_func = prob_func,
                     verbose = verbose)
}

#' @keywords internal
master_plot <- function(parameters,
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
                        steps = 10,
                        prob_func = NULL,
                        verbose = FALSE) {

  if (is.null(prob_func)) {
    stop("need to set a probability function, check description to how")
  }

  if (verbose) message("collecting all states on nodes")
  ll1 <-  master_loglik(parameter = parameters,
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

  if (verbose) message("collecting branch likelihoods\n")
  eval_res <- master_eval(parameter = parameters,
                          phy = focal_tree,
                          traits = traits,
                          num_concealed_states =
                            num_concealed_states,
                          ancestral_states = ll1$states,
                          cond = cond,
                          root_state_weight = root_state_weight,
                          num_steps = steps,
                          sampling_fraction = sampling_fraction,
                          is_complete_tree = is_complete_tree,
                          atol = atol,
                          rtol = rtol,
                          method = method,
                          verbose = verbose)

  if (verbose) message("\nconverting collected likelihoods
                       to graph positions:\n")

  xs <- ape::node.depth.edgelength(focal_tree)
  ys <- ape::node.height(focal_tree)
  num_tips <- length(focal_tree$tip.label)
  num_nodes <- (1 + num_tips):length(ys)

  nodes <- data.frame(x = xs, y = ys, n = c(1:num_tips, num_nodes))

  to_plot <- eval_res
  if (is.list(parameters[[1]]))  to_plot[, c(1, 2)] <- to_plot[, c(1, 2)] + 1

  for_plot <- collect_branches(to_plot, nodes, prob_func, verbose)

  node_bars <- collect_node_bars(to_plot, nodes, prob_func, ll1)

  if (verbose) message("\ngenerating ggplot object\n")

  focal_plot <- make_ggplot(for_plot, node_bars)
  return(focal_plot)
}


#' @importFrom rlang .data
#' @keywords internal
make_ggplot <- function(for_plot, node_bars) {
  ggplot_plot <- ggplot2::ggplot(for_plot) +
    ggplot2::geom_segment(ggplot2::aes(x = .data[["x0"]],
                                       y = .data[["y"]],
                                       xend = .data[["x1"]],
                                       yend = .data[["y"]],
                                       col = .data[["prob"]])) +
    ggplot2::geom_segment(data = node_bars,
                          ggplot2::aes(x = .data[["x"]],
                                       y = .data[["y0"]],
                                       yend = .data[["y1"]],
                                       xend = .data[["x"]],
                                       col = .data[["prob"]])
    ) +
    ggplot2::theme_classic() +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank())

  return(ggplot_plot)
}

#' @keywords internal
collect_branches <- function(to_plot,
                             nodes,
                             prob_func,
                             verbose) {
  num_rows <- length(to_plot[, 1])

  for_plot <- matrix(nrow = num_rows, ncol = 6)
  for_plot_cnt <- 1
  if (verbose) pb <- utils::txtProgressBar(max = length(unique(to_plot[, 1])),
                                           style = 3)
  cnt <- 1
  for (parent in unique(to_plot[, 1])) {
    if (verbose) utils::setTxtProgressBar(pb, cnt)
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

        probs <- apply(focal_branch[, 4:length(focal_branch[1, ])],
                       1,
                       prob_func)
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
  colnames(for_plot) <- c("x0", "x1", "y", "prob", "p", "d")
  for_plot <- tibble::as_tibble(for_plot)
  return(for_plot)
}

#' @keywords internal
collect_node_bars <- function(to_plot,
                              nodes,
                              prob_func,
                              ll) {
  node_bars <- matrix(nrow = length(unique(to_plot[, 1])), ncol = 4)
  node_bars_cnt <- 1
  for (parent in unique(to_plot[, 1])) {
    focal_data <- subset(to_plot, to_plot[, 1] == parent)
    daughters <- unique(focal_data[, 2])
    start_x <- nodes$x[which(nodes$n == parent)]
    y <- c()
    for (i in seq_along(daughters)) {
      y <- c(y, nodes$y[nodes$n == daughters[i]])
    }
    y <- sort(y)

    probs <- ll$states[parent, ]
    rel_prob <- prob_func(probs)
    node_bars[node_bars_cnt, ] <- c(start_x, y, rel_prob)
    node_bars_cnt <- node_bars_cnt + 1
  }

  colnames(node_bars) <- c("x", "y0", "y1", "prob")
  node_bars <- tibble::as_tibble(node_bars)
  return(node_bars)
}
