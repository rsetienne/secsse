#' @keywords internal
construct_row <- function(v, state_names) {
  out <- c(0, 0, 0, as.numeric(v[4]))
  out[1] <- which(state_names == v[1])
  out[2] <- which(state_names == v[2])
  out[3] <- which(state_names == v[3])
  return(out)
}

#' @keywords internal
convert_transition_list <- function(transition_list, state_names) {
  res <- apply(transition_list, 1, construct_row, state_names)
  return(t(res))
}

#' @keywords internal
construct_row_q <- function(v, state_names) {
  out <- c(0, 0, as.numeric(v[3]))
  out[1] <- which(state_names == v[1])
  out[2] <- which(state_names == v[2])
  return(out)
}

#' @keywords internal
convert_transition_list_q <- function(transition_list, state_names) {
  res <- apply(transition_list, 1, construct_row_q, state_names)
  return(t(res))
}



#' helper function to automatically create lambda matrices, based on input
#' @param state_names vector of names of all observed states
#' @param num_concealed_states number of hidden states
#' @param transition_list a matrix containing a description of all speciation
#' events, where the first column indicates the source state, the second and
#' third column indicate the two daughter states, and the fourth column gives 
#' the rate indicator used. E.g.: ["SA", "S", "A", 1] for a trait state "SA"
#' which upon speciation generates two daughter species with traits "S" and "A",
#' where the number 1 is used as indicator for optimization of the likelihood.
#' @param model used model, choice of "ETD" (Examined Traits Diversification) or
#' "CTD" (Concealed Traits Diversification).
#' @param concealed_spec_rates vector specifying the rate indicators for each
#' concealed state, length should be identical to num_concealed_states. If left
#' empty when using the CTD model, it is assumed that all available speciation
#' rates are distributed uniformly over the concealed states.
#' @export
create_lambda_matrices <- function(state_names,
                                   num_concealed_states,
                                   transition_list,
                                   model = "ETD",
                                   concealed_spec_rates = NULL) {
  
  if (!(model %in% c("ETD", "CTD"))) {
    stop("only ETD or CTD are specified")
  }
  
  num_obs_states <- length(state_names)
  total_num_states <- num_obs_states * num_concealed_states
  
  lambdas <- list()
  for (i in 1:total_num_states) {
    lambdas[[i]] <- matrix(0, nrow = total_num_states,
                              ncol = total_num_states)
  }
  
  transition_list <- convert_transition_list(transition_list, state_names)

  
  if (model == "CTD") {
    if (is.null(concealed_spec_rates)) {
      spec_rates <- unique(transition_list[, 4])
      spec_rates <- sort(spec_rates)
      num_times <- ceiling(num_concealed_states / length(spec_rates))
      concealed_spec_rates <- rep(spec_rates, num_times)
      concealed_spec_rates <- concealed_spec_rates[1:num_concealed_states]
      concealed_spec_rates <- sort(concealed_spec_rates)
    }
  }
  
  # ETD settings
  for (i in 1:nrow(transition_list)) {
    focal_state <- transition_list[i, 1]
    daughter1   <- transition_list[i, 2]
    daughter2   <- transition_list[i, 3]
    target_rate <- transition_list[i, 4]
    
    for (j in seq_len(num_concealed_states)) {
      incr <- (j - 1) * num_obs_states
      focal_rate <- target_rate
      if (model == "CTD") focal_rate <- concealed_spec_rates[j]
      lambdas[[focal_state + incr]][daughter1 + incr,
                                    daughter2 + incr] <- focal_rate
      lambdas[[focal_state + incr]][daughter2 + incr, 
                                    daughter1 + incr] <- focal_rate
    }
  }
  
  return(lambdas)
}

#' helper function to neatly setup a transition matrix
#' @param state_names names of observed states
#' @param transition_list matrix of transitions, indicating in order: 1) 
#' starting state (typically the column in the transition matrix), 2) ending 
#' state (typically the row in the transition matrix) and 3) associated rate 
#' indicator
#' @return transition matrix
#' @export
create_transition_matrix <- function(state_names,
                                     transition_list) {
  num_obs_states <- length(state_names)
 
  trans_matrix <- matrix(0, ncol = num_obs_states,
                            nrow = num_obs_states)
  
  transition_list <- convert_transition_list_q(transition_list, state_names)
  
  for (i in 1:nrow(transition_list)) {
    parent_state <- transition_list[i, 1]
    daughter_state <- transition_list[i, 2]
    focal_rate <- transition_list[i, 3]
    trans_matrix[daughter_state, parent_state] <- focal_rate
  }
  
  colnames(trans_matrix) <- state_names
  rownames(trans_matrix) <- state_names
  diag(trans_matrix) <- NA
  return(trans_matrix)
}
