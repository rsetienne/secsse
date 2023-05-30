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
  if (!(model %in% c("CR", "ETD", "CTD"))) {
    stop("only CR, ETD or CTD are specified")
  }
  
  num_obs_states <- length(state_names)
  total_num_states <- num_obs_states * num_concealed_states
  
  concealed_state_names <- LETTERS[1:num_concealed_states]
  all_state_names <- c()
  cnt <- 1
  for (j in 1:num_concealed_states) {
    for (i in 1:num_obs_states) {
    
      all_state_names[cnt] <- paste0(state_names[i], 
                                     concealed_state_names[j])
      cnt <- cnt + 1
    }
  }
  
  lambdas <- list()
  for (i in 1:total_num_states) {
    lambdas[[i]] <- matrix(0, nrow = total_num_states,
                              ncol = total_num_states)
    rownames(lambdas[[i]]) <- all_state_names
    colnames(lambdas[[i]]) <- all_state_names
    
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
      if (model == "CR") focal_rate <- 1
      
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


#' helper function to create a default q_matrix list
#' @param state_names names of the observed states
#' @description
#' This function generates a generic transition list
#' @export
create_default_q_list <- function(state_names = c("0", "1"),
                                  num_concealed_states,
                                  lambda_matrices = NULL) {
  
  lm <- unlist(lambda_matrices)
  focal_rate <- max(lm) + 1
  
  num_obs_states <- length(state_names)
  concealed_state_names <- LETTERS[1:num_concealed_states]
  
  transition_list <- c()
  # now, we need to find those entries with rate A->B, or rate B->A
  for (j in 1:num_concealed_states) {
      for (k in 1:num_concealed_states) {    
        if (j != k) {
           # transition of concealed state j -> k
          for (i in 1:num_obs_states) {
            start_state <- paste0(state_names[i], concealed_state_names[j])
            end_state   <- paste0(state_names[i], concealed_state_names[k])
            to_add <- c(start_state, end_state, focal_rate)
            transition_list <- rbind(transition_list, to_add)
          }
          focal_rate <- focal_rate + 1
        }
      }
  }
  
  return(transition_list)
}


#' helper function to create a default transition list
#' @param state_names names of the observed states
#' @param consider_combinations should there be extra 
#' states considering combinations of the observed states?
#' @description
#' This function generates a generic transition list, assuming no transitions
#' between states, e.g. a species of observed state 0 generates two daughter
#' species with state 0 as well. 
#' If consider_combinations is set to TRUE, the function automatically
#' generates transitions for the combined states as well, for example if there
#' are two states A and B, it adds a third state AB, which generates two unique
#' daughter species, in this case each having either trait A, or B.
#' @export
create_default_transition_list <- function(state_names = c("0", "1"),
                                           consider_combinations = FALSE) {
  transition_list <- c()
  
  for (i in 1:length(state_names)) {
    transition_list <- rbind(transition_list, c(state_names[i], state_names[i], state_names[i], i))
  }
  cnt <- length(state_names)
  if (consider_combinations) {
    for (i in 1:length(state_names)) {
      for (j in i:length(state_names)) {
        focal_state_name <- paste0(state_names[i], state_names[j])
        to_add <- c(focal_state_name, state_names[i], state_names[j], cnt + 1)
        cnt <- cnt + 1
        transition_list <- rbind(transition_list, to_add)
      }
    }
  }
  
  return(transition_list)
}

#' function to generate generic mus vector
#' @param state_names full state names, including concealed states, for example
#' c("0A", "1A", "0B", "1B")
#' @param num_concealed_states number of concealed states
#' @param model model replicated, available are "CR", "ETD" and "CTD"
#' @param q_matrix previously generated q_matrix, used to infer the rate
#' number to start with
#' @return mu vector
#' @export
create_mus <- function(state_names, 
                       num_concealed_states,
                       model = "CR",
                       q_matrix) {
  focal_rate <- 1 + max(q_matrix, na.rm = TRUE)
  if (!(model %in% c("CR", "ETD", "CTD"))) {
    stop("only CR, ETD or CTD are specified")
  }
  
  mus <- rep(focal_rate, length(state_names))
  
  num_obs_states <- length(state_names) / num_concealed_states
  if (model == "ETD") {
    for (i in 1:num_obs_states) {
      indices <- seq(i, length(mus), by = num_concealed_states)
      mus[indices] <- focal_rate
      focal_rate <- focal_rate + 1
    }     
  }
  if (model == "CTD") {
    mus <- c()
    for (i in 1:num_obs_states) {
      mus <- c(mus, rep(focal_rate, num_concealed_states))
      focal_rate <- focal_rate + 1
    }
  }
  
  names(mus) <- state_names
  return(mus)
}

#' @keywords internal
replace_matrix <- function(focal_matrix,
                           params) {
  for (i in 1:nrow(focal_matrix)) {
    for (j in 1:ncol(focal_matrix)) {
      if (focal_matrix[i, j] != 0 && !is.na(focal_matrix[i, j]))  {
        index <- focal_matrix[i, j]
        focal_matrix[i, j] <- params[index]
      }
    }
  }
  
  return(focal_matrix)
}



#' helper function to enter parameter value on their right place
#' @param object lambda matrices, q_matrix or mu vector
#' @param params parameters in order, where each value reflects the value
#' of the parameter at that position, e.g. c(0.3, 0.2, 0.1) will fill out
#' the value 0.3 for the parameter with rate indentifier 1, 0.2 for the
#' parameter with rate identifier 2 and 0.1 for the parameter with rate
#' identifier 3
#' @export
fill_in <- function(object,
                    params) {
  
  if (is.list(object)) { # lambda matrix
    for (k in 1:length(object)) {
      object[[k]] <- replace_matrix(object[[k]], params)
    }
  } else if (is.matrix(object)) {
    object <- replace_matrix(object, params)
  } else if (is.vector(object)) {
    for (k in 1:length(object)) {
      if (object[k] != 0 && !is.na(object[k])) {
        object[k] <- params[object[k]]
      }
    }
  }
  return(object)
}
