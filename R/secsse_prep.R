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

#' @keywords internal
get_state_names <- function(state_names, num_concealed_states) {
  num_obs_states <- length(state_names)

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
  return(all_state_names)
}

#' Helper function to automatically create lambda matrices, based on input.
#' When choosing the CTD model, rates associated with observed states are now
#' re-distributed to concealed states. This implicitly assumes that the number
#' of observed and concealed states is identical.
#' 
#' @inheritParams default_params_doc
#' 
#' @examples
#' trans_matrix <- c(0, 0, 0, 1)
#' trans_matrix <- rbind(trans_matrix, c(1, 1, 1, 2))
#' lambda_list <- create_lambda_list(state_names = c(0, 1),
#'                                   num_concealed_states = 2,
#'                                   transition_matrix = trans_matrix,
#'                                   model = "ETD")
#'
#' @export
create_lambda_list <- function(state_names = c(0, 1),
                               num_concealed_states = 2,
                               transition_matrix,
                               model = "ETD") {
  if (!(model %in% c("CR", "ETD", "CTD"))) {
    stop("only CR, ETD or CTD are specified")
  }

  num_obs_states <- length(state_names)
  total_num_states <- num_obs_states * num_concealed_states

  all_state_names <- get_state_names(state_names, num_concealed_states)

  lambdas <- list()
  for (i in 1:total_num_states) {
    lambdas[[i]] <- matrix(0, nrow = total_num_states,
                           ncol = total_num_states)
    rownames(lambdas[[i]]) <- all_state_names
    colnames(lambdas[[i]]) <- all_state_names
  }
  names(lambdas) <- all_state_names

  transition_list <- convert_transition_list(transition_matrix, state_names)

  if (model == "ETD" || model == "CR") {
  
    # ETD settings
    for (i in seq_len(nrow(transition_list))) {
      focal_state <- transition_list[i, 1]
      daughter1   <- transition_list[i, 2]
      daughter2   <- transition_list[i, 3]
      target_rate <- transition_list[i, 4]
  
      for (j in seq_len(num_concealed_states)) {
        incr <- (j - 1) * num_obs_states
        focal_rate <- target_rate
        #if (model == "CTD") focal_rate <- concealed_spec_rates[j]
  
        lambdas[[focal_state + incr]][daughter1 + incr,
                                      daughter2 + incr] <- focal_rate
        lambdas[[focal_state + incr]][daughter2 + incr,
                                      daughter1 + incr] <- focal_rate
      }
    }
  }
  
  if (model == "CTD") {
    for (i in seq_len(nrow(transition_list))) {
      focal_state <- (transition_list[i, 1] - 1) * num_obs_states
      daughter1   <- (transition_list[i, 2] - 1) * num_obs_states
      daughter2   <- (transition_list[i, 3] - 1) * num_obs_states
      target_rate <- transition_list[i, 4]
      
      for (incr in 1:num_concealed_states) {
        
        lambdas[[focal_state + incr]][daughter1 + incr,
                                      daughter2 + incr] <- target_rate
        lambdas[[focal_state + incr]][daughter2 + incr,
                                      daughter1 + incr] <- target_rate
      }
    }
  }
  return(lambdas)
}

#' Helper function to neatly setup a Q matrix, without transitions to
#' concealed states (only observed transitions shown)
#' 
#' @inheritParams default_params_doc
#' 
#' @return transition matrix
#' @examples
#' shift_matrix <- c(0, 1, 5)
#' shift_matrix <- rbind(shift_matrix, c(1, 0, 6))
#' q_matrix <- secsse::create_q_matrix(state_names = c(0, 1),
#'                                     num_concealed_states = 2,
#'                                     shift_matrix = shift_matrix,
#'                                     diff.conceal = TRUE)
#' @export
create_q_matrix <- function(state_names,
                            num_concealed_states,
                            shift_matrix,
                            diff.conceal = FALSE) {

  total_num_states <- length(state_names)
  if (total_num_states != num_concealed_states) {
    stop("number of concealed states has to be equal to the number of observed states")
  }
  trans_matrix <- matrix(0, ncol = total_num_states,
                         nrow = total_num_states)

  transition_list <- convert_transition_list_q(shift_matrix, state_names)
  for (i in seq_len(nrow(transition_list))) {
    parent_state <- transition_list[i, 1]
    daughter_state <- transition_list[i, 2]
    focal_rate <- transition_list[i, 3]
    trans_matrix[parent_state, daughter_state] <- focal_rate
  }

  diag(trans_matrix) <- NA

  trans_matrix <- secsse::q_doubletrans(traits = state_names,
                                        masterBlock = trans_matrix,
                                        diff.conceal = diff.conceal)
  
  all_state_names <- get_state_names(state_names, num_concealed_states)
  colnames(trans_matrix) <- all_state_names
  rownames(trans_matrix) <- all_state_names
  diag(trans_matrix) <- NA

  return(trans_matrix)
}

#' Function to expand an existing q_matrix to a number of concealed states
#' 
#' @inheritParams default_params_doc
#' 
#' @note This is highly similar to [q_doubletrans()].
#' 
#' @return updated q matrix
#' @export
expand_q_matrix <- function(q_matrix,
                            num_concealed_states,
                            diff.conceal = FALSE) {
  
  warning("This function is deprecated, please use q_doubletrans, piping result
          to q_doubletrans, this may introduce inaccuracies.")
  
  traits <- get_state_names(names(q_matrix), num_concealed_states)
  
  return(secsse::q_doubletrans(traits = traits,
                               masterBlock = q_matrix,
                               diff.conceal = diff.conceal))
}

#' Helper function to create a default `shift_matrix` list
#' 
#' This function generates a generic shift matrix to be used with the function
#' [create_q_matrix()].
#' 
#' @inheritParams default_params_doc
#' 
#' @examples
#' shift_matrix <- create_default_shift_matrix(state_names = c(0, 1),
#'                                             num_concealed_states = 2,
#'                                             mu_vector = c(1, 2, 1, 2))
#' q_matrix <- create_q_matrix(state_names = c(0, 1),
#'                             num_concealed_states = 2,
#'                             shift_matrix = shift_matrix,
#'                             diff.conceal = FALSE)
#' @export
create_default_shift_matrix <- function(state_names = c("0", "1"),
                                        num_concealed_states = 2,
                                        mu_vector = NULL) {
  lm <- unlist(mu_vector)
  focal_rate <- max(lm) + 1
  num_obs_states <- length(state_names)
  transition_list <- c()
  for (i in 1:num_obs_states) {
    for (j in 1:num_obs_states) {
      if (i != j) {
        start_state <- state_names[i]
        end_state <- state_names[j]
        to_add <- c(start_state, end_state, focal_rate)
        transition_list <- rbind(transition_list, to_add)
        focal_rate <- focal_rate + 1
      }
    }
  }

  rownames(transition_list) <- rep("", nrow(transition_list))
  return(transition_list)
}

#' Helper function to create a default lambda list
#' 
#' This function generates a generic lambda list, assuming no transitions
#' between states, e.g. a species of observed state 0 generates daughter
#' species with state 0 as well.
#'
#' @inheritParams default_params_doc
#' 
#' @examples
#' lambda_matrix <-
#'      create_default_lambda_transition_matrix(state_names = c(0, 1),
#'                                              model = "ETD")
#' lambda_list <- create_lambda_list(state_names = c(0, 1),
#'                                   num_concealed_states = 2,
#'                                   transition_matrix = lambda_matrix,
#'                                   model = "ETD")
#' @export
create_default_lambda_transition_matrix <- function(state_names = c("0", "1"),
                                                    model = "ETD") {
  transition_list <- c()
  for (i in seq_along(state_names)) {
    focal_rate <- i
    if (model == "CR") focal_rate <- 1
    transition_list <- rbind(transition_list,
                             c(state_names[i],
                               state_names[i],
                               state_names[i],
                               focal_rate))
  }
  rownames(transition_list) <- rep("", nrow(transition_list))
  return(transition_list)
}

#' Generate mus vector
#'
#' @inheritParams default_params_doc
#'
#' @return mu vector
#' @export
create_mu_vector <- function(state_names,
                             num_concealed_states,
                             model = "CR",
                             lambda_list) {
  focal_rate <- 1 + max(unlist(lambda_list), na.rm = TRUE)

  if (!(model %in% c("CR", "ETD", "CTD"))) {
    stop("only CR, ETD or CTD are specified")
  }

  all_names <- get_state_names(state_names, num_concealed_states)

  mus <- rep(focal_rate, length(all_names))

  num_obs_states <- length(state_names)

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

  names(mus) <- all_names
  return(mus)
}

#' @keywords internal
replace_matrix <- function(focal_matrix,
                           params,
                           is_lambda = FALSE) {
  if (is_lambda) {
    entries <- table(focal_matrix)
    entries <- entries[names(entries) != 0]
    for (i in seq_len(nrow(focal_matrix))) {
      for (j in seq_len(ncol(focal_matrix))) {
        if (focal_matrix[i, j] != 0 && !is.na(focal_matrix[i, j]))  {
          index <- focal_matrix[i, j]
          
          new_val <- params[index]
          if (index %in% names(entries)) {
            mult <- 1 / entries[names(entries) == index]
            new_val <- new_val * mult
          }
          
          focal_matrix[i, j] <- new_val
        }
      }
    }
  } else {
    for (i in seq_len(nrow(focal_matrix))) {
      for (j in seq_len(ncol(focal_matrix))) {
        if (focal_matrix[i, j] != 0 && !is.na(focal_matrix[i, j]))  {
          index <- focal_matrix[i, j]
          focal_matrix[i, j] <- params[index]
        }
      }
    }
  }
  return(focal_matrix)
}

#' Helper function to enter parameter value on their right place
#' 
#' @inheritParams default_params_doc
#' @return lambda matrices, `q_matrix` or mu vector with the correct values in
#'   their right place.
#' @export
fill_in <- function(object,
                    params) {
  if (is.list(object)) { # lambda matrix
    for (k in seq_along(object)) {
      object[[k]] <- replace_matrix(object[[k]], params, is_lambda = FALSE)
    }
  } else if (is.matrix(object)) {
    object <- replace_matrix(object, params)
  } else if (is.vector(object)) {
    for (k in seq_along(object)) {
      if (object[k] != 0 && !is.na(object[k])) {
        object[k] <- params[object[k]]
      }
    }
  }
  return(object)
}

#' @keywords internal
extract_answ <- function(indic_mat,
                         param_mat,
                         answ) {
  for (i in seq_len(nrow(indic_mat))) {
    for (j in seq_len(ncol(param_mat))) {
      if (indic_mat[i, j] > 0 && !is.na(indic_mat[i, j])) {
        index <- indic_mat[i, j]
        answ[index] <- param_mat[i, j]
      }
    }
  }
  return(answ)
}


#' Extract parameter values out of the result of a maximum likelihood inference 
#' run
#' 
#' @inheritParams default_params_doc

#' @return Vector of parameter estimates.
#' @export
extract_par_vals <- function(param_posit,
                             ml_pars) {
  if (length(param_posit) != length(ml_pars)) {
    stop("param posit doesn't match ml_pars in structure")
  }

  answ <- c()
  for (i in seq_along(param_posit[[1]])) {
    answ <- extract_answ(param_posit[[1]][[i]],
                         ml_pars[[1]][[i]],
                         answ)
  }

  answ <- extract_answ(param_posit[[3]], # Q matrix
                       ml_pars[[3]],
                       answ)
  for (i in seq_along(param_posit[[2]])) {
    if (param_posit[[2]][i] > 0 && !is.na(param_posit[[2]][i])) {
      answ[param_posit[[2]][i]] <- ml_pars[[2]][i]
    }
  }
  return(answ)
}
