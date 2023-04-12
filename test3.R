transform_parameters <- function(param_posit, param_id,
                                 param_values) {
  
  new_param_posit <- param_posit
  
  for (j in 1:nrow(new_param_posit[[3]])) {
    new_param_posit[[1]][[j]][, ] <- 0
  }
  
  for (j in 2:3) {
    new_param_posit[[j]][] <- 0
  }
  
  for (i in seq_along(param_id)) {
    
    for (j in 1:nrow(param_posit[[3]])) {
      id <- which(param_posit[[1]][[j]] == param_id[i])
      new_param_posit[[1]][[j]][id] <- param_values[i]
    }
    for (j in 2:3) {
      id <- which(param_posit[[j]] == param_id[i])
      new_param_posit[[j]][id] <- param_values[i]
    }
  }
  return(new_param_posit)
}

b <- 1.0
d <- 0.01
qr <- 0.0
params <- c(b, d, qr,
            b, d, qr) # transition rate

param_numbers <- 1:length(params) 

traits <- c(0, 1)

param_posit <- secsse::cla_id_paramPos(traits, num_concealed_states = 1)

namez <- names(param_posit$mus)

param_posit$mus[1] <- d; param_posit$mus[2] <- d

full_lambdas <- list()
for (i in 1:2) {
  full_lambdas[[i]] <- matrix(0, 2, 2)
}

full_lambdas[[1]][1, 1] <- b
full_lambdas[[2]][2, 2] <- b

param_posit$lambdas <- full_lambdas

param_posit$Q[1, 2] <- qr
param_posit$Q[2, 1] <- qr

diag(param_posit$Q) <- 0

total_tried <- 0

for (r in 1:1000) {

tree <- secsse::secsse_sim(lambdas = param_posit$lambdas,
                           mus = param_posit$mus,
                           qs = param_posit$Q,
                           crown_age = 5,
                           pool_init_states = c("0A", "1A"),
                           maxSpec = 999,
                           conditioning = "none",
                           max_tries = 1e7,
                           verbose = FALSE)
}
