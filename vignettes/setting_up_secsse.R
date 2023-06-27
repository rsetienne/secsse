## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----default_trans_list-------------------------------------------------------
focal_list <- secsse::create_default_transition_list(state_names = c("S", "N"))
focal_list

## ----default lambda matrices--------------------------------------------------
lambda_matrices <- secsse::create_lambda_matrices(state_names = c("S", "N"),
                                                  num_concealed_states = 2, 
                                                  transition_list = focal_list,
                                                  model = "CR")
lambda_matrices

## ----adding extinction--------------------------------------------------------
mus <- secsse::create_mus(state_names = colnames(lambda_matrices[[1]]),
                          num_concealed_states = 2,
                          model = "CR",
                          lambdas = lambda_matrices)
mus

## ----default_trans------------------------------------------------------------
q_list <- secsse::create_default_q_list(state_names = c("S", "N"),
                                        num_concealed_states = 2,
                                        mus = mus)

q_list

trans_matrix <- secsse::create_transition_matrix(state_names = colnames(lambda_matrices[[1]]),
                                                 transition_list = q_list)
trans_matrix

## ----fill in parameters-------------------------------------------------------

speciation_S <- 0.5
speciation_N <- 0.3
extinction <- 0.0
q_ab <- 0.5
q_ba <- 0.5
sp_ns <- 0.2
sp_sn <- 0.2

params <- c(speciation_S, 
            speciation_N,
            extinction, 
            sp_ns, sp_sn, q_ab, q_ba)

lambda_matrices_p <- secsse::fill_in(lambda_matrices,
                                     params)
trans_matrix_p <- secsse::fill_in(trans_matrix,
                                  params)
mus_p <- secsse::fill_in(mus,
                         params)

## ----simulate tree------------------------------------------------------------
simulated_tree <- secsse::secsse_sim(lambdas = lambda_matrices_p,
                                     mus = mus_p,
                                     qs = trans_matrix_p,
                                     num_concealed_states = 2,
                                     crown_age = 5,
                                     conditioning = "obs_states",
                                     verbose = TRUE)
sim_traits <- simulated_tree$obs_traits
focal_tree <- simulated_tree$phy

## ----maximum likelihood-------------------------------------------------------
param_posit <- list()
param_posit[[1]] <- lambda_matrices
param_posit[[2]] <- mus
param_posit[[3]] <- trans_matrix

initpars <- params
initpars[initpars == 0.0] <- 1e-5

answ <- secsse::cla_secsse_ml(phy = focal_tree,
                              traits = sim_traits,
                              num_concealed_states = 2,
                              idparslist = param_posit,
                              idparsopt = c(1, 2, 3, 4, 5, 6, 7),
                              initparsopt = initpars,
                              idparsfix = c(0),
                              parsfix = c(0.0),
                              sampling_fraction = c(1, 1),
                              optimmethod = "subplex",
                              verbose = FALSE,
                              num_threads = 6)

## ----extract_pars-------------------------------------------------------------
found_pars_vals <- secsse::extract_par_vals(param_posit, answ$MLpars)
found_pars_vals

## ----define_model_function----------------------------------------------------
fit_model <- function(tree, traits, model) {
  focal_list <- secsse::create_default_transition_list(state_names = c("S", "N"))
  lambda_matrices <- secsse::create_lambda_matrices(state_names = c("S", "N"),
                                                    num_concealed_states = 2, 
                                                    transition_list = focal_list,
                                                    model = model)
  mus <- secsse::create_mus(state_names = colnames(trans_matrix),
                            num_concealed_states = 2,
                            model = model,
                            lambdas = lambda_matrices)
  q_list <- secsse::create_default_q_list(state_names = c("S", "N"),
                                          num_concealed_states = 2,
                                          mus = mus)
  
  trans_matrix <- secsse::create_transition_matrix(state_names = 
                                                     colnames(lambda_matrices[[1]]),
                                                   transition_list = q_list)
  
  param_posit <- list()
  param_posit[[1]] <- lambda_matrices
  param_posit[[2]] <- mus
  param_posit[[3]] <- trans_matrix
  
  max_indicator <- max(trans_matrix, na.rm = TRUE)
  
  # we cheat a bit by setting extinction to zero -
  # in a real analysis this should be avoided.
  extinct_rates <- unique(mus)
  idparsopt <- 1:max_indicator
  idparsopt <- idparsopt[-extinct_rates]
  idparsfix <- c(0, extinct_rates)
  parsfix <- rep(0.0, length(idparsfix))
  
  initpars <- runif(n = length(idparsopt))

  testthat::expect_output( # suppress output
  answ <- secsse::cla_secsse_ml(phy = focal_tree,
                                traits = traits,
                                num_concealed_states = 2,
                                idparslist = param_posit,
                                idparsopt = idparsopt,
                                initparsopt = initpars,
                                idparsfix = idparsfix,
                                parsfix = parsfix,
                                sampling_fraction = c(1, 1),
                                optimmethod = "subplex",
                                verbose = FALSE,
                                num_threads = 6)
  )
  found_pars_vals <- secsse::extract_par_vals(param_posit, answ$MLpars)
  aic <- 2 * max_indicator - 2 * as.numeric(answ$ML)
  return(list(pars = found_pars_vals,
              ml = as.numeric(answ$ML),
              aic = aic))
}

## ----model looping------------------------------------------------------------
found <- c()
for (focal_model in c("CR", "CTD", "ETD")) {
  local_answ <- fit_model(tree = focal_tree, 
                          traits = sim_traits,
                          model = focal_model)
  to_add <- c(focal_model, local_answ$ml, local_answ$aic)
  found <- rbind(found, to_add)
}
colnames(found) <- c("model", "LL", "AIC")
found <- as.data.frame(found)
found$LL <- as.numeric(found$LL)
found$AIC <- as.numeric(found$AIC)
found 

