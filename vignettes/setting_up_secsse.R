## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----default_trans_list-------------------------------------------------------
used_states <- c("S", "N")
lambda_transition_matrix <- secsse::create_default_lambda_transition_matrix(state_names = used_states,
                                                 model = "CR")
lambda_transition_matrix

## ----default lambda matrices--------------------------------------------------
num_hidden_states <- 2
lambda_list <- secsse::create_lambda_list(state_names = used_states,
                                                  num_concealed_states = num_hidden_states,
                                                  transition_matrix = lambda_transition_matrix,
                                                  model = "CR")
lambda_list

## ----adding extinction--------------------------------------------------------
mus <- secsse::create_mu_vector(state_names = used_states,
                                num_concealed_states = num_hidden_states,
                                model = "CR",
                                lambda_list = lambda_list)
mus

## ----default_trans------------------------------------------------------------
shift_matrix <- secsse::create_default_shift_matrix(state_names = used_states,
                                        num_concealed_states = num_hidden_states,
                                        mus = mus)

shift_matrix

q_matrix <- secsse::create_q_matrix(state_names = used_states,
                                    num_concealed_states = num_hidden_states,
                                    shift_matrix = shift_matrix,
                                    diff.conceal = TRUE)
q_matrix

## ----fill in parameters-------------------------------------------------------

speciation <- 0.5
extinction <- 0.0
sp_sn <- 0.2
sp_ns <- 0.2
q_ab <- 0.5
q_ba <- 0.5

params <- c(speciation,
            extinction,
            sp_sn, sp_ns,
            q_ab, q_ba)

# we use the suffix p to signal that these are filled in with the params
lambda_list_p <- secsse::fill_in(lambda_list,
                                 params)
q_matrix_p <- secsse::fill_in(q_matrix,
                              params)
mus_p <- secsse::fill_in(mus,
                         params)

## ----simulate tree------------------------------------------------------------
simulated_tree <- secsse::secsse_sim(lambdas = lambda_list_p,
                                     mus = mus_p,
                                     qs = q_matrix_p,
                                     num_concealed_states = num_hidden_states,
                                     crown_age = 5,
                                     conditioning = "obs_states",
                                     verbose = TRUE,
                                     seed = 26)
sim_traits <- simulated_tree$obs_traits
focal_tree <- simulated_tree$phy

## ----maximum likelihood-------------------------------------------------------
param_posit <- list()
param_posit[[1]] <- lambda_list
param_posit[[2]] <- mus
param_posit[[3]] <- q_matrix

initpars <- params
initpars <- initpars[-2]

answ <- secsse::cla_secsse_ml(phy = focal_tree,
                              traits = sim_traits,
                              num_concealed_states = num_hidden_states,
                              idparslist = param_posit,
                              idparsopt = c(1, 3, 4, 5, 6),
                              initparsopt = initpars,
                              idparsfix = c(0, 2),
                              parsfix = c(0.0, 0.0),
                              sampling_fraction = c(1, 1),
                              optimmethod = "subplex",
                              verbose = FALSE,
                              num_threads = 6,
                              atol = 0.1, # high values for demonstration 
                              rtol = 0.1) # purposes, don't use at home!

## ----extract_pars-------------------------------------------------------------
found_pars_vals <- secsse::extract_par_vals(param_posit, answ$MLpars)
found_pars_vals

## ----define_model_function----------------------------------------------------
fit_model <- function(focal_tree, traits, model) {
  focal_list <- secsse::create_default_lambda_transition_matrix(state_names = used_states,
                                                   model = model)
  lambda_matrices <- secsse::create_lambda_list(state_names = used_states,
                                                    num_concealed_states = num_hidden_states,
                                                    transition_matrix =
                                                        focal_list,
                                                    model = model)
  mus <- secsse::create_mu_vector(state_names = used_states,
                                  num_concealed_states = num_hidden_states,
                                  model = model,
                                  lambda_list = lambda_matrices)
  q_list <- secsse::create_default_shift_matrix(state_names = used_states,
                                          num_concealed_states = num_hidden_states,
                                          mus = mus)

  trans_matrix <- secsse::create_q_matrix(state_names = used_states,
                                          num_concealed_states = num_hidden_states,
                                          shift_matrix = q_list,
                                          diff.conceal = TRUE)

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

  initpars <- c(rep(params[1], min(extinct_rates) - 1),
                params[-c(1, 2)])

  answ <- secsse::cla_secsse_ml(phy = focal_tree,
                                traits = traits,
                                num_concealed_states = num_hidden_states,
                                idparslist = param_posit,
                                idparsopt = idparsopt,
                                initparsopt = initpars,
                                idparsfix = idparsfix,
                                parsfix = parsfix,
                                sampling_fraction = c(1, 1),
                                optimmethod = "subplex",
                                verbose = FALSE,
                                num_threads = 6,
                                atol = 0.1, # high values for demonstration 
                                rtol = 0.1) # purposes, don't use at home!
  found_pars_vals <- secsse::extract_par_vals(param_posit, answ$MLpars)
  aic <- 2 * max_indicator - 2 * as.numeric(answ$ML)
  return(list(pars = found_pars_vals,
              ml = as.numeric(answ$ML),
              aic = aic))
}

## ----model looping------------------------------------------------------------

found <- c()
for (focal_model in c("CR", "CTD", "ETD")) {
  local_answ <- fit_model(focal_tree = focal_tree,
                          traits = sim_traits,
                          model = focal_model)
  found <- rbind(found, c(focal_model, local_answ$ml, local_answ$aic))
}
colnames(found) <- c("model", "LL", "AIC")
found <- as.data.frame(found)
found$LL <- as.numeric(found$LL)
found$AIC <- as.numeric(found$AIC)
found

