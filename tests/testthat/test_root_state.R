test_that("check root state", {
  
  num_concealed_states <- 2
  
  focal_matrix <-
    secsse::create_default_lambda_transition_matrix(state_names = c("S", "N"),
                                                    model = "CR")
  lambda_list_CR <- secsse::create_lambda_list(state_names = c("S", "N"),
                                               num_concealed_states =
                                                 num_concealed_states,
                                               transition_matrix = focal_matrix,
                                               model = "CR")
  mus_CR <- secsse::create_mu_vector(state_names = c("S", "N"),
                                     num_concealed_states =
                                       num_concealed_states,
                                     model = "CR",
                                     lambda_list = lambda_list_CR)
  t_CR <- secsse::create_default_shift_matrix(state_names = c("S", "N"),
                                              num_concealed_states =
                                                num_concealed_states,
                                              mu_vector = mus_CR)
  q_CR <- secsse::create_q_matrix(state_names = c("S", "N"),
                                  num_concealed_states = num_concealed_states,
                                  shift_matrix = t_CR,
                                  diff.conceal = FALSE)
  
  idparslist <- list()
  idparslist[[1]] <- lambda_list_CR
  idparslist[[2]] <- mus_CR
  idparslist[[3]] <- q_CR
  
  pars <- c(0.3, 0.0, 0.2, 0.2)
  params <- list()
  for(i in 1:3) {
    params[[i]] <- secsse::fill_in(idparslist[[i]], pars)
  }
  
  phy <- secsse::secsse_sim(lambdas = params[[1]],
                            mus = params[[2]],
                            qs = params[[3]],
                            crown_age = 5,
                            num_concealed_states = num_concealed_states,
                            sampling_fraction = c(1, 1))
  
  
  
  ances_res <- secsse::cla_secsse_loglik(parameter = params,
                                         phy = phy$phy,
                                         traits = phy$obs_traits,
                                         num_concealed_states =
                                           num_concealed_states,
                                         sampling_fraction = c(1, 1),
                                         see_ancestral_states = TRUE,
                                         return_root_state = FALSE,
                                         display_warning = FALSE)
  
  root_res <- secsse::cla_secsse_loglik(parameter = params,
                                         phy = phy$phy,
                                         traits = phy$obs_traits,
                                         num_concealed_states =
                                          num_concealed_states,
                                         sampling_fraction = c(1, 1),
                                         see_ancestral_states = FALSE,
                                         return_root_state = TRUE,
                                         display_warning = FALSE)
  
  

   root_s <- ances_res$ancestral_states[1, ]
   rs <- as.vector(root_res$root_state)
   testthat::expect_true(all.equal(root_s, rs))
})