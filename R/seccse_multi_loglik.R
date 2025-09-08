#' @keywords internal
multi_loglik <- function(parameter,
                         phy,
                         traits,
                         num_concealed_states,
                         cond = "proper_cond",
                         root_state_weight = "proper_weights",
                         sampling_fraction,
                         setting_calculation = NULL,
                         see_ancestral_states = FALSE,
                         loglik_penalty = 0,
                         is_complete_tree = FALSE,
                         take_into_account_root_edge = FALSE,
                         num_threads = 1,
                         atol = 1e-8,
                         rtol = 1e-7,
                         method = "odeint::runge_kutta_cash_karp54",
                         display_warning = FALSE,
                         use_normalization = TRUE,
                         return_root_state = FALSE) {
  
  res <- list()
  root_states <- list()
  for (i in 1:length(phy)) {
    if (is.list(sampling_fraction)) {
      focal_sampling_fraction <- sampling_fraction[[i]]
    } else {
      focal_sampling_fraction <- sampling_fraction
    }
    
    if (is.list(root_state_weight)) {
      if (sum(is.na(root_state_weight[[i]])) || 
          length(root_state_weight[[i]]) < 1) {
        focal_root_state_weight <- "proper_weights"
      } else {
        focal_root_state_weight <- root_state_weight[[i]]
      }
    } else {
      focal_root_state_weight <- root_state_weight
    }
    
    focal_setting_calculation <- NULL
    if (is.list(setting_calculation)) {
      focal_setting_calculation <- setting_calculation[[i]]
    }
    
    if (length(phy[[i]]$tip.label) == 1) {
      res[[i]] <- secsse::secsse_single_branch_loglik(parameter = parameter,
                                                      phy = phy[[i]],
                                                      traits = traits[[i]],
                                                      num_concealed_states =
                                                        num_concealed_states,
                                                      cond = cond,
                                                      root_state_weight = 
                                                        focal_root_state_weight,
                                                      sampling_fraction = 
                                                        focal_sampling_fraction,
                                                      setting_calculation = 
                                                        focal_setting_calculation,
                                                      see_ancestral_states = FALSE,
                                                      loglik_penalty = loglik_penalty,
                                                      is_complete_tree = 
                                                        is_complete_tree,
                                                      take_into_account_root_edge = 
                                                        take_into_account_root_edge,
                                                      num_threads = num_threads,
                                                      atol = atol,
                                                      rtol = rtol,
                                                      method = method,
                                                      display_warning = display_warning,
                                                      use_normalization = use_normalization)$loglik
    } else {
      local_res <- secsse_loglik(parameter = parameter,
                                phy = phy[[i]],
                                traits = traits[[i]],
                                num_concealed_states = num_concealed_states,
                                cond = cond,
                                root_state_weight = focal_root_state_weight,
                                sampling_fraction = focal_sampling_fraction,
                                setting_calculation = focal_setting_calculation,
                                see_ancestral_states = FALSE,
                                loglik_penalty = loglik_penalty,
                                is_complete_tree = is_complete_tree,
                                take_into_account_root_edge = 
                                  take_into_account_root_edge,
                                num_threads = num_threads,
                                atol = atol,
                                rtol = rtol,
                                method = method,
                                display_warning = display_warning,
                                use_normalization = use_normalization,
                                return_root_state = return_root_state) 
      
      if (return_root_state) {
        root_states[[i]] <- local_res$root_state
      } else {
        res[[i]] <- local_res
      }
    }
  }
  
  ll <- do.call(sum, res)
  
  if (return_root_state) {
    return(list(LL = ll,
                root_state = root_states))
  }
  
  return(ll) 
}
