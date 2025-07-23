test_that("plot idparlist", {
  idparslist <- list()
  focal_matrix <-
    secsse::create_default_lambda_transition_matrix(state_names = c("1", "2"),
                                                    model = "CR")
  idparslist[[1]] <- 
    secsse::create_lambda_list(state_names = c("1", "2"),
                               num_concealed_states = 2,
                               transition_matrix = focal_matrix,
                               model = "CR")
  idparslist[[2]] <- secsse::create_mu_vector(state_names = c("1", "2"),
                                              num_concealed_states = 2,
                                              model = "CR",
                                              lambda_list = idparslist[[1]])
  shift_mat <- secsse::create_default_shift_matrix(state_names = c("1", "2"),
                                                   num_concealed_states = 2,
                                                   mu_vector = idparslist[[2]])
  idparslist[[3]] <- secsse::create_q_matrix(state_names = c("1", "2"),
                                             num_concealed_states = 2,
                                             shift_matrix = shift_mat,
                                             diff.conceal = FALSE)
  focal_plot <- secsse::plot_idparslist(idparslist, 
                                        state_names = names(idparslist[[1]])) 
  
  testthat::expect_is(focal_plot$plot_qmat, "ggplot")
  testthat::expect_is(focal_plot$plot_lambda, "ggplot")

  # test if the plot can be rendered:
  testthat::expect_error(print(focal_plot$plot_qmat), NA)
  testthat::expect_error(print(focal_plot$plot_lambda), NA)
})