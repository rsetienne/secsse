test_that("the loglik for the complete tree under cla_secsse", {
  Sys.unsetenv("R_TESTS")
  utils::data("example_phy_GeoSSE", package = "secsse")
  traits <- as.numeric(example_phy_GeoSSE$tip.state)
  lambdas <- list()
  lambdas[[1]] <- matrix(0, ncol = 9, nrow = 9, byrow = TRUE)
  lambdas[[1]][2, 1] <- 1.5
  lambdas[[1]][3, 1] <- 0.5
  lambdas[[1]][3, 2] <- 1
  for (i in 2:9) {
    lambdas[[i]] <- lambdas[[1]]
  }
  mus <- rep(0, 9)
  Q <- matrix(stats::runif(81), ncol = 9, nrow = 9, byrow = TRUE)
  parameter <- list()
  parameter[[1]] <- lambdas
  parameter[[2]] <- mus
  parameter[[3]] <- Q

  num_concealed_states <- 3
  sampling_fraction <- c(1, 1, 1)

  secsse_cla_LL3 <- secsse_loglik(parameter = parameter,
                                  phy = example_phy_GeoSSE,
                                  traits = traits,
                                  num_concealed_states =
                                    num_concealed_states,
                                  cond = "maddison_cond",
                                  root_state_weight = "maddison_weights",
                                  sampling_fraction = sampling_fraction,
                                  setting_calculation = NULL,
                                  see_ancestral_states = FALSE,
                                  loglik_penalty = 0,
                                  is_complete_tree = FALSE)

  secsse_cla_LL4 <- secsse_loglik(parameter = parameter,
                                  phy = example_phy_GeoSSE,
                                  traits = traits,
                                  num_concealed_states =
                                    num_concealed_states,
                                  cond = "maddison_cond",
                                  root_state_weight = "maddison_weights",
                                  sampling_fraction = sampling_fraction,
                                  setting_calculation = NULL,
                                  see_ancestral_states = FALSE,
                                  loglik_penalty = 0,
                                  is_complete_tree = TRUE)
  testthat::expect_equal(secsse_cla_LL3, secsse_cla_LL4)

  skip_on_cran()
  secsse_cla_LL5 <- secsse_loglik(parameter = parameter,
                                  phy = example_phy_GeoSSE,
                                  traits = traits,
                                  num_concealed_states =
                                    num_concealed_states,
                                  cond = "maddison_cond",
                                  root_state_weight = "maddison_weights",
                                  sampling_fraction = sampling_fraction,
                                  setting_calculation = NULL,
                                  see_ancestral_states = FALSE,
                                  loglik_penalty = 0,
                                  is_complete_tree = TRUE,
                                  num_threads = 4)
  testthat::expect_equal(secsse_cla_LL5, secsse_cla_LL4,
                         tolerance = 1e-5)

  parenthesis <- "((t1:13.27595158,(((t7:3.890382947,t44:3.890382947):1.853160984,((t28:1.711947644,t52:0.4956923013):1.025240512,t49:2.737188156):3.006355775):8.137718231,t8:0.505931684):0.03852050838):1.080217329,(((((((t2:1.223724296,t54:1.223724296):2.937627297,(t43:1.877801583,t51:1.477270763):2.283550009):0.3267835885,t39:4.488135181):3.978299002,(t20:5.332776925,t33:1.090685514):3.133657257):0.6198399825,(t17:2.592728197,t21:8.418528959):0.6677452056):0.5788113411,((t13:9.543568307,t15:4.657699849):0.03128867016,(((t14:0.2753485556,((t27:1.893882667,t34:4.969412207):0.4876873725,t31:5.45709958):0.2968375929):2.956689195,((t18:3.089806926,t47:3.089806926):3.812406896,(t23:4.616705952,t37:3.696779257):2.28550787):1.808412546):0.6634713591,t16:4.343870947):0.2007592503):0.09022852898):5.130443554,((t3:3.025694309,(((t5:0.6527575809,((t10:8.190240586,t22:4.624901141):1.973824751,((t12:4.230710001,(t42:0.2233137827,t55:0.2233137827):4.007396218):4.263802978,((((t19:4.431551413,t40:4.431551413):1.104239624,t30:0.1129381496):1.083744321,t26:1.989902921):0.2782431807,t24:0.2097131009):1.596734441):1.669552358):1.61638294):1.700092275,((t9:1.444919643,t53:1.444919643):5.416788797,(((t25:4.956186112,(t35:0.07136896428,((t41:2.961601359,(t48:0.04657504123,t56:0.04657504123):2.915026317):0.6168912293,t45:3.578492588):0.7569031841):0.6207903395):0.4454730422,(t32:3.460649902,t46:3.460649902):1.941009252):0.3114551734,t29:4.364985142):1.148594113):6.618832112):0.9318119344,((((t6:2.605426467,t50:0.4317387896):2.002392571,t38:4.607819038):0.207438208,t36:4.815257246):6.619291453,t11:11.4345487):2.977803786):0.1895024879):0.1670130749,t4:0.903839228):0.026661011):0.20447094):0;" # nolint
  phy <- ape::read.tree(file = "", parenthesis)
  traits <- c(0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 1, 1,
              0, 0, 0, 0, 1, 1, 0, 0, 2, 0, 0, 0, 1, 1, 0,
              0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0,
              0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0)
  lambdas <- list()
  lambdas[[1]] <- matrix(0, ncol = 9, nrow = 9, byrow = TRUE)
  lambdas[[1]][2, 1] <- 0.5
  lambdas[[1]][3, 1] <- 0.3
  lambdas[[1]][3, 2] <- 0.2
  for (i in 2:9) {
    lambdas[[i]] <- lambdas[[1]]
  }
  mus <- rep(0.2, 9)
  set.seed(42)
  Q <- matrix(stats::runif(81), ncol = 9, nrow = 9, byrow = TRUE)
  diag(Q) <- 0
  parameter <- list()
  parameter[[1]] <- lambdas
  parameter[[2]] <- mus
  parameter[[3]] <- Q

  num_concealed_states <- 3
  sampling_fraction <- c(1, 1, 1)

  secsse_cla_LL6 <- secsse_loglik(parameter = parameter,
                                  phy = phy,
                                  traits = traits,
                                  num_concealed_states =
                                    num_concealed_states,
                                  cond = "maddison_cond",
                                  root_state_weight = "maddison_weights",
                                  sampling_fraction = sampling_fraction,
                                  setting_calculation = NULL,
                                  see_ancestral_states = FALSE,
                                  loglik_penalty = 0,
                                  is_complete_tree = TRUE)

  # hardcoded LL, don't know where the value comes from!
  # pauze this test until reply from Rampal, seems to be carry-over from
  # accidental coding placement.
  # testthat::expect_equal(secsse_cla_LL6, -439.7388, tol = 1E-3)
})
