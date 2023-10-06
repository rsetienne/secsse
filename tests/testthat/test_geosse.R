test_that("secsse gives the same result as GeoSSE", {
  Sys.unsetenv("R_TESTS")

  if (requireNamespace("diversitree")) {
    #geosse
    pars <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
    names(pars) <- c("sA", "sB", "sAB", "xA", "xB", "dA", "dB")
    utils::data("example_phy_GeoSSE", package = "secsse")
    traits <- as.numeric(example_phy_GeoSSE$tip.state)
    lik.g <- diversitree::make.geosse(example_phy_GeoSSE,
                                      example_phy_GeoSSE$tip.state)
    pars.g <- c(1.5, 0.5, 1.0, 0.7, 0.7, 1.4, 1.3)
    names(pars.g) <- diversitree::argnames(lik.g)
    lik.c <- diversitree::make.classe(example_phy_GeoSSE,
                                      example_phy_GeoSSE$tip.state + 1, 3)
    pars.c <- 0 * diversitree::starting.point.classe(example_phy_GeoSSE, 3)
    pars.c["lambda222"] <- pars.c["lambda112"] <- pars.g["sA"]
    pars.c["lambda333"] <- pars.c["lambda113"] <- pars.g["sB"]
    pars.c["lambda123"] <- pars.g["sAB"]
    pars.c["mu2"] <- pars.c["q13"] <- pars.g["xA"]
    pars.c["mu3"] <- pars.c["q12"] <- pars.g["xB"]
    pars.c["q21"] <- pars.g["dA"]
    pars.c["q31"] <- pars.g["dB"]
    lik.g(pars.g) # -175.7685
    classe_diversitree_LL <- lik.c(pars.c) # -175.7685

    ## Secsse part
    lambdas <- list()
    lambdas[[1]] <- matrix(0, ncol = 3, nrow = 3, byrow = TRUE)
    lambdas[[1]][2, 1] <- 1.5
    lambdas[[1]][3, 1] <- 0.5
    lambdas[[1]][3, 2] <- 1
    lambdas[[2]] <- matrix(0, ncol = 3, nrow = 3, byrow = TRUE)
    lambdas[[2]][2, 2] <- 1.5
    lambdas[[3]] <- matrix(0,
                           ncol = 3,
                           nrow = 3,
                           byrow = TRUE)
    lambdas[[3]][3, 3] <- 0.5

    mus <- c(0, 0.7, 0.7)

    q <- matrix(0, ncol = 3, nrow = 3, byrow = TRUE)
    q[2, 1] <- 1.4
    q[3, 1] <- 1.3
    q[1, 2] <- 0.7
    q[1, 3] <- 0.7

    parameter <- list()
    parameter[[1]] <- lambdas
    parameter[[2]] <- mus
    parameter[[3]] <- q

    num_concealed_states <- 3

    num_modeled_traits <- ncol(q) / floor(num_concealed_states)

    setting_calculation <- build_initStates_time(example_phy_GeoSSE,
                                                 traits,
                                                 num_concealed_states,
                                                 sampling_fraction = c(1, 1, 1),
                                                 is_complete_tree = FALSE,
                                                 mus,
                                                 num_modeled_traits,
                                                 first_time = TRUE)
    states <- setting_calculation$states
    d <- ncol(states) / 2
    new_states <- states[, c(1, 2, 3, 10, 11, 12)]
    states <- new_states
    
    setting_calculation$states <-
         states

    # -191.9567
    secsse_cla_LL <- secsse_loglik(parameter,
                                   example_phy_GeoSSE,
                                   traits,
                                   num_concealed_states,
                                   cond = "maddison_cond",
                                   root_state_weight = "maddison_weights",
                                   sampling_fraction = c(1, 1, 1),
                                   setting_calculation = setting_calculation,
                                   see_ancestral_states = FALSE,
                                   loglik_penalty = 0)

    testthat::expect_equal(classe_diversitree_LL,  secsse_cla_LL,
                           tolerance = 1e-5)

    # Parallel code doesn't work on CI
    testthat::skip_on_cran()
    secsse_cla_LL3 <- secsse_loglik(parameter,
                                    example_phy_GeoSSE,
                                    traits,
                                    num_concealed_states,
                                    cond = "maddison_cond",
                                    root_state_weight = "maddison_weights",
                                    sampling_fraction = c(1, 1, 1),
                                    setting_calculation = setting_calculation,
                                    see_ancestral_states = FALSE,
                                    loglik_penalty = 0,
                                    num_threads = 4)
    testthat::expect_equal(classe_diversitree_LL, secsse_cla_LL3,
                           tolerance = 1e-5)
  }
})
