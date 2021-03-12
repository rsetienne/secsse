context("improve speed")


test_that("test cla_loglik_cpp", {
  phy <- NULL; rm(phy);
  utils::data('example_phy_GeoSSE', package = 'secsse');
  traits <- as.numeric(phy$tip.state)
  lambdas <- list()
  lambdas[[1]] <- matrix(0,ncol = 9,nrow = 9,byrow = TRUE)
  lambdas[[1]][2,1] <- 1.5
  lambdas[[1]][3,1] <- 0.5
  lambdas[[1]][3,2] <- 1
  for (i in 2:9) {
    lambdas[[i]] <- lambdas[[1]]
  }
  mus <- rep(0,9)
  Q <- matrix(stats::runif(81),ncol = 9,nrow = 9,byrow = TRUE)
  #diag(Q) <- NA
  parameter <- list()
  parameter[[1]] <- lambdas
  parameter[[2]] <- mus
  parameter[[3]] <- Q

  num_concealed_states <- 3
  sampling_fraction <- c(1,1,1)

  secsse_cla_LL_CPP <- secsseCPP::cla_secsse_loglik_cpp(parameter = parameter,
                                          phy = phy,
                                          traits = traits,
                                          num_concealed_states = num_concealed_states,
                                          use_fortran = FALSE,
                                          methode = "ode45",
                                          cond = "maddison_cond",
                                          root_state_weight = "maddison_weights",
                                          sampling_fraction = sampling_fraction,
                                          run_parallel = FALSE,
                                          setting_calculation = NULL,
                                          setting_parallel = NULL,
                                          see_ancestral_states = FALSE,
                                          loglik_penalty = 0,
                                          is_complete_tree = FALSE)

  secsse_cla_LL_R <- secsse::cla_secsse_loglik(parameter = parameter,
                                      phy = phy,
                                      traits = traits,
                                      num_concealed_states = num_concealed_states,
                                      use_fortran = TRUE,
                                      methode = "ode45",
                                      cond = "maddison_cond",
                                      root_state_weight = "maddison_weights",
                                      sampling_fraction = sampling_fraction,
                                      run_parallel = FALSE,
                                      setting_calculation = NULL,
                                      setting_parallel = NULL,
                                      see_ancestral_states = FALSE,
                                      loglik_penalty = 0,
                                      is_complete_tree = FALSE)
  secsse_cla_LL_R
  secsse_cla_LL_CPP


  testthat::expect_equal(secsse_cla_LL_R, secsse_cla_LL_CPP, tolerance = 0.001)

  r_func <- function() {
    secsse::cla_secsse_loglik(parameter = parameter,
                      phy = phy,
                      traits = traits,
                      num_concealed_states = num_concealed_states,
                      use_fortran = TRUE,
                      methode = "ode45",
                      cond = "maddison_cond",
                      root_state_weight = "maddison_weights",
                      sampling_fraction = sampling_fraction,
                      run_parallel = FALSE,
                      setting_calculation = NULL,
                      setting_parallel = NULL,
                      see_ancestral_states = FALSE,
                      loglik_penalty = 0,
                      is_complete_tree = FALSE)
  }

  cpp_func <- function() {
    cla_secsse_loglik_cpp(parameter = parameter,
                          phy = phy,
                          traits = traits,
                          num_concealed_states = num_concealed_states,
                          use_fortran = FALSE,
                          methode = "ode45",
                          cond = "maddison_cond",
                          root_state_weight = "maddison_weights",
                          sampling_fraction = sampling_fraction,
                          run_parallel = FALSE,
                          setting_calculation = NULL,
                          setting_parallel = NULL,
                          see_ancestral_states = FALSE,
                          loglik_penalty = 0,
                          is_complete_tree = FALSE)
  }

  testthat::expect_equal(r_func(), cpp_func(), tolerance = 0.001)

  repl <- 100
  timez <- c()
  for (i in 1:repl) {
    t0 <- Sys.time()
    cpp_func()
    t1 <- Sys.time()
    timez[i] <- difftime(t1, t0, units = "secs")[[1]]
  }

  repl <- 3
  timez2 <- c()
  for (i in 1:repl) {
    t0 <- Sys.time()
    r_func()
    t1 <- Sys.time()
    timez2[i] <- difftime(t1, t0, units = "secs")[[1]]
  }

  cat("mean time cpp: ", mean(timez), "\n")
  cat("mean time r  : ", mean(timez2), "\n")
  cat("speedup :", mean(timez2) / mean(timez), "\n")
})

test_that("secsse new gives same result as secsse old", {

  set.seed(42)
  out <- DDD::dd_sim(pars = c(0.4, 0.1, 40), age = 15)
  # out <- DDD::dd_sim(pars = c(0.4, 0.1, 4000), age = 25)
  phy <- out$tes
  traits <- sample(c(0,1), ape::Ntip(phy),replace = T)
  b <- c(0.04, 0.04)  # lambda
  d <- rep(0, 2)
  userTransRate <- 0.2 # transition rate among trait states
  num_concealed_states <- 2
  sampling_fraction <- c(1, 1)
  toCheck <- secsse::id_paramPos(traits, num_concealed_states)
  toCheck[[1]][] <- b
  toCheck[[2]][] <- d
  toCheck[[3]][,] <- userTransRate
  diag(toCheck[[3]]) <- NA
  root_state_weight <- "maddison_weights"
  cond <- "noCondit"

  fortran_func <- function() {
    secsse::secsse_loglik(parameter = toCheck,
                          phy = phy,
                          traits = traits,
                          num_concealed_states = num_concealed_states,
                          use_fortran = TRUE,
                          methode = "ode45",
                          cond = cond,
                          root_state_weight = root_state_weight,
                          sampling_fraction = sampling_fraction,
                          is_complete_tree = TRUE,
                          func = "secsse_runmod_ct")
  }

  cpp_func <- function() {
    secsseCPP::secsse_loglik_cpp(parameter = toCheck,
                                 phy = phy,
                                 traits = traits,
                                 num_concealed_states = num_concealed_states,
                                 cond = cond,
                                 root_state_weight = root_state_weight,
                                 sampling_fraction = sampling_fraction,
                                 is_complete_tree = TRUE)
  }
  testthat::expect_equal(fortran_func(), cpp_func())

  vv <- microbenchmark::microbenchmark(fortran_func,
                                       cpp_func,
                                       times = 10000)

  # require(ggplot2)
  # ggplot2::autoplot(vv)
})

