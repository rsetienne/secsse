test_that("normal plotting", {
   set.seed(5)
   phy <- ape::rphylo(n = 4, birth = 1, death = 0)
   traits <- c(0, 1, 1, 0)
   params <- secsse::id_paramPos(c(0, 1), 2)
   params[[1]][] <- c(0.2, 0.2, 0.1, 0.1)
   params[[2]][] <- 0.01
   params[[3]][, ] <- 0.1
   diag(params[[3]]) <- NA
   #  Thus, we have for both, rates
   # 0A, 1A, 0B and 1B. If we are interested in the posterior probability of
   # trait 0,we have to provide a helper function that sums the probabilities of
   # 0A and 0B, e.g.:
   helper_function <- function(x) {
     return(sum(x[c(5, 7)]) / sum(x)) # normalized by total sum, just in case.
   }
  testthat::expect_warning(
   px <- plot_state_exact(parameters = params,
                    phy = phy,
                    traits = traits,
                    num_concealed_states = 2,
                    sampling_fraction = c(1, 1),
                    num_steps = 10,
                    prob_func = helper_function)
  )
   testthat::expect_true(inherits(px, "ggplot"))
})

test_that("cla plotting", {
  skip_on_cran()
  parenthesis <- "(((6:0.2547423371,(1:0.0496153503,4:0.0496153503):0.2051269868):0.1306304758,(9:0.2124135406,5:0.2124135406):0.1729592723):1.151205247,(((7:0.009347664296,3:0.009347664296):0.2101416075,10:0.2194892718):0.1035186448,(2:0.2575886319,8:0.2575886319):0.06541928469):1.213570144);" #nolint
  phylotree <- ape::read.tree(file = "", parenthesis)
  traits <- c(2, 0, 1, 0, 2, 0, 1, 2, 2, 0)
  num_concealed_states <- 3
  idparslist <- cla_id_paramPos(traits, num_concealed_states)
  idparslist$lambdas[2, ] <- rep(1, 9)
  idparslist[[2]][] <- 4
  masterBlock <- matrix(5, ncol = 3, nrow = 3, byrow = TRUE)
  diag(masterBlock) <- NA
  diff.conceal <- FALSE

  idparslist[[3]] <- q_doubletrans(traits, masterBlock, diff.conceal)

  testthat::expect_output(
    startingpoint <- DDD::bd_ML(brts = ape::branching.times(phylotree))
  )
  intGuessLamba <- startingpoint$lambda0
  intGuessMu <- startingpoint$mu0
  idparsopt <- c(1)
  initparsopt <- c(rep(intGuessLamba, 1))
  idparsfix <- c(0, 4, 5)
  parsfix <- c(0, 0, 0.01)
  tol <- c(1e-04, 1e-05, 1e-07)
  maxiter <- 1000 * round((1.25) ^ length(idparsopt))
  optimmethod <- "simplex"
  cond <- "proper_cond"
  root_state_weight <- "proper_weights"
  sampling_fraction <- c(1, 1, 1)

  testthat::expect_warning(
    model_R <- cla_secsse_ml(
      phy = phylotree,
      traits = traits,
      num_concealed_states = num_concealed_states,
      idparslist = idparslist,
      idparsopt = idparsopt,
      initparsopt = initparsopt,
      idparsfix = idparsfix,
      parsfix = parsfix,
      cond = cond,
      root_state_weight = root_state_weight,
      sampling_fraction = sampling_fraction,
      tol = tol,
      maxiter = maxiter,
      optimmethod = optimmethod,
      num_cycles = 1,
      verbose = FALSE)
  )

  helper_function <- function(x) {
    return(sum(x[c(10, 13, 16)]) / sum(x))
  }

  testthat::expect_warning(
    px <- secsse::plot_state_exact(parameters = model_R$MLpars,
                                   phy = phylotree,
                                   traits = traits,
                                   num_concealed_states =
                                   num_concealed_states,
                                   sampling_fraction = sampling_fraction,
                                   cond = cond,
                                   root_state_weight = root_state_weight,
                                   prob_func = helper_function)
  )

  testthat::expect_true(inherits(px, "ggplot"))
})
