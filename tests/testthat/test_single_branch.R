test_that("single branch check", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 2, birth = 0.3 ,death = 0)
  
  focal_tree$root.edge <- NULL
  traits <- c(1, 1)
  
  num_concealed_states <- 2
  idparslist <- cla_id_paramPos(c(1, 2), num_concealed_states)
  idparslist$lambdas[1, ] <- rep(1, 2)
  idparslist[[2]][] <- 2
  masterBlock <- matrix(3, ncol = 2, nrow = 2, byrow = TRUE)
  diag(masterBlock) <- NA
  diff.conceal <- FALSE
  idparslist[[3]] <- q_doubletrans(c(1, 2), masterBlock, diff.conceal)
  idparslist[[1]] <- secsse::prepare_full_lambdas(c(1, 2),
                                                  num_concealed_states,
                                                  idparslist[[1]])
  
 # Expect warning because some transitions are set to be impossible

  params <- c(0.3, 0.0, 0.0) # extinction and shifts to zero, to allow direct
                             # comparison
  
  lambdas <- secsse::fill_in(idparslist[[1]], params)
  mus <- secsse::fill_in(idparslist[[2]], params)
  q_mat <- secsse::fill_in(idparslist[[3]], params)
  
  parslist <- list()
  parslist[[1]] <- lambdas
  parslist[[2]] <- mus
  parslist[[3]] <- q_mat
  
  sf <- c(1, 1)
  
  testthat::expect_warning(
    res1 <- secsse::cla_secsse_loglik(parameter = parslist,
                                    phy = focal_tree,
                                    traits = traits,
                                    num_concealed_states = num_concealed_states,
                                    sampling_fraction = sf,
                                    cond = "no_cond")
  )
  # create phylogeny with single branch:
  
  phy <- focal_tree
  phy$edge <- phy$edge[-2, ]
  phy$edge.length <- phy$edge.length[-2]
  phy$tip.label <- phy$tip.label[-2]
  traits <- traits[-2]
  
  testthat::expect_warning(
  res2 <- secsse::secsse_single_branch_loglik(parameter = parslist,
                                              phy = phy,
                                              traits = traits,
                                              num_concealed_states = 
                                                num_concealed_states,
                                              sampling_fraction = sf,
                                              cond = "no_cond")
  )
  sz <- res2$nodeM * params[1]
  prefact <- log(sum(abs(sz)))
  answ_normal <- (res1 - prefact) / 2
  answ_single_branch <- res2$loglik
  testthat::expect_equal(answ_normal, answ_single_branch)
})

test_that("root branch check", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 2, birth = 0.3 ,death = 0)
  
  focal_tree$root.edge <- NULL
  traits <- c(1, 1)
  
  num_concealed_states <- 2
  idparslist <- cla_id_paramPos(c(1, 2), num_concealed_states)
  idparslist$lambdas[1, ] <- rep(1, 2)
  idparslist[[2]][] <- 2
  masterBlock <- matrix(3, ncol = 2, nrow = 2, byrow = TRUE)
  diag(masterBlock) <- NA
  diff.conceal <- FALSE
  idparslist[[3]] <- q_doubletrans(c(1, 2), masterBlock, diff.conceal)
  idparslist[[1]] <- secsse::prepare_full_lambdas(c(1, 2),
                                                  num_concealed_states,
                                                  idparslist[[1]])
  
  # Expect warning because some transitions are set to be impossible
  
  params <- c(0.3, 0.0, 0.0) # extinction and shifts to zero, to allow direct
  # comparison
  
  lambdas <- secsse::fill_in(idparslist[[1]], params)
  mus <- secsse::fill_in(idparslist[[2]], params)
  q_mat <- secsse::fill_in(idparslist[[3]], params)
  
  parslist <- list()
  parslist[[1]] <- lambdas
  parslist[[2]] <- mus
  parslist[[3]] <- q_mat
  
  sf <- c(1, 1)
  
  testthat::expect_warning(
    res1 <- secsse::cla_secsse_loglik(parameter = parslist,
                                      phy = focal_tree,
                                      traits = traits,
                                      num_concealed_states =
                                        num_concealed_states,
                                      sampling_fraction = sf,
                                      cond = "no_cond")
  )
  # create phylogeny with a root branch:
  
  focal_tree$root.edge <- 3
  testthat::expect_warning(
    res2 <- secsse::cla_secsse_loglik(parameter = parslist,
                                      phy = focal_tree,
                                      traits = traits,
                                      num_concealed_states =
                                        num_concealed_states,
                                      sampling_fraction = sf,
                                      cond = "no_cond")
  )
  
  testthat::expect_warning(
    res3 <- secsse::cla_secsse_loglik(parameter = parslist,
                                      phy = focal_tree,
                                      traits = traits,
                                      num_concealed_states =
                                        num_concealed_states,
                                      sampling_fraction = sf,
                                      cond = "no_cond",
                                      take_into_account_root_edge = TRUE)
  )
  
  testthat::expect_equal(res1, res2)
  testthat::expect_lt(res3, res1)
})