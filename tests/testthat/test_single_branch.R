test_that("single branch check", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 2, birth = 0.3 ,death = 0)
  focal_tree$root.edge <- NULL
  traits <- c(1, 1)
  num_concealed_states <- 2
  idparslist <- secsse::cla_id_paramPos(c(1, 2), num_concealed_states)
  idparslist$lambdas[1, ] <- rep(1, 2)
  idparslist[[2]][] <- 2
  masterBlock <- matrix(3, ncol = 2, nrow = 2, byrow = TRUE)
  diag(masterBlock) <- NA
  diff.conceal <- FALSE
  idparslist[[3]] <- secsse::q_doubletrans(c(1, 2), masterBlock, diff.conceal)
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
  
  res1 <- secsse::cla_secsse_loglik(parameter = parslist,
                                    phy = focal_tree,
                                    traits = traits,
                                    num_concealed_states = num_concealed_states,
                                    sampling_fraction = sf,
                                    cond = "no_cond",
                                    display_warning = FALSE)
  
  # create phylogeny with single branch:
  
  phy <- focal_tree
  phy$edge <- phy$edge[-2, ]
  phy$edge.length <- phy$edge.length[-2]
  phy$tip.label <- phy$tip.label[-2]
  traits <- traits[-2]
  
  res2 <- secsse::secsse_single_branch_loglik(parameter = parslist,
                                              phy = phy,
                                              traits = traits,
                                              num_concealed_states = 
                                                num_concealed_states,
                                              sampling_fraction = sf,
                                              cond = "no_cond",
                                              display_warning = FALSE)
  
  d <- length(mus)
  sz <- res2$nodeM[1:(d + d)] * params[1]
  prefact <- log(sum(abs(sz)))
  answ_normal <- (res1$LL - prefact) / 2
  answ_single_branch <- res2$loglik
  testthat::expect_equal(answ_normal, answ_single_branch)
  
  # now check that the single branch loglik is the same as DDD::bd_loglik
  
  secsse_ll <- res2$loglik
  bd_ll <- DDD::bd_loglik(pars1 = c(0.3,0.0),
                          pars2 = c(0,0,0,0,1),
                          brts = phy$edge.length,
                          missnumspec = 0)
  testthat::expect_equal(bd_ll,secsse_ll)
  
  # Do the same with nonzero extinction
  
  params <- c(0.3, 0.1, 0.0)
  lambdas <- secsse::fill_in(idparslist[[1]], params)
  mus <- secsse::fill_in(idparslist[[2]], params)
  q_mat <- secsse::fill_in(idparslist[[3]], params)
  parslist <- list()
  parslist[[1]] <- lambdas
  parslist[[2]] <- mus
  parslist[[3]] <- q_mat
  sf <- c(1, 1)
  res3 <- secsse::secsse_single_branch_loglik(parameter = parslist,
                                              phy = phy,
                                              traits = traits,
                                              num_concealed_states = 
                                                num_concealed_states,
                                              sampling_fraction = sf,
                                              cond = "no_cond",
                                              display_warning = FALSE)
  secsse_ll <- res3$loglik
  bd_ll <- DDD::bd_loglik(pars1 = c(0.3,0.1),
                          pars2 = c(0,0,1,0,1),
                          brts = phy$edge.length,
                          missnumspec = 0)
  testthat::expect_equal(bd_ll,secsse_ll)
  
  # go back to tree with two tips and check loglik with stem age
  
  set.seed(42)
  focal_tree <- ape::rphylo(n = 2, birth = 0.3, death = 0)
  focal_tree$root.edge <- 0.3
  traits <- c(1, 1)
  secsse_ll <- secsse::cla_secsse_loglik(parameter = parslist,
                                    phy = focal_tree,
                                    traits = traits,
                                    num_concealed_states = num_concealed_states,
                                    sampling_fraction = sf,
                                    take_into_account_root_edge = TRUE,
                                    cond = "no_cond",
                                    display_warning = FALSE)
  brts <- ape::branching.times(focal_tree)
  bd_ll <- DDD::bd_loglik(pars1 = c(0.3,0.1),
                     pars2 = c(0,0,1,0,1),
                     brts = c(focal_tree$root.edge + max(brts), brts),
                     missnumspec = 0)
  testthat::expect_equal(bd_ll,secsse_ll$LL)
  
  # now the penultimate test of a tree with three tips and check loglik 
  # with stem age

  set.seed(42)
  focal_tree <- ape::rphylo(n = 3, birth = 0.3, death = 0)
  focal_tree$root.edge <- 0.3
  traits <- c(1, 1, 1)
  secsse_ll <- secsse::cla_secsse_loglik(parameter = parslist,
                                         phy = focal_tree,
                                         traits = traits,
                                         num_concealed_states =
                                           num_concealed_states,
                                         sampling_fraction = sf,
                                         take_into_account_root_edge = TRUE,
                                         cond = "no_cond",
                                         display_warning = FALSE)
  brts <- ape::branching.times(focal_tree)
  bd_ll <- DDD::bd_loglik(pars1 = c(0.3,0.1),
                     pars2 = c(0,0,1,0,1),
                     brts = c(focal_tree$root.edge + max(brts), brts),
                     missnumspec = 0)
  testthat::expect_equal(bd_ll,secsse_ll$LL)
  
  # now the ultimate test of a tree with four tips and check loglik
  # with stem age
  
  set.seed(42)
  focal_tree <- ape::rphylo(n = 4, birth = 0.3, death = 0)
  focal_tree$root.edge <- 0.3
  traits <- c(1, 1, 1, 1)
  secsse_ll <- secsse::cla_secsse_loglik(parameter = parslist,
                                         phy = focal_tree,
                                         traits = traits,
                                         num_concealed_states =
                                            num_concealed_states,
                                         sampling_fraction = sf,
                                         take_into_account_root_edge = TRUE,
                                         cond = "no_cond",
                                         display_warning = FALSE)
  brts <- ape::branching.times(focal_tree)
  bd_ll <- DDD::bd_loglik(pars1 = c(0.3,0.1),
                          pars2 = c(0,0,1,0,1),
                          brts = c(focal_tree$root.edge + max(brts), brts),
                          missnumspec = 0)
  testthat::expect_equal(bd_ll,secsse_ll$LL)
  
  # Now try on a multi-phylo data set
  
  set.seed(42)
  focal_tree <- list()
  class(focal_tree) <- "multiPhylo"
  focal_tree[[1]] <- ape::rphylo(n = 4, birth = 0.3, death = 0)
  focal_tree[[2]] <- ape::rphylo(n = 4, birth = 0.3, death = 0)
  focal_tree[[1]]$root.edge <- 0.3
  focal_tree[[2]]$root.edge <- 0.2
  traits <- list()
  traits[[1]] <- c(1, 1, 1, 1)
  traits[[2]] <- c(1, 1, 1, 1)
  secsse_ll <- secsse::cla_secsse_loglik(parameter = parslist,
                                         phy = focal_tree,
                                         traits = traits,
                                         num_concealed_states =
                                           num_concealed_states,
                                         sampling_fraction = sf,
                                         take_into_account_root_edge = TRUE,
                                         cond = "no_cond",
                                         display_warning = FALSE)
  bd_ll <- 0
  for(i in 1:2) {
    brts <- ape::branching.times(focal_tree[[i]])
    bd_ll <- bd_ll + DDD::bd_loglik(pars1 = c(0.3,0.1),
                          pars2 = c(0,0,1,0,1),
                          brts = c(focal_tree[[i]]$root.edge + max(brts), brts),
                          missnumspec = 0)
  }
  testthat::expect_equal(bd_ll,secsse_ll)
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
  
  res1 <- secsse::cla_secsse_loglik(parameter = parslist,
                                    phy = focal_tree,
                                    traits = traits,
                                    num_concealed_states =
                                      num_concealed_states,
                                    sampling_fraction = sf,
                                    cond = "no_cond",
                                    display_warning = FALSE)
  # create phylogeny with a root branch:
  
  focal_tree$root.edge <- 3
  res2 <- secsse::cla_secsse_loglik(parameter = parslist,
                                    phy = focal_tree,
                                    traits = traits,
                                    num_concealed_states =
                                      num_concealed_states,
                                    sampling_fraction = sf,
                                    cond = "no_cond",
                                    display_warning = FALSE)
  
  res3 <- secsse::cla_secsse_loglik(parameter = parslist,
                                    phy = focal_tree,
                                    traits = traits,
                                    num_concealed_states =
                                      num_concealed_states,
                                    sampling_fraction = sf,
                                    cond = "no_cond",
                                    take_into_account_root_edge = TRUE,
                                    display_warning = FALSE)

  testthat::expect_equal(res1, res2)
  testthat::expect_lt(res3$LL, res1$LL)
})
