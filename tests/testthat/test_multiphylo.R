test_that("multi phylo", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 3, numbsim = 1, 
                                     lambda = 1, mu = 0)[[1]]
  focal_tree$root.edge <- NULL
  traits <- c(1, 1, 1)
  
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
  
  focal_tree$root.edge <- NULL
  testthat::expect_warning(
  res1 <- secsse::cla_secsse_loglik(parameter = parslist,
                                    phy = focal_tree,
                                    traits = traits,
                                    num_concealed_states = num_concealed_states,
                                    sampling_fraction = sf,
                                    cond = "no_cond")
  )
  trees <- list()
  trees[[1]]<- focal_tree
  trees[[2]]<- focal_tree
  
  class(trees) <- "multiPhylo"
  
  trait_list <- list()
  trait_list[[1]] <- traits
  trait_list[[2]] <- traits
  testthat::expect_warning(
  res2 <- secsse::cla_secsse_loglik(parameter = parslist,
                                    phy = trees,
                                    traits = trait_list,
                                    num_concealed_states = num_concealed_states,
                                    sampling_fraction = sf,
                                    cond = "no_cond")
  )
  testthat::expect_equal(2*res1, res2)
})