test_that("build states basic", {
  phy <- ape::rphylo(3, 1, 0)
  traits <- c(1, 2, 2)
  
  
  res <- secsse:::build_states(phy = phy,
                               traits = traits,
                               num_concealed_states = 2,
                               sampling_fraction = c(1, 1),
                               is_complete_tree = FALSE)
  testthat::expect_equal(dim(res), c(5, 12))
  E_part <- res[1:3, 1:4]
  testthat::expect_equal(sum(E_part), 0)
  D_part <- res[1:3, 5:8]
  testthat::expect_equal(D_part[1, ], c(1, 0, 1, 0))
  testthat::expect_equal(D_part[2, ], c(0, 1, 0, 1))
  testthat::expect_equal(D_part[3, ], c(0, 1, 0, 1))
})

test_that("build states multi info", {
  phy <- ape::rphylo(3, 1, 0)
  traits <- c(1, 2, 2)
  traits2 <- c(2, 2, 2)
  traits <- cbind(traits, traits2)
  
  
  res <- secsse:::build_states(phy = phy,
                               traits = traits,
                               num_concealed_states = 2,
                               sampling_fraction = c(1, 1),
                               is_complete_tree = FALSE)
  testthat::expect_equal(dim(res), c(5, 12))
  E_part <- res[1:3, 1:4]
  testthat::expect_equal(sum(E_part), 0)
  D_part <- res[1:3, 5:8]
  testthat::expect_equal(D_part[1, ], c(1, 1, 1, 1))
  testthat::expect_equal(D_part[2, ], c(0, 1, 0, 1))
  testthat::expect_equal(D_part[3, ], c(0, 1, 0, 1))
})

test_that("build states NAs", {
  phy <- ape::rphylo(3, 1, 0)
  traits <- c(1, NA, 2)
  
  
  res <- secsse:::build_states(phy = phy,
                               traits = traits,
                               num_concealed_states = 2,
                               sampling_fraction = c(1, 1),
                               is_complete_tree = FALSE)
  testthat::expect_equal(dim(res), c(5, 12))
  E_part <- res[1:3, 1:4]
  testthat::expect_equal(sum(E_part), 0)
  D_part <- res[1:3, 5:8]
  testthat::expect_equal(D_part[1, ], c(1, 0, 1, 0))
  testthat::expect_equal(D_part[2, ], c(1, 1, 1, 1))
  testthat::expect_equal(D_part[3, ], c(0, 1, 0, 1))
  
  
  
  sf <- c(0.1, 0.5)
  res <- secsse:::build_states(phy = phy,
                               traits = traits,
                               num_concealed_states = 2,
                               sampling_fraction = sf,
                               is_complete_tree = FALSE)
  testthat::expect_equal(dim(res), c(5, 12))
  E_part <- res[1:3, 1:4]
  testthat::expect_equal(sum(E_part), sum((1-sf) * 3 * 2))
  D_part <- res[1:3, 5:8]
  testthat::expect_equal(D_part[1, ], c(sf[1], 0, sf[1], 0))
  testthat::expect_equal(D_part[2, ], c(sf, sf))
  testthat::expect_equal(D_part[3, ], c(0, sf[2], 0, sf[2]))
  
  S_part <- res[1:3, 9:12]
  testthat::expect_equal(E_part, 1 - S_part)
})
