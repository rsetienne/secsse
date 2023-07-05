context("lambda_and_qmat_setup")

test_that("lambda setup", {
  # Islandness, ETD model
  full_lambdas <- list()
  
  for (i in 1:6) {
    full_lambdas[[i]] <- matrix(0, 6, 6)
    colnames(full_lambdas[[i]]) <- c("MA", "IA", "CA", "MB", "IB", "CB")
    rownames(full_lambdas[[i]]) <- c("MA", "IA", "CA", "MB", "IB", "CB")
  }
  
  full_lambdas[[1]][1, 1] <- 1  # MA, lambda_Mainland_sympatric
  full_lambdas[[2]][2, 2] <- 2  # IA, lambda_Island_sympatric
  full_lambdas[[3]][1, 2] <- 3  # CA, lambda_CA->MA,IA
  full_lambdas[[3]][2, 1] <- 3  # CA, lambda_CA->MA,IA
  
  full_lambdas[[4]][4, 4] <- 1  # MB,  lambda_Mainland_sympatric
  full_lambdas[[5]][5, 5] <- 2  # IB, lambda_Island_sympatric
  full_lambdas[[6]][4, 5] <- 3  # CB, lambda_CB->MB,IB
  full_lambdas[[6]][5, 4] <- 3  # CB, lambda_CB->MB,IB
  
  states <- c("M", "I", "C")
  
  transition_matrix <- c()
  transition_matrix <- rbind(transition_matrix, c("M", "M", "M", 1))
  transition_matrix <- rbind(transition_matrix, c("I", "I", "I", 2))
  transition_matrix <- rbind(transition_matrix, c("C", "M", "I", 3))
  
  lambdas <- secsse::create_lambda_list(state_names = states,
                                        num_concealed_states = 2,
                                        transition_matrix = transition_matrix)
  
  testthat::expect_equal(length(lambdas), length(full_lambdas))
  for (i in seq_along(lambdas)) {
    testthat::expect_equal(lambdas[[i]], full_lambdas[[i]])
  }
})

test_that("q_matrix", {
  q_mat <- matrix(data = NA, nrow = 2, ncol = 2)
  q_mat[1, 2] <- 1
  q_mat[2, 1] <- 2
  
  # first, we test on a 2x2 matrix
  for (dd in c(TRUE, FALSE)) {
    q1 <- secsse::q_doubletrans(traits = c(1, 2),
                                masterBlock = q_mat,
                                diff.conceal = dd)
    
    q2 <- secsse::expand_q_matrix(q_matrix = q_mat,
                                  num_concealed_states = 2,
                                  diff.conceal = dd)
    v1 <- as.vector(q1)
    v2 <- as.vector(q2)
    testthat::expect_true(all.equal(v1, v2))
  }
})

test_that("setup", {
  focal_matrix <- 
    secsse::create_default_lambda_transition_matrix(state_names = c("S", "N"),
                                                    model = "CR")
  lambda_list_CR <- secsse::create_lambda_list(state_names = c("S", "N"),
                                               num_concealed_states = 2,
                                               transition_matrix = focal_matrix,
                                               model = "CR")
  
  for (i in 1:4) {
    testthat::expect_equal(lambda_list_CR[[i]][i, i], 1)
  }
  
  focal_matrix <-
    secsse::create_default_lambda_transition_matrix(state_names = c("S", "N"),
                                                    model = "CTD")
  # now for the CTD model:
  lambda_list_CTD <- secsse::create_lambda_list(state_names = c("S", "N"),
                                                num_concealed_states = 2,
                                              transition_matrix = focal_matrix,
                                                model = "CTD")
  
  for (i in 1:4) {
    testthat::expect_equal(lambda_list_CTD[[i]][i, i], ceiling(i / 2))
  }
  
  # and the ETD model:
  lambda_list_ETD <- secsse::create_lambda_list(state_names = c("S", "N"),
                                                num_concealed_states = 2,
                                              transition_matrix = focal_matrix,
                                                model = "ETD")
  
  for (i in 1:4) {
    testthat::expect_equal(lambda_list_ETD[[i]][i, i], 2 - i %% 2)
  }
  
  # and now the mu vector
  mus_CR <- secsse::create_mu_vector(state_names = c("S", "N"),
                                     num_concealed_states = 2,
                                     model = "CR",
                                     lambda_list = lambda_list_CR)
  for (i in 1:4) {
    testthat::expect_equal(mus_CR[[i]], 2)
  }
  
  mus_CTD <- secsse::create_mu_vector(state_names = c("S", "N"),
                                      num_concealed_states = 2,
                                      model = "CTD",
                                      lambda_list = lambda_list_CTD)
  for (i in 1:4) {
    testthat::expect_equal(mus_CTD[[i]], 3 + floor(i / 3))
  }
  
  mus_ETD <- secsse::create_mu_vector(state_names = c("S", "N"),
                                      num_concealed_states = 2,
                                      model = "ETD",
                                      lambda_list = lambda_list_ETD)
  for (i in 1:4) {
    testthat::expect_equal(mus_ETD[[i]], 4 - i %% 2)
  }
  
  # and the q matrices
  t_CR <- secsse::create_default_shift_matrix(state_names = c("S", "N"),
                                              num_concealed_states = 2,
                                              mus = mus_CR)
  q_CR <- secsse::create_q_matrix(state_names = c("S", "N"),
                                  num_concealed_states = 2,
                                  shift_matrix = t_CR,
                                  diff.conceal = TRUE)
  testthat::expect_equal(6, max(q_CR, na.rm = TRUE))
  
  t_CTD <- secsse::create_default_shift_matrix(state_names = c("S", "N"),
                                               num_concealed_states = 2,
                                               mus = mus_CTD)
  q_CTD <- secsse::create_q_matrix(state_names = c("S", "N"),
                                   num_concealed_states = 2,
                                   shift_matrix = t_CTD,
                                   diff.conceal = TRUE)
  
  testthat::expect_equal(8, max(q_CTD, na.rm = TRUE))
  
  t_ETD <- secsse::create_default_shift_matrix(state_names = c("S", "N"),
                                               num_concealed_states = 2,
                                               mus = mus_ETD)
  q_ETD <- secsse::create_q_matrix(state_names = c("S", "N"),
                                   num_concealed_states = 2,
                                   shift_matrix = t_ETD,
                                   diff.conceal = TRUE)
  testthat::expect_equal(8, max(q_ETD, na.rm = TRUE))
})
