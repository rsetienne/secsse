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
  
  lambdas <- secsse::create_lambda_matrices(state_names = states,
                                            num_concealed_states = 2,
                                            transition_list = transition_matrix)
  
  testthat::expect_equal(length(lambdas), length(full_lambdas))
  for (i in seq_along(lambdas)) {
    testthat::expect_equal(lambdas[[i]], full_lambdas[[i]])
  }
})

test_that("qmat setup", {
  q_mat <- matrix(0, nrow = 6, ncol = 6)
  diag(q_mat) <- NA
  q_mat[1, 3] <- 5 # mu_I
  q_mat[1, 4] <- 7 # Q_BA
  q_mat[2, 3] <- 4 # mu_m 
  q_mat[2, 5] <- 7 # Q_BA
  q_mat[3, 1] <- 8 # gamma
  q_mat[3, 2] <- 9 # delta
  q_mat[3, 6] <- 7 # Q_BA
  q_mat[4, 1] <- 6 # Q_AB
  q_mat[4, 6] <- 5 # mu_I
  q_mat[5, 2] <- 6 # Q_AB
  q_mat[5, 6] <- 4 # mu_m
  q_mat[6, 3] <- 6 # Q_AB
  q_mat[6, 4] <- 8 # gamma
  q_mat[6, 5] <- 9 # delta
  
  namez <- c("MA", "IA", "CA", "MB", "IB", "CB")
  
  colnames(q_mat) <- namez
  rownames(q_mat) <- namez
  
  q_mat <- t(q_mat)
  
  trans_list <- c()
  trans_list <- rbind(trans_list, c("CA", "MA", 5))
  trans_list <- rbind(trans_list, c("CA", "IA", 4))
  trans_list <- rbind(trans_list, c("CA", "CB", 6))
  
  trans_list <- rbind(trans_list, c("CB", "MB", 5))
  trans_list <- rbind(trans_list, c("CB", "IB", 4))
  trans_list <- rbind(trans_list, c("CB", "CA", 7))
  
  trans_list <- rbind(trans_list, c("MA", "CA", 8))
  trans_list <- rbind(trans_list, c("MA", "MB", 6))
  
  trans_list <- rbind(trans_list, c("MB", "CB", 8))
  trans_list <- rbind(trans_list, c("MB", "MA", 7))
  
  trans_list <- rbind(trans_list, c("IA", "CA", 9))
  trans_list <- rbind(trans_list, c("IA", "IB", 6))
  trans_list <- rbind(trans_list, c("IB", "CB", 9))
  trans_list <- rbind(trans_list, c("IB", "IA", 7))
  
  q_mat2 <- secsse::create_transition_matrix(namez, trans_list)
  
  testthat::expect_equal(nrow(q_mat), nrow(q_mat2))
  testthat::expect_true(all.equal(q_mat, q_mat2))
})

test_that("q_matrix", {
  
  q_mat <- matrix(data = NA, nrow = 2, ncol = 2)
  q_mat[1, 2] <- 1
  q_mat[2, 1] <- 2
  
  # first, we test on a 2x2 matrix
  for (dd in c(TRUE, FALSE)) {
    testthat::expect_output(
    q1 <- secsse::q_doubletrans(traits = c(1, 2),
                                masterBlock = q_mat,
                                diff.conceal = dd))
    q2 <- secsse::expand_q_matrix(q_matrix = q_mat,
                                  num_concealed_states = 2,
                                  diff.conceal = dd)
    testthat::expect_true(all.equal(q1, q2))
  }
})