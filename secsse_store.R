library(secsse)
library(RcppParallel)

rm(list = ls())
set.seed(42)
#set.seed(51)
out <- DDD::dd_sim(pars = c(0.5 , 0.3, 10000), age = 50)
num_steps = 100

phy <- out$tes
cat("this tree has:", phy$Nnode + 1, "tips and", phy$Nnode, "internal nodes, num_steps =", num_steps, "\n")

num_concealed_states <- 3

traits <- sample(c(0,1, 2), ape::Ntip(phy),replace = TRUE)

sampling_fraction = c(1, 1, 1)
idparlist <- cla_id_paramPos(traits, num_concealed_states)
lambda_and_modeSpe <- idparlist$lambdas
lambda_and_modeSpe[1,] <- c(0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 0.01, 0.01, 0.01)

parameter <- list()
parameter[[1]] <- prepare_full_lambdas(traits, num_concealed_states,
                                       lambda_and_modeSpe)

parameter[[2]] <- rep(0.05,9)

masterBlock <- matrix(0.07, ncol = 3, nrow = 3, byrow = TRUE)
diag(masterBlock) <- NA
parameter[[3]] <- q_doubletrans(traits, masterBlock, diff.conceal = FALSE)


run_secsse <- function(num_threads) {
  secsse_loglik_eval(parameter = parameter,
                     phy = phy,
                     traits = traits,
                     num_concealed_states = num_concealed_states,
                     sampling_fraction = sampling_fraction,
                     is_complete_tree = FALSE,
                     num_threads = num_threads,
                     method = "odeint::runge_kutta_fehlberg78",
                     atol = 1e-8,
                     rtol = 1e-6,
                     num_steps = num_steps)
}

rr <- microbenchmark::microbenchmark("single thr." = run_secsse(1),
                                     "2 threads" = run_secsse(2),
                                     "4 threads" = run_secsse(4),
                                     "8 threads" = run_secsse(8),
                                     times = 10)
print(rr)

