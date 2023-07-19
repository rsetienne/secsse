library(secsse)
library(RcppParallel)

set.seed(42)
#set.seed(51)
out <- DDD::dd_sim(pars = c(0.5, 0.3, 10000), age = 40)
phy <- out$tes
#plot(phy)
cat("this tree has: ", phy$Nnode + 1, " tips and ", phy$Nnode, " internal nodes\n")


traits <- sample(c(0,1),ape::Ntip(phy),replace = T)
b <- c(0.04,0.04)  # lambda
d <- rep(1,2)
userTransRate <- 0.2 # transition rate among trait states
num_concealed_states <- 2
sampling_fraction <- c(1,1)
toCheck <- secsse::id_paramPos(traits,num_concealed_states)
toCheck[[1]][] <- b
toCheck[[2]][] <- d
toCheck[[3]][,] <- userTransRate
diag(toCheck[[3]]) <- NA
root_state_weight <- "maddison_weights"
use_fortran <- TRUE
methode <- "odeint::bulirsch_stoer"
cond <- "noCondit"

run_secsse <- function(nt) {
  as.numeric(secsse_loglik(parameter = toCheck,
                           phy = phy,
                           traits = traits,
                           num_concealed_states = num_concealed_states,
                           cond = cond,
                           root_state_weight = root_state_weight,
                           sampling_fraction = sampling_fraction,
                           num_threads = nt,
                           is_complete_tree = FALSE))
}

cat("Starting bench...\n")
control <- list("inorder", 2)
names(control) <- c("order", "warmup")
rr <- microbenchmark::microbenchmark(
                                     "1 thread" = run_secsse(1),
                                     "2 threads" = run_secsse(2), 
                                     "4 threads" = run_secsse(4), 
                                     "8 threads" = run_secsse(8), 
                                     control = control,
                                     times = 100)
print(rr)
