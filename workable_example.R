
# generating rates:


# depend_first <- "ETD"
# the_lambdas <- c(0.3,0.1)
# the_q <- 0.1# the transition rates
# 


library(secsse)
library(DDD)
#library(deSolve)
#library(apTreeshape)
#library(foreach)
#library(vegan)
library(stringr)
#library(gtools)
library(Matrix)
source("/Users/thijsjanzen/Downloads/workable_example/fitting_secsse2.R")

sim_dataset <- readRDS("/Users/thijsjanzen/Downloads/workable_example/simulated_data.RDS")
phylotree <- sim_dataset[[1]]
traits <- sim_dataset[[2]]



max(ape::branching.times(phylotree))
table(traits)
# fitting secsse part  
depend_fit_first <- "ETD"

depend_fit_second <- "ETD"


table_rates_ll <- NULL



# crown age is 20.critical_t
scenarios_timecritical <- c(1,4,6,8,10,12,14,16,18,100) # the last number 100 is to try the case with simple secsse_ml 

#for(i in 1:length(scenarios_timecritical)){
  
  
 # critical_t <- scenarios_timecritical[i]
  critical_t <- 10
  # time_stratified_on is whether or not the timezone function is used
  if(critical_t != 100){
    time_stratified_on <- TRUE
  } else {
    time_stratified_on <- FALSE
  }
 
  estimated_rates <- NULL
  estimated_rates <- fit_secsse_now(phylotree,traits,critical_t,depend_fit_first,depend_fit_second, time_stratified_on)

  

  
  table_rates_ll <- rbind(table_rates_ll,
                     c(estimated_rates[1],
                       estimated_rates[2],
                         estimated_rates[3],
                          estimated_rates[4],
                             estimated_rates[5],
                               estimated_rates[6],
                       critical_t
                     ) )
colnames(table_rates_ll) <- c("lambda1_first",
                           "lambda2_first",
                           "q",
                           "lambda1_second",
                           "lambda2_second",
                           "ll",
                           "time_critical"
)



