context("test_secsse")

test_that("secsse gives the same result as hisse", {
  
  # library(deSolve)
  ## Test to check that our approach reaches the same likelihood than HiSSE.
  # to calculate likelihood of a trait with 2 states using Hisse
  #####pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
  set.seed(4)
  phy <- NULL
  #while (is.null(phy)) {
    ######phy <- diversitree::tree.bisse(pars, max.t=30, x0=0)
  #}
    
    phy <- ape::rcoal(52)
    
  ######traits <- as.numeric(phy$tip.state)
  
  traits<-sample(c(0,1),52,replace = TRUE)
  testit::assert(!is.null(phy))
  testit::assert(length(traits)==52)

  f <- c(1,1)
   
  b <- c(0.04,0.02,0.03,0.04)# lambda
  d <- c(0.03,0.01,0.01,0.02)  # Mu
  userTransRate <- 0.2 # transition rate among trait states
  condition.on.survival <- FALSE
  #x <- as.numeric(secsse::hisse_loglik(phy,traits,b,d,userTransRate,condition.on.survival,f))
  
  condition.on.survival <- TRUE
  #x1 <- as.numeric(secsse::hisse_loglik(phy,traits,b,d,userTransRate,condition.on.survival,f))  

  
  ## Now with different  sampling_fraction 
  f <- c(0.8,1)
  #x2 <- as.numeric(secsse::hisse_loglik(phy,traits,b,d,userTransRate,condition.on.survival,f))  

  ## Now our equivalent version, with only 2 states
  num_concealed_states <- 2

  sampling_fraction <- c(1,1)
  toCheck <- secsse::id_paramPos(traits,num_concealed_states)
  toCheck[[1]][] <- b
  toCheck[[2]][] <- d
  toCheck[[3]][,]<-userTransRate
  diag(toCheck[[3]]) <- NA
  root_state_weight <- "maddison_weights"
  use_fortran<-TRUE
  methode<-"ode45"
  cond<-"noCondit"

  y <- secsse_loglik(parameter=toCheck,phy=phy,
    traits=traits,num_concealed_states=num_concealed_states,
    use_fortran=use_fortran,methode=methode,
    cond=cond,root_state_weight=root_state_weight,sampling_fraction=sampling_fraction)
  
  cond <- "maddison_cond"
  y1 <- round(as.numeric(secsse_loglik(toCheck,phy,traits,num_concealed_states,use_fortran,methode,cond,root_state_weight,sampling_fraction)),4)
  
  ## Now with different  sampling_fraction 
  
  sampling_fraction <- c(0.8,1)
  
  #setwd()<-path.package("secsse")
  y2 <- round(as.numeric(secsse_loglik(toCheck,phy,traits,num_concealed_states,use_fortran,methode,cond,root_state_weight,sampling_fraction)),4)
  
  
  # Test to compare different solvers: Fortran vs Rsolver
  #z1<-round(as.numeric(secsse_loglik(parameter=toCheck,phy,traits,num_concealed_states,use_fortran="Rsolver",methode="ode45",cond,sampling_fraction)),4)
  #dyn.unload(directoryplusfilename)
  
  #setwd()<-path.package("secsse")
  #z2<-round(as.numeric(secsse_loglik(parameter=toCheck,phy,traits,num_concealed_states,use_fortran="Fortran",methode="ode45",cond,sampling_fraction)),4)
   
  testthat::expect_equal(-237.8611,y1)##-237.8611 is the right one, 
  testthat::expect_equal(-243.8611,y2)
  #testthat::expect_equal(z1, z2) 
  
 
})


