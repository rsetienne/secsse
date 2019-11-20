secsse_test_hisse <- function(){
  ## Test to check that our approach reaches the same likelihood than HiSSE.
  # to calculate likelihood of a trait with 2 states using Hisse
  #####pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
  #set.seed(4); phy <- ape::rcoal(52)
  newickphy <- "((((t15:0.03654175604,t36:0.03654175604):0.1703092581,(((t42:0.01312768801,t23:0.01312768801):0.01026551964,(((t19:0.006565648042,t5:0.006565648042):0.000589637007,t35:0.007155285049):0.0075478055,t51:0.01470309055):0.008690117099):0.1040593382,(t20:0.05092066659,t16:0.05092066659):0.07653187925):0.07939846827):0.6519637868,(((((t43:0.006616860045,t3:0.006616860045):0.08611719299,(t48:0.004896235936,t40:0.004896235936):0.0878378171):0.1515206506,((t44:0.09487672192,t2:0.09487672192):0.07712689077,((t37:0.006132013467,t32:0.006132013467):0.1177191576,((t46:0.01830302153,t21:0.01830302153):0.03858278382,((t25:0.02071187578,t24:0.02071187578):0.02799215338,t47:0.04870402916):0.008181776188):0.06696536571):0.04815244163):0.07225109099):0.03049659492,((t6:0.02021971253,t45:0.02021971253):0.1267950773,t18:0.1470147899):0.1277365087):0.5391698492,(((((t27:0.008082361089,t17:0.008082361089):0.00456225043,t39:0.01264461152):0.103375347,(t7:0.06545659749,((t26:0.005452218586,t12:0.005452218586):0.03594003265,((t13:0.0001294122247,t9:0.0001294122247):0.01141726784,t31:0.01154668006):0.02984557118):0.02406434625):0.05056336106):0.04543362477,((t34:0.0748070545,t11:0.0748070545):0.01677840675,(((t38:0.01479762241,(t41:0.004213712966,t14:0.004213712966):0.01058390944):0.000225587269,t4:0.01502320968):0.06205778867,((t49:0.01206564111,(t10:0.00350505531,t52:0.00350505531):0.008560585803):0.03485629493,(t28:0.04155870788,((t8:0.01119536676,t22:0.01119536676):0.02493294048,t50:0.03612830725):0.005430400635):0.005363228164):0.0301590623):0.01450446291):0.06986812207):0.1092343488,(t1:0.1156934975,t30:0.1156934975):0.1549944346):0.5432332157):0.04489365312):1.400701854,(t29:0.04276331213,t33:0.04276331213):2.216753343);"
  phy <- phytools::read.newick(text = newickphy)
  testit::assert(!is.null(phy))
  traits <- c(0,1,1,0,1,0,0,0,1,1,0,1,0,1,1,0,1,0,0,0,0,1,0,1,1,1,1,1,0,1,0,1,1,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0)
  
  #f <- c(1,1)
  b <- c(0.04,0.02,0.03,0.04)# lambda
  d <- c(0.03,0.01,0.01,0.02)  # Mu
  userTransRate <- 0.2 # transition rate among trait states
  #condition.on.survival <- FALSE
  #x <- as.numeric(secsse::hisse_loglik(phy,traits,b,d,userTransRate,condition.on.survival,f))
  #condition.on.survival <- TRUE
  #x1 <- as.numeric(secsse::hisse_loglik(phy,traits,b,d,userTransRate,condition.on.survival,f))  
  ## Now with different  sampling_fraction 
  #f <- c(0.8,1)
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
  
  y <- secsse_loglik(parameter = toCheck,
                     phy = phy,
                     traits = traits,
                     num_concealed_states = num_concealed_states,
                     use_fortran = use_fortran,
                     methode = methode,
                     cond = cond,
                     root_state_weight = root_state_weight,
                     sampling_fraction = sampling_fraction)
  
  cond <- "maddison_cond"
  y1 <- round(as.numeric(secsse_loglik(parameter = toCheck,
                                      phy = phy,
                                      traits = traits,
                                      num_concealed_states = num_concealed_states,
                                      use_fortran = TRUE,
                                      methode = "ode45",
                                      cond = cond,
                                      root_state_weight = root_state_weight,
                                      sampling_fraction = sampling_fraction)
  ),4)
  
  ## Now with different  sampling_fraction 
  
  sampling_fraction <- c(0.8,1)
  
  y2 <- round(as.numeric(secsse_loglik(parameter = toCheck,
                                       phy = phy,
                                       traits = traits,
                                       num_concealed_states = num_concealed_states,
                                       use_fortran = TRUE,
                                       methode = "ode45",
                                       cond = cond,
                                       root_state_weight = root_state_weight,
                                       sampling_fraction = sampling_fraction)
  ),4)
  
  # Test to compare different solvers: Fortran vs Rsolver
  z1 <- round(as.numeric(secsse_loglik(parameter = toCheck,
                                       phy = phy,
                                       traits = traits,
                                       num_concealed_states = num_concealed_states,
                                       use_fortran = FALSE,
                                       methode = "ode45",
                                       cond = cond,
                                       root_state_weight = root_state_weight,
                                       sampling_fraction = sampling_fraction)
                         ),4)
  z2 <- round(as.numeric(secsse_loglik(parameter = toCheck,
                                       phy = phy,
                                       traits = traits,
                                       num_concealed_states = num_concealed_states,
                                       use_fortran = TRUE,
                                       methode = "ode45",
                                       cond = cond,
                                       root_state_weight = root_state_weight,
                                       sampling_fraction = sampling_fraction)
                         ),4)
  testthat::expect_equal(-237.8611,y1)##-237.8611 is the right one, 
  testthat::expect_equal(-243.8611,y2)
  testthat::expect_equal(z1, z2) 
}

secsse_test_geosse <- function(){
  #geosse
  pars <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
  names(pars) <- c("sA", "sB", "sAB", "xA", "xB", "dA", "dB")
  #set.seed(5)
  #phy <- diversitree::tree.geosse(pars, max.t=4, x0=0)
  phy <- NULL; rm(phy);
  utils::data('example_phy_GeoSSE', package = 'secsse');
  traits <- as.numeric(phy$tip.state)
  testit::assert(!is.null(phy))
  lik.g <- diversitree::make.geosse(phy, phy$tip.state)
  pars.g <- c(1.5, 0.5, 1.0, 0.7, 0.7, 1.4, 1.3)
  names(pars.g) <- diversitree::argnames(lik.g)
  lik.c <- diversitree::make.classe(phy, phy$tip.state+1, 3)
  pars.c <- 0 * diversitree::starting.point.classe(phy, 3)
  pars.c['lambda222'] <- pars.c['lambda112'] <- pars.g['sA']
  pars.c['lambda333'] <- pars.c['lambda113'] <- pars.g['sB']
  pars.c['lambda123'] <- pars.g['sAB']
  pars.c['mu2'] <- pars.c['q13'] <- pars.g['xA']
  pars.c['mu3'] <- pars.c['q12'] <- pars.g['xB']
  pars.c['q21'] <- pars.g['dA']
  pars.c['q31'] <- pars.g['dB']
  lik.g(pars.g) # -175.7685
  classe_diversitree_LL <- lik.c(pars.c) # -175.7685
  
  ## Secsse part 
  lambdas<-list()
  lambdas[[1]] <- matrix(0,ncol = 3,nrow = 3,byrow = TRUE)
  #lambdas[[1]][1,1] <- 1.5
  lambdas[[1]][2,1] <- 1.5
  lambdas[[1]][3,1] <- 0.5
  lambdas[[1]][3,2] <- 1
  lambdas[[2]]<-matrix(0,ncol = 3,nrow = 3,byrow = TRUE)
  lambdas[[2]][2,2] <- 1.5
  #lambdas[[2]][2,2] <- 1.1
  lambdas[[3]] <- matrix(0,ncol=3,nrow=3,byrow=TRUE)
  lambdas[[3]][3,3] <- 0.5
  #lambdas[[3]][3,3] <- 1
  
  mus<-c(0,0.7,0.7)
  
  q<-matrix(0,ncol = 3,nrow = 3,byrow = TRUE)
  q[2,1] <- 1.4
  q[3,1] <- 1.3
  q[1,2] <- 0.7
  q[1,3] <- 0.7
  
  parameter <- list()
  parameter[[1]] <- lambdas
  parameter[[2]] <- mus
  parameter[[3]] <- q
  
  num_concealed_states <- 3.1
  
  secsse_cla_LL <- cla_secsse_loglik(parameter,
                                     phy,
                                     traits,
                                     num_concealed_states,
                                     use_fortran = FALSE,
                                     methode = "ode45",
                                     cond = "maddison_cond",
                                     root_state_weight = "maddison_weights",
                                     sampling_fraction = c(1,1,1),
                                     run_parallel = FALSE,
                                     setting_calculation = NULL,
                                     setting_parallel = NULL,
                                     see_ancestral_states = FALSE,
                                     loglik_penalty = 0)
  
  testthat::expect_equal(classe_diversitree_LL,secsse_cla_LL)
  
  secsse_cla_LL2 <- cla_secsse_loglik(parameter,
                                      phy,
                                      traits,
                                      num_concealed_states,
                                      use_fortran = TRUE,
                                      methode = "ode45",
                                      cond = "maddison_cond",
                                      root_state_weight = "maddison_weights",
                                      sampling_fraction = c(1,1,1),
                                      run_parallel = FALSE,
                                      setting_calculation = NULL,
                                      setting_parallel = NULL,
                                      see_ancestral_states = FALSE,
                                      loglik_penalty = 0)
  #print(secsse_cla_LL)
  #print(secsse_cla_LL2)
  testthat::expect_equal(secsse_cla_LL,secsse_cla_LL2)
}

secsse_test_ml <- function(){
  parenthesis<-"(((6:0.2547423371,(1:0.0496153503,4:0.0496153503):0.2051269868):0.1306304758,(9:0.2124135406,5:0.2124135406):0.1729592723):1.151205247,(((7:0.009347664296,3:0.009347664296):0.2101416075,10:0.2194892718):0.1035186448,(2:0.2575886319,8:0.2575886319):0.06541928469):1.213570144);"
  phylotree<-ape::read.tree(file="",parenthesis)
  traits <-  c(2, 0, 1, 0, 2, 0, 1, 2, 2, 0)
  num_concealed_states <- 3
  idparslist <- id_paramPos(traits, num_concealed_states)
  idparslist[[1]][c(1,4,7)] <- 1
  idparslist[[1]][c(2,5,8)] <- 2
  idparslist[[1]][c(3,6,9)] <- 3
  idparslist[[2]][] <- 4
  masterBlock<-matrix(5,ncol = 3,nrow = 3,byrow = TRUE) 
  diag(masterBlock) <- NA
  diff.conceal <- FALSE
  idparslist[[3]] <- q_doubletrans(traits,masterBlock,diff.conceal)
  startingpoint<-DDD::bd_ML(brts = ape::branching.times(phylotree))
  intGuessLamba <- startingpoint$lambda0
  intGuessMu <- startingpoint$mu0
  idparsopt <- c(1,2,3)
  initparsopt <- c(rep(intGuessLamba,3))
  idparsfix <- c(0,4,5)
  parsfix <- c(0,0,0.1)
  tol = c(1e-04, 1e-05, 1e-07)
  maxiter = 1000 * round((1.25)^length(idparsopt))
  use_fortran = TRUE
  methode = "ode45"
  optimmethod = "simplex"
  run_parallel = TRUE
  cond<-"proper_cond"
  root_state_weight <- "proper_weights"
  sampling_fraction <- c(1,1,1)
  model<-secsse_ml(
    phylotree,
    traits,
    num_concealed_states,
    idparslist,
    idparsopt,
    initparsopt,
    idparsfix,
    parsfix,
    cond,
    root_state_weight,
    sampling_fraction,
    tol,
    maxiter,
    use_fortran,
    methode,
    optimmethod,
    num_cycles = 1,
    run_parallel
  )
  
  testthat::expect_equal(model$ML,-29.89993)
}
  
secsse_test_ml2 <- function(){
  parenthesis <- "(((6:0.2547423371,(1:0.0496153503,4:0.0496153503):0.2051269868):0.1306304758,(9:0.2124135406,5:0.2124135406):0.1729592723):1.151205247,(((7:0.009347664296,3:0.009347664296):0.2101416075,10:0.2194892718):0.1035186448,(2:0.2575886319,8:0.2575886319):0.06541928469):1.213570144);"
  phylotree <- ape::read.tree(file="",parenthesis)
  traits <- c(2, 0, 1, 0, 2, 0, 1, 2, 2, 0)
  startingpoint <- DDD::bd_ML(brts = ape::branching.times(phylotree))
  intGuessLamba <- startingpoint$lambda0
  intGuessMu <- startingpoint$mu0
  num_concealed_states <- 3
  idparslist <- id_paramPos(traits, num_concealed_states)
  idparslist[[1]][] <- 1
  idparslist[[2]][] <- 2
  masterBlock <- matrix(c(3,4,3,4,3,4,3,4,3),ncol = 3,nrow = 3,byrow = TRUE)
  diag(masterBlock) <- NA
  diff.conceal <- FALSE
  idparslist[[3]] <- q_doubletrans(traits,masterBlock,diff.conceal)
  idparsfuncdefpar <- c(3)
  idparsopt <- c(1,4)
  idparsfix <- c(0,2)
  initparsopt <- c(rep(intGuessLamba,1),intGuessLamba/5)
  parsfix <- c(0,5)
  idfactorsopt <- 1
  initfactors <- 1
  functions_defining_params <- list()
  par_4 <- NA 
  factor_1 <- NA
  functions_defining_params[[1]] <- function(){
    par_3 <- par_4 * factor_1
  }
  tol <- c(1e-04, 1e-05, 1e-07)
  maxiter <- 1000 * round((1.25)^length(idparsopt))
  use_fortran <- TRUE
  methode <- "ode45"
  optimmethod <- "simplex"
  run_parallel <- FALSE
  cond <- "proper_cond"
  root_state_weight <- "proper_weights"
  sampling_fraction <- c(1,1,1)
  model <- secsse_ml_func_def_pars(phy = phylotree,
                                   traits = traits,
                                   num_concealed_states = num_concealed_states,
                                   idparslist = idparslist,
                                   idparsopt = idparsopt,
                                   initparsopt = initparsopt,
                                   idfactorsopt = idfactorsopt,
                                   initfactors = initfactors,
                                   idparsfix = idparsfix,
                                   parsfix = parsfix,
                                   idparsfuncdefpar = idparsfuncdefpar,
                                   functions_defining_params = functions_defining_params,
                                   cond = cond,
                                   root_state_weight = root_state_weight,
                                   sampling_fraction = sampling_fraction,
                                   tol = tol,
                                   maxiter = maxiter,
                                   use_fortran = use_fortran,
                                   methode = methode,
                                   optimmethod = optimmethod,
                                   num_cycles = 1,
                                   run_parallel = run_parallel)
  testthat::expect_equal(model$ML,-12.87974)
}  

secsse_test_ml3 <- function(){
  parenthesis <- "(((6:0.2547423371,(1:0.0496153503,4:0.0496153503):0.2051269868):0.1306304758,(9:0.2124135406,5:0.2124135406):0.1729592723):1.151205247,(((7:0.009347664296,3:0.009347664296):0.2101416075,10:0.2194892718):0.1035186448,(2:0.2575886319,8:0.2575886319):0.06541928469):1.213570144);"
  phylotree <- ape::read.tree(file = "",parenthesis)
  traits <- c(2, 0, 1, 0, 2, 0, 1, 2, 2, 0)
  num_concealed_states <- 3
  idparslist <- cla_id_paramPos(traits,num_concealed_states)
  idparslist$lambdas[2,] <- rep(1,9)
  idparslist[[2]][] <- 4
  masterBlock <- matrix(5,ncol = 3,nrow = 3,byrow = TRUE) 
  diag(masterBlock) <- NA
  diff.conceal <- FALSE
  idparslist[[3]] <- q_doubletrans(traits,masterBlock,diff.conceal)
  startingpoint <- DDD::bd_ML(brts = ape::branching.times(phylotree))
  intGuessLamba <- startingpoint$lambda0
  intGuessMu <- startingpoint$mu0
  idparsopt <- c(1)
  initparsopt <- c(rep(intGuessLamba,1))
  idparsfix <- c(0,4,5)
  parsfix <- c(0,0,0.01)
  tol <- c(1e-04, 1e-05, 1e-07)
  maxiter <- 1000 * round((1.25)^length(idparsopt))
  use_fortran <- FALSE
  methode <- "ode45"
  optimmethod <- "simplex"
  run_parallel <- FALSE
  cond <- "proper_cond"
  root_state_weight <- "proper_weights"
  sampling_fraction <- c(1,1,1)
  
  model<-cla_secsse_ml(
    phylotree,
    traits,
    num_concealed_states,
    idparslist,
    idparsopt,
    initparsopt,
    idparsfix,
    parsfix,
    cond,
    root_state_weight,
    sampling_fraction,
    tol,
    maxiter,
    use_fortran,
    methode,
    optimmethod,
    num_cycles = 1,
    run_parallel)
  testthat::expect_equal(model$ML,-16.1342246206186)
}  
