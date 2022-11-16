context("test_timezones")

test_that("secsse gives the same result with one or two identical timezones", {
  
  newickphy <- "((((t15:0.03654175604,t36:0.03654175604):0.1703092581,(((t42:0.01312768801,t23:0.01312768801):0.01026551964,(((t19:0.006565648042,t5:0.006565648042):0.000589637007,t35:0.007155285049):0.0075478055,t51:0.01470309055):0.008690117099):0.1040593382,(t20:0.05092066659,t16:0.05092066659):0.07653187925):0.07939846827):0.6519637868,(((((t43:0.006616860045,t3:0.006616860045):0.08611719299,(t48:0.004896235936,t40:0.004896235936):0.0878378171):0.1515206506,((t44:0.09487672192,t2:0.09487672192):0.07712689077,((t37:0.006132013467,t32:0.006132013467):0.1177191576,((t46:0.01830302153,t21:0.01830302153):0.03858278382,((t25:0.02071187578,t24:0.02071187578):0.02799215338,t47:0.04870402916):0.008181776188):0.06696536571):0.04815244163):0.07225109099):0.03049659492,((t6:0.02021971253,t45:0.02021971253):0.1267950773,t18:0.1470147899):0.1277365087):0.5391698492,(((((t27:0.008082361089,t17:0.008082361089):0.00456225043,t39:0.01264461152):0.103375347,(t7:0.06545659749,((t26:0.005452218586,t12:0.005452218586):0.03594003265,((t13:0.0001294122247,t9:0.0001294122247):0.01141726784,t31:0.01154668006):0.02984557118):0.02406434625):0.05056336106):0.04543362477,((t34:0.0748070545,t11:0.0748070545):0.01677840675,(((t38:0.01479762241,(t41:0.004213712966,t14:0.004213712966):0.01058390944):0.000225587269,t4:0.01502320968):0.06205778867,((t49:0.01206564111,(t10:0.00350505531,t52:0.00350505531):0.008560585803):0.03485629493,(t28:0.04155870788,((t8:0.01119536676,t22:0.01119536676):0.02493294048,t50:0.03612830725):0.005430400635):0.005363228164):0.0301590623):0.01450446291):0.06986812207):0.1092343488,(t1:0.1156934975,t30:0.1156934975):0.1549944346):0.5432332157):0.04489365312):1.400701854,(t29:0.04276331213,t33:0.04276331213):2.216753343);"
  phy <- phytools::read.newick(text = newickphy)
  testit::assert(!is.null(phy))
  traits <- c(0,1,1,0,1,0,0,0,1,1,0,1,0,1,1,0,1,0,0,0,0,1,0,1,1,1,1,1,0,1,0,1,1,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0)
  
  #f <- c(1,1)
  b <- c(0.04,0.02,0.03,0.04)# lambda
  d <- c(0.03,0.01,0.01,0.02)  # Mu
  userTransRate <- 0.2 # transition rate among trait states
  
  ## Now our equivalent version, with only 2 states
  num_concealed_states <- 2
  sampling_fraction <- c(1,1)
  toCheck <- id_paramPos(traits,num_concealed_states)
  toCheck[[1]][] <- b
  toCheck[[2]][] <- d
  toCheck[[3]][,] <- userTransRate
  diag(toCheck[[3]]) <- NA
  root_state_weight <- "maddison_weights"
  cond <- "noCondit"
  
  y <- secsse::secsse_loglik(parameter = toCheck,
                             phy = phy,
                             traits = traits,
                             num_concealed_states = num_concealed_states,
                             cond = cond,
                             root_state_weight = root_state_weight,
                             sampling_fraction = sampling_fraction)
  
  toCheck2 <- list()
  toCheck2[[1]] <- toCheck$lambdas
  toCheck2[[2]] <- toCheck$mus
  toCheck2[[3]] <- toCheck$Q
  toCheck2[[4]] <- toCheck2[[1]]
  toCheck2[[5]] <- toCheck2[[2]]
  toCheck2[[6]] <- toCheck2[[3]]
  
  y2 <- secsse::secsse_loglik_timezones(parameter = toCheck2,
                                        phy = phy,
                                        traits = traits,
                                        critical_t = 0.5,
                                        num_concealed_states = num_concealed_states,
                                        cond = cond,
                                        root_state_weight = root_state_weight,
                                        sampling_fraction = sampling_fraction)
  y2
  testthat::expect_equal(y, y2)
  
  
  for (crit_t in ape::branching.times(phy)[1:5]) {
    y2 <- secsse::secsse_loglik_timezones(parameter = toCheck2,
                                          phy = phy,
                                          traits = traits,
                                          critical_t = crit_t,
                                          num_concealed_states = num_concealed_states,
                                          cond = cond,
                                          root_state_weight = root_state_weight,
                                          sampling_fraction = sampling_fraction)
    
    testthat::expect_equal(y, y2)
  }
  
  toCheck2[[5]] <- c(0.0, 0.0, 0.0, 0.0)
  found <- c()
  for (crit_t in ape::branching.times(phy)[2:10]) {
    y2 <- secsse::secsse_loglik_timezones(parameter = toCheck2,
                                          phy = phy,
                                          traits = traits,
                                          critical_t = crit_t,
                                          num_concealed_states = num_concealed_states,
                                          cond = cond,
                                          root_state_weight = root_state_weight,
                                          sampling_fraction = sampling_fraction)
    cat(crit_t, y - y2, y, y2, "\n")
    found <- rbind(found, c(crit_t, y - y2))
    testthat::expect_false(y == y2)
  }
  A <- lm(found[,2]~log(found[,1]))  
  ax <- summary(A)  
  ay <- ax$coefficients
  testthat::expect_true(ay[1, 2] > 0)
})

test_that("secsse gives the same result with one or two identical timezones", {
  
  newickphy <- "((((t15:0.03654175604,t36:0.03654175604):0.1703092581,(((t42:0.01312768801,t23:0.01312768801):0.01026551964,(((t19:0.006565648042,t5:0.006565648042):0.000589637007,t35:0.007155285049):0.0075478055,t51:0.01470309055):0.008690117099):0.1040593382,(t20:0.05092066659,t16:0.05092066659):0.07653187925):0.07939846827):0.6519637868,(((((t43:0.006616860045,t3:0.006616860045):0.08611719299,(t48:0.004896235936,t40:0.004896235936):0.0878378171):0.1515206506,((t44:0.09487672192,t2:0.09487672192):0.07712689077,((t37:0.006132013467,t32:0.006132013467):0.1177191576,((t46:0.01830302153,t21:0.01830302153):0.03858278382,((t25:0.02071187578,t24:0.02071187578):0.02799215338,t47:0.04870402916):0.008181776188):0.06696536571):0.04815244163):0.07225109099):0.03049659492,((t6:0.02021971253,t45:0.02021971253):0.1267950773,t18:0.1470147899):0.1277365087):0.5391698492,(((((t27:0.008082361089,t17:0.008082361089):0.00456225043,t39:0.01264461152):0.103375347,(t7:0.06545659749,((t26:0.005452218586,t12:0.005452218586):0.03594003265,((t13:0.0001294122247,t9:0.0001294122247):0.01141726784,t31:0.01154668006):0.02984557118):0.02406434625):0.05056336106):0.04543362477,((t34:0.0748070545,t11:0.0748070545):0.01677840675,(((t38:0.01479762241,(t41:0.004213712966,t14:0.004213712966):0.01058390944):0.000225587269,t4:0.01502320968):0.06205778867,((t49:0.01206564111,(t10:0.00350505531,t52:0.00350505531):0.008560585803):0.03485629493,(t28:0.04155870788,((t8:0.01119536676,t22:0.01119536676):0.02493294048,t50:0.03612830725):0.005430400635):0.005363228164):0.0301590623):0.01450446291):0.06986812207):0.1092343488,(t1:0.1156934975,t30:0.1156934975):0.1549944346):0.5432332157):0.04489365312):1.400701854,(t29:0.04276331213,t33:0.04276331213):2.216753343);"
  phy <- phytools::read.newick(text = newickphy)
  testit::assert(!is.null(phy))
  traits <- c(0,1,1,0,1,0,0,0,1,1,0,1,0,1,1,0,1,0,0,0,0,1,0,1,1,1,1,1,0,1,0,1,1,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0)
  
  #f <- c(1,1)
  b <- c(0.04,0.02,0.03,0.04)# lambda
  d <- c(0.03,0.01,0.01,0.02)  # Mu
  userTransRate <- 0.2 # transition rate among trait states
  
  ## Now our equivalent version, with only 2 states
  num_concealed_states <- 2
  sampling_fraction <- c(1,1)
  toCheck <- id_paramPos(traits,num_concealed_states)
  toCheck[[1]][] <- b
  toCheck[[2]][] <- d
  toCheck[[3]][,] <- userTransRate
  diag(toCheck[[3]]) <- NA
  root_state_weight <- "maddison_weights"
  cond <- "noCondit"
  
  y <- secsse::secsse_loglik(parameter = toCheck,
                             phy = phy,
                             traits = traits,
                             num_concealed_states = num_concealed_states,
                             cond = cond,
                             root_state_weight = root_state_weight,
                             sampling_fraction = sampling_fraction)
  
  toCheck2 <- list()
  toCheck2[[1]] <- toCheck$lambdas
  toCheck2[[2]] <- toCheck$mus
  toCheck2[[3]] <- toCheck$Q
  toCheck2[[4]] <- toCheck2[[1]]
  toCheck2[[5]] <- toCheck2[[2]]
  toCheck2[[6]] <- toCheck2[[3]]
  
  y2 <- secsse::secsse_loglik_timezones(parameter = toCheck2,
                                        phy = phy,
                                        traits = traits,
                                        critical_t = 0.5,
                                        num_concealed_states = num_concealed_states,
                                        cond = cond,
                                        root_state_weight = root_state_weight,
                                        sampling_fraction = sampling_fraction)
  y2
  testthat::expect_equal(y, y2)
  
  toCheck2[[2]] <- c(0.0, 0.0, 0.0, 0.0)
  y2_2 <- secsse::secsse_loglik_timezones(parameter = toCheck2,
                                        phy = phy,
                                        traits = traits,
                                        critical_t = 0.5,
                                        num_concealed_states = num_concealed_states,
                                        cond = cond,
                                        root_state_weight = root_state_weight,
                                        sampling_fraction = sampling_fraction)
  
  testthat::expect_false(y2 == y2_2)
  
  toCheck3 <- list()
  toCheck3[[1]] <- toCheck
  toCheck3[[2]] <- toCheck
 
  y3 <- secsse::secsse_loglik_timezones_compound(parameter = toCheck3,
                                                 phy = phy,
                                                 traits = traits,
                                                 critical_t = c(0.5),
                                                 num_concealed_states = num_concealed_states,
                                                 cond = cond,
                                                 root_state_weight = root_state_weight,
                                                 sampling_fraction = sampling_fraction,
                                                 num_threads = 1);
  
  testthat::expect_equal(y2, y3)
  
  toCheck3 <- list()
  toCheck3[[1]] <- toCheck
  toCheck3[[2]] <- toCheck
  toCheck3[[1]]$mus <- c(0.0, 0.0, 0.0, 0.0)
  
  y3 <- secsse::secsse_loglik_timezones_compound(parameter = toCheck3,
                                                 phy = phy,
                                                 traits = traits,
                                                 critical_t = c(0.5),
                                                 num_concealed_states = num_concealed_states,
                                                 cond = cond,
                                                 root_state_weight = root_state_weight,
                                                 sampling_fraction = sampling_fraction,
                                                 num_threads = 1);
  
  testthat::expect_false(y2 == y3)
  
  
  
  
  toCheck3 <- list()
  toCheck3[[1]] <- toCheck
  toCheck3[[2]] <- toCheck
  toCheck3[[3]] <- toCheck
  
  y3 <- secsse::secsse_loglik_timezones_compound(parameter = toCheck3,
                                                 phy = phy,
                                                 traits = traits,
                                                 critical_t = c(0.5, 1.0),
                                                 num_concealed_states = num_concealed_states,
                                                 cond = cond,
                                                 root_state_weight = root_state_weight,
                                                 sampling_fraction = sampling_fraction,
                                                 num_threads = 1);
  
  testthat::expect_equal(y2, y3)
  
  
  
  toCheck3 <- list()
  toCheck3[[1]] <- toCheck
  toCheck3[[2]] <- toCheck
  toCheck3[[1]]$mus <- c(0.0, 0.0, 0.0, 0.0)
  toCheck3[[3]] <- toCheck
  
  y3 <- secsse::secsse_loglik_timezones_compound(parameter = toCheck3,
                                                 phy = phy,
                                                 traits = traits,
                                                 critical_t = c(0.5, 1.0),
                                                 num_concealed_states = num_concealed_states,
                                                 cond = cond,
                                                 root_state_weight = root_state_weight,
                                                 sampling_fraction = sampling_fraction,
                                                 num_threads = 1);
  
  testthat::expect_false(y2 == y3)
})