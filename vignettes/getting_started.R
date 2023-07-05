## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----data---------------------------------------------------------------------
focal_tree <- ape::read.tree(text = "((t1:2.586692974,((t13:0.8873239777,(t34:0.7528802927,t44:0.7528802927):0.134443685):1.04055583,t21:1.927879808):0.6588131659):2.413307026,(((((t2:0.5479660671,t49:0.5479660671):1.365717301,((t22:0.6512648145,(t47:0.2049372686,t61:0.2049372686):0.446327546):0.2120207488,((t36:0.1639000795,t62:0.1639000795):0.6477110608,(t40:0.7930189281,t41:0.7930189281):0.01859221218):0.05167442306):1.050397805):0.3472868921,t17:2.26097026):0.06999807465,((t15:0.3772986189,t55:0.3772986189):1.943391575,(t16:1.196137061,((t26:0.9825148624,((t31:0.02350286365,t67:0.02350286365):0.9171104913,(t32:0.7837716587,t42:0.7837716587):0.1568416962):0.04190150739):0.10830946,t28:1.090824322):0.1053127385):1.124553133):0.01027814129):2.38725674,((((t3:1.083671921,(((t29:0.01474991842,t68:0.01474991842):0.07607893188,t66:0.0908288503):0.1671480643,t59:0.2579769146):0.8256950068):1.890217226,((((t8:0.5490733351,(t48:0.297289906,t58:0.297289906):0.2517834291):0.513828506,t30:1.062901841):1.550593154,(((t12:0.7280698912,(((t45:0.2153519103,t60:0.2153519103):0.1016221231,t57:0.3169740334):0.2050447992,t51:0.5220188326):0.2060510586):0.6068040975,(t23:0.09157856985,t65:0.09157856985):1.243295419):0.8042199558,(t18:1.142634904,((t27:0.4159759965,(t54:0.1286255831,t63:0.1286255831):0.2873504134):0.3623222946,t43:0.7782982911):0.3643366127):0.9964590407):0.4744010501):0.01756278112,(t11:1.233211103,(t25:0.09389678922,t64:0.09389678922):1.139314314):1.397846673):0.3428313717):1.179614451,(((t5:1.939201203,((t20:0.510099777,t53:0.510099777):0.008496491915,t52:0.5185962689):1.420604934):0.8473155895,(t9:2.345942038,((t14:0.8336233238,t38:0.8336233238):0.09720732437,t33:0.9308306482):1.41511139):0.4405747545):0.9289638311,((t6:2.048542159,(((t19:0.5446778552,t50:0.5446778552):0.3019110389,(t37:0.8329466556,t39:0.8329466556):0.01364223839):0.02536522407,(t35:0.7210670433,(t46:0.3394466403,t56:0.3394466403):0.381620403):0.1508870748):1.176588041):1.260277226,(t7:2.68471232,(t10:1.267107488,t24:1.267107488):1.417604832):0.6241070654):0.4066612378):0.4380229751):0.1283519611,t4:4.28185556):0.4363695154):0.2817749251):0;")
focal_traits <- data.frame(trait = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                     1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                                     1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0),
  row.names = focal_tree$tip.label)
if (requireNamespace("diversitree")) {
diversitree::trait.plot(focal_tree, dat = focal_traits, 
                        cols = list("trait" = c("blue", "red")),
                        type = "p")
}

## ----ETD_lambda---------------------------------------------------------------
spec_matrix <- c()
spec_matrix <- rbind(spec_matrix, c(0, 0, 0, 1))
spec_matrix <- rbind(spec_matrix, c(1, 1, 1, 2))
lambda_list <- secsse::create_lambda_list(state_names = c(0, 1),
                                          num_concealed_states = 2,
                                          transition_matrix = spec_matrix,
                                          model = "ETD")
lambda_list

## ----ETD_mu-------------------------------------------------------------------
mu_vec <- secsse::create_mu_vector(state_names = c(0, 1),
                                   num_concealed_states = 2,
                                   model = "ETD",
                                   lambda_list = lambda_list)
mu_vec

## ----ETD_Q--------------------------------------------------------------------
shift_matrix <- c()
shift_matrix <- rbind(shift_matrix, c(0, 1, 5))
shift_matrix <- rbind(shift_matrix, c(1, 0, 6))

q_matrix <- secsse::create_q_matrix(state_names = c(0, 1),
                                    num_concealed_states = 2,
                                    shift_matrix = shift_matrix,
                                    diff.conceal = TRUE)
q_matrix

## ----ETD_ML_init--------------------------------------------------------------
idparsopt <- 1:8 # our maximum rate parameter was 8
idparsfix <- c(0) # we want to keep al zeros at zero
initparsopt <- rep(0.1, 8) 
initparsfix <- c(0.0) # all zeros remain at zero.
sampling_fraction <- c(1, 1)

## ----ETD_ML-------------------------------------------------------------------

idparslist <- list()
idparslist[[1]] <- lambda_list
idparslist[[2]] <- mu_vec
idparslist[[3]] <- q_matrix

answ <- secsse::cla_secsse_ml(phy = focal_tree,
                              traits = focal_traits$trait,
                              num_concealed_states = 2,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = sampling_fraction,
                              verbose = FALSE,
                              num_threads = 4)

## ----ETD_res------------------------------------------------------------------
ML_ETD <- answ$ML
ETD_par <- secsse::extract_par_vals(idparslist, answ$MLpars)
ML_ETD
ETD_par
spec_rates <- ETD_par[1:2]
ext_rates <- ETD_par[3:4]
Q_Examined <- ETD_par[5:6]
Q_Concealed <- ETD_par[7:8]
spec_rates
ext_rates
Q_Examined
Q_Concealed

## ----CTD_lambda---------------------------------------------------------------
spec_matrix <- c()
spec_matrix <- rbind(spec_matrix, c(0, 0, 0, 1))
spec_matrix <- rbind(spec_matrix, c(1, 1, 1, 2))
lambda_list <- secsse::create_lambda_list(state_names = c(0, 1),
                                          num_concealed_states = 2,
                                          transition_matrix = spec_matrix,
                                          model = "CTD")
lambda_list

## ----CTD_mu-------------------------------------------------------------------
mu_vec <- secsse::create_mu_vector(state_names = c(0, 1),
                                   num_concealed_states = 2,
                                   model = "CTD",
                                   lambda_list = lambda_list)
mu_vec

## ----CTD_Q--------------------------------------------------------------------
shift_matrix <- c()
shift_matrix <- rbind(shift_matrix, c(0, 1, 5))
shift_matrix <- rbind(shift_matrix, c(1, 0, 6))

q_matrix <- secsse::create_q_matrix(state_names = c(0, 1),
                                    num_concealed_states = 2,
                                    shift_matrix = shift_matrix,
                                    diff.conceal = TRUE)
q_matrix

## ----CTD_ML-------------------------------------------------------------------
idparsopt <- 1:8 # our maximum rate parameter was 8
idparsfix <- c(0) # we want to keep al zeros at zero
initparsopt <- rep(0.1, 8)
initparsfix <- c(0.0) # all zeros remain at zero.
sampling_fraction <- c(1, 1)

idparslist <- list()
idparslist[[1]] <- lambda_list
idparslist[[2]] <- mu_vec
idparslist[[3]] <- q_matrix

answ <- secsse::cla_secsse_ml(phy = focal_tree,
                              traits = focal_traits$trait,
                              num_concealed_states = 2,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = sampling_fraction,
                              verbose = FALSE,
                              num_threads = 4)
ML_CTD <- answ$ML
CTD_par <- secsse::extract_par_vals(idparslist, answ$MLpars)
ML_CTD
CTD_par
spec_rates <- CTD_par[1:2]
ext_rates <- CTD_par[3:4]
Q_Examined <- CTD_par[5:6]
Q_Concealed <- CTD_par[7:8]
spec_rates
ext_rates
Q_Examined
Q_Concealed

## ----CR_lambda----------------------------------------------------------------
spec_matrix <- c()
spec_matrix <- rbind(spec_matrix, c(0, 0, 0, 1))
spec_matrix <- rbind(spec_matrix, c(1, 1, 1, 1))
lambda_list <- secsse::create_lambda_list(state_names = c(0, 1),
                                          num_concealed_states = 2,
                                          transition_matrix = spec_matrix,
                                          model = "CR")
lambda_list

## ----CR_mu--------------------------------------------------------------------
mu_vec <- secsse::create_mu_vector(state_names = c(0, 1),
                                   num_concealed_states = 2,
                                   model = "CR",
                                   lambda_list = lambda_list)
mu_vec

## ----CR_Q---------------------------------------------------------------------
shift_matrix <- c()
shift_matrix <- rbind(shift_matrix, c(0, 1, 3))
shift_matrix <- rbind(shift_matrix, c(1, 0, 4))

q_matrix <- secsse::create_q_matrix(state_names = c(0, 1),
                                    num_concealed_states = 2,
                                    shift_matrix = shift_matrix,
                                    diff.conceal = TRUE)
q_matrix

## ----CR_ML--------------------------------------------------------------------
idparsopt <- 1:6 # our maximum rate parameter was 6
idparsfix <- c(0) # we want to keep al zeros at zero
initparsopt <- rep(0.1, 6)
initparsfix <- c(0.0) # all zeros remain at zero.
sampling_fraction <- c(1, 1)

idparslist <- list()
idparslist[[1]] <- lambda_list
idparslist[[2]] <- mu_vec
idparslist[[3]] <- q_matrix

answ <- secsse::cla_secsse_ml(phy = focal_tree,
                              traits = focal_traits$trait,
                              num_concealed_states = 2,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = sampling_fraction,
                              verbose = FALSE,
                              num_threads = 4)
ML_CR <- answ$ML
CR_par <- secsse::extract_par_vals(idparslist, answ$MLpars)
ML_CR
CR_par
spec_rate <- CR_par[1]
ext_rate <-  CR_par[2]
Q_Examined <- CR_par[3:4]
Q_Concealed <- CR_par[5:6]
spec_rate
ext_rate
Q_Examined
Q_Concealed

## ----AIC----------------------------------------------------------------------
res <- data.frame(ll = c(ML_ETD, ML_CTD, ML_CR),
                  k  = c(8, 8, 6),
                  model = c("ETD", "CTD", "CR"))
res$AIC <- 2 * res$k - 2 * res$ll
res

