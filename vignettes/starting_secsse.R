## -----------------------------------------------------------------------------
library(secsse)
data(traits)
tail(traits)

## -----------------------------------------------------------------------------
data("phylo_vignette")

## -----------------------------------------------------------------------------
sorted_traits <- sortingtraits(traits, phylo_vignette)

## -----------------------------------------------------------------------------
library(geiger)
#pick out all elements that do not agree between tree and data
mismat <- name.check(phylo_vignette, traits)
#this will call all taxa that are in the tree, but not the data file
#mismat$tree_not_data
#and conversely,
#mismat$data_not_tree

## ----plot_tree----------------------------------------------------------------
if (requireNamespace("diversitree")) {
  for_plot <- data.frame(trait = traits$trait, row.names = phylo_vignette$tip.label)
diversitree::trait.plot(phylo_vignette, dat = for_plot, 
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

answ <- secsse::cla_secsse_ml(phy = phylo_vignette,
                              traits = traits$trait,
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

answ <- secsse::cla_secsse_ml(phy = phylo_vignette,
                              traits = traits$trait,
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

answ <- secsse::cla_secsse_ml(phy = phylo_vignette,
                              traits = traits$trait,
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

