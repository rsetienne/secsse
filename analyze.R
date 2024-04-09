setwd("/Users/thijsjanzen/MEGAsync2/lucas/sim_mar_2024/complete_stage2/")
require(secsse)
source("secsse_prep.R")
source("setup_models_15.R")
#args = commandArgs(trailingOnly = TRUE)
sim_number = 23 #as.numeric(args[[1]])

phy <- readRDS(paste0("tree_complete_23.rds"))
phylotree = phy$phy
traits <- phy$obs_traits

new_states <- c("S", "N", "E", "A",
                "SN", "NE", "EA", "SE", "NA", "SA",
                "SNE", "SEA", "SNA", "NEA",
                "SNEA")

old_states <-  c("N", "S", "E", "A", "SN", "NE", "EA", "NEA", "NEA")

sim_traits <- c()
for (i in 1:length(phy$obs_traits)) {
  sim_traits[i] <- which(new_states == phy$obs_traits[i])
}

num_concealed_states <- 15

file_name <- paste0("complete_", sim_number, ".txt")

model_order <- c(33, 28)

require(readr)
if (file.exists(file_name)) {
  cat("found previous results\n")
  require(readr)
  vx <- read_delim(file_name,
                   delim = " ",
                   show_col_types = FALSE,
                   col_names = FALSE)

  vv <- table(vx$X1)
  to_remove <- which(vv > 2)
  to_remove <- to_remove[which(names(to_remove) %in% model_order)]
  if (length(to_remove > 0)) {
    for (i in 1:length(to_remove)) {
      model_order <- model_order[-which(model_order == names(to_remove)[i])]
    }
  }
  if (length(model_order) < 1) {
    cat("This is done")
    stop("Found results from two previous fits")
  }
}

model_order <- rev(model_order)

for (model in model_order) {
  lambdas <- create_lambda_model(model = model,
                                 num_concealed_states = num_concealed_states)

  mus_info <- create_mu_vector_local(model, lambdas,
                                     num_concealed_states = num_concealed_states)
  mus <- mus_info$mu_vec

  mu_vals <- mus_info$start_rate:mus_info$end_rate
  mu1 <- mu_vals[1]
  mu2 <- mu_vals[2]
  mu3 <- mu_vals[3]
  mu4 <- mu_vals[4] # not all of these are actually used!

  q_vals <- mus_info$end_rate + 1:4
  qu1 <- q_vals[1]
  qu2 <- q_vals[2]
  qu3 <- q_vals[3]
  qu4 <- q_vals[4] # not all of these are actually used!

  q_mat <- create_q_mat_local(model,
                              q1 = qu1, q2 = qu2, q3 = qu3, q4 = qu4,
                              m1 = mu1, m3 = mu3, m4 = mu4,
                              num_concealed_states = num_concealed_states,
                              diff_conceal = FALSE)

  current_max_rate <- max(q_mat, na.rm = TRUE)
  # here we have to add the double transitions!

  # SE -> SNE
  x <- which(rownames(q_mat) == "SE")
  y <- which(rownames(q_mat) == "SNE")
  rate1 <- current_max_rate + 1
  q_mat[x, y] <- rate1

  x <- which(rownames(q_mat) == "SEA")
  y <- which(rownames(q_mat) == "SNEA")
  q_mat[x, y] <- rate1

  # we model q_SN + q_EN
  s_index <- which(rownames(q_mat) == "S")
  sn_index <- which(rownames(q_mat) == "SN")
  q_SN <- q_mat[s_index, sn_index]

  e_index <- which(rownames(q_mat) == "E")
  en_index <- which(rownames(q_mat) == "NE")
  q_EN <- q_mat[e_index, en_index]

  functions_defining_params <- list()
  functions_defining_params[[1]] <- function() {
    tt <-  paste0("par_", eval(rate1), " <- par_", eval(q_SN), " + par_", eval(q_EN))
    eval(parse(text = tt))
  }

  rate2 <- rate1 + 1
  # rate2 = q_NE + q_AE
  x <- which(rownames(q_mat) == "NA")
  y <- which(rownames(q_mat) == "NEA")
  q_mat[x, y] <- rate2

  x <- which(rownames(q_mat) == "SNA")
  y <- which(rownames(q_mat) == "SNEA")
  q_mat[x, y] <- rate2

  n_index <- which(rownames(q_mat) == "N")
  ne_index <- which(rownames(q_mat) == "NE")
  q_NE <- q_mat[n_index, ne_index]

  a_index <- which(rownames(q_mat) == "A")
  ea_index <- which(rownames(q_mat) == "EA")
  q_AE <- q_mat[a_index, ea_index]

  functions_defining_params[[2]] <- function() {
    tt <- paste0("par_", eval(rate2), " <- par_", eval(q_NE), " + par_", eval(q_AE))
    eval(parse(text = tt))
  }

  # and now for the lambdas:
  rate3 <- rate2 + 1
  rate4 <- rate3 + 1
  for (i in 1:length(lambdas)) {
    for (j in 1:nrow(lambdas[[i]])) {
      for (k in 1:ncol(lambdas[[i]])) {
        if (lambdas[[i]][j, k] == 102) lambdas[[i]][j, k] <- rate3
        if (lambdas[[i]][j, k] == 103) lambdas[[i]][j, k] <- rate4
      }
    }
  }

  states <- c("S", "N", "E", "A",
              "SN", "NE", "EA", "SE", "NA", "SA",
              "SNE", "SEA", "SNA", "NEA",
              "SNEA")

  index_SN <- which(states == "SN")
  index_S <- which(states == "S")
  index_N <- which(states == "N")
  l2 <- lambdas[[index_SN]][index_S, index_N]

  # and now the custom functions:
  functions_defining_params[[3]] <- function() {
    tt <- paste0("par_", eval(rate3), " <- par_", eval(l2), " * 2")
    eval(parse(text = tt))
  }

  # and now the custom functions:
  functions_defining_params[[4]] <- function() {
    tt <- paste0("par_", eval(rate4), " <- par_", eval(l2), " * 3")
    eval(parse(text = tt))
  }

  idparsfuncdefpar <- c(rate1, rate2, rate3, rate4)

  states <- c("S", "N", "E", "A",
              "SN", "NE", "EA", "SE", "NA", "SA",
              "SNE", "SEA", "SNA", "NEA",
              "SNEA")

  q_mat_double <- secsse::q_doubletrans(traits = states,
                                        masterBlock = q_mat,
                                        diff.conceal = FALSE)

  root_state_weight <- "proper_weights"
  num_traits <- length(new_states)
  sampling_fraction <- rep(0, num_traits) # c(1, 1, 1, 1, 1, 1, 1, 1)
  for (i in 1:length(traits)) {
    index <- which(traits[i] == states)
    sampling_fraction[index] <- 1
  }

  cond <- "proper_cond"

  param_posit <- list()
  param_posit[[1]] <- lambdas
  param_posit[[2]] <- mus
  param_posit[[3]] <- q_mat_double

  idparslist <- param_posit

  max_par_val <- current_max_rate

  idparsopt <- 1:max_par_val
  idparsopt_orig <- idparsopt
  idparsopt <- idparsopt[-which(idparsopt == mu2)]

  idparsfix <- c(0, mu2)
  parsfix <- c(0.0, 0.0)

  found_res <- c()
  for (r in 1:3) {
    initpars <- runif(n = length(idparsopt), min = 0, max = 1)

    answ <- try(
     secsse::cla_secsse_ml_func_def_pars(
      phy = phylotree,
      traits = traits,
      num_concealed_states = num_concealed_states,
      idparslist = idparslist,
      idparsopt = idparsopt,
      initparsopt = initpars,
      idparsfix = idparsfix,
      parsfix = parsfix,
      idparsfuncdefpar = idparsfuncdefpar,
      functions_defining_params = functions_defining_params,
      idfactorsopt = "noFactor",
      initfactors = c(0),
      cond = cond,
      root_state_weight = root_state_weight,
      sampling_fraction = sampling_fraction,
     # tol = c(1e-3, 1e-4, 1e-5),
      optimmethod = "subplex",
      num_cycles = 3,
      loglik_penalty = 0,
      is_complete_tree = TRUE,
      num_threads = 10,
      verbose = TRUE,
      atol = 1e-7,
      rtol = 1e-6),
      silent = FALSE)

    if (class(answ) != "try-error") {

      found_pars <- secsse::extract_par_vals(param_posit, answ$MLpars)
      aic <- 2 * (1 + length(idparsopt)) - 2 * as.numeric(answ$ML)
      found_pars <- found_pars[idparsopt_orig]

      to_add <- c(model, found_pars, answ$ML, aic)
      file_name <- paste0("complete_", sim_number, ".txt")

      cat(to_add, file = file_name, "\n", append = TRUE)
      cat(to_add, "\n")
      found_res <- rbind(found_res, to_add)
    }
  }
}
