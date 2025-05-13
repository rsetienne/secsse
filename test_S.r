library(secsse)
library(DDD)
setwd("/Users/thijsjanzen/Downloads/3_state_biwoody")
source("2b_woody_prep_biwoody_TJ.R")

all_trees <- readRDS("multi_trees_woody_3s_biwoody.rds")
all_traits <- readRDS("multi_traits_woody_3s_biwoody.rds")
all_sf <- readRDS("multi_sf_woody_3s_biwoody.rds")
all_rs <- readRDS("multi_rs_woody_3s_biwoody.rds")

for (i in 1:length(all_traits)) {
  v <- all_traits[[i]]
  v[v == 0] <- "CW"
  v[v == 1] <- "HE"
  v[v == 2] <- "IW"
  all_traits[[i]] <- v
}


param_grid <- expand.grid(model = c("ETD", "CTD", "CR", "M1", "M2", "M3"),
                          repl = 1:20)

#args <- commandArgs(trailingOnly = TRUE)

#sim_number <- as.numeric(args[[1]])

sim_number <- 60

focal_model <- param_grid$model[[sim_number]]

set.seed(sim_number)


res <- do_analysis(phy_set = all_trees[[10]],
                   trait_set = all_traits[[10]],
                   sampling_set = all_sf[[10]],
                   root_state_set = all_rs[[10]],
                   optimmethod = "simplex",
                   model = focal_model,
                   num_cycles = 1,
                   num_threads = 10,
                   verbose = TRUE)


for (i in 1:length(all_trees)) {
  focal_tree <- all_trees[[i]]
  cat(i, length(focal_tree$tip.label), "\n")
}
