## ------------------------------------------------------------------------
rm(list = ls())
library(secsse)

## ------------------------------------------------------------------------
data(traitinfo)
trait <- traitinfo
tail(trait)

## ------------------------------------------------------------------------
data("phylo_Vign")

## ------------------------------------------------------------------------
traits <- sortingtraits(trait,phylo_Vign)

## ------------------------------------------------------------------------
library(geiger)
#making sure that the first line is identified as containing header info:
rownames(trait) <- trait[,1]
#pick out all elements that do not agree between tree and data
mismat <- name.check(phylo_Vign,trait)
#this will call all taxa that are in the tree, but not the data file
#mismat$tree_not_data
#and conversely,
#mismat$data_not_tree

## ------------------------------------------------------------------------
#First we have to define idparslist, as well as, again, a user-specified value for the number of concealed states to be assessed by SecSSE.

idparslist <- id_paramPos(traits, num_concealed_states = 3)

#Let's take a look at the full all-free model by now simply typing

idparslist

## ------------------------------------------------------------------------
#idparslist[[1]][c(5,6)] <- 5 

## ------------------------------------------------------------------------
#idparslist[[2]][c(1:9)] <- 7

## ------------------------------------------------------------------------
diag(idparslist[[3]]) <- NA

## ------------------------------------------------------------------------
idparslist[[3]][1,c(5,6,8,9)] <- 0
idparslist[[3]][2,c(4,6,7,9)] <- 0
idparslist[[3]][3,c(4,5,7,8)] <- 0
idparslist[[3]][4,c(2,3,8,9)] <- 0
idparslist[[3]][5,c(1,3,7,9)] <- 0
idparslist[[3]][6,c(1,2,7,8)] <- 0
idparslist[[3]][7,c(2,3,5,6)] <- 0
idparslist[[3]][8,c(1,3,4,6)] <- 0
idparslist[[3]][9,c(1,2,4,5)] <- 0

## ------------------------------------------------------------------------
idparslist

## ------------------------------------------------------------------------
idparslist[[3]][1,c(2)] <- 19
idparslist[[3]][1,c(3)] <- 20
idparslist[[3]][1,c(4)] <- 21
idparslist[[3]][1,c(7)] <- 22
idparslist[[3]][1,c(5,6,8,9)] <- 0
idparslist[[3]][2,c(1)] <- 23
idparslist[[3]][2,c(3)] <- 24
idparslist[[3]][2,c(5)] <- 25
idparslist[[3]][2,c(8)] <- 26
idparslist[[3]][2,c(4,6,7,9)] <- 0
idparslist[[3]][3,c(1)] <- 27
idparslist[[3]][3,c(2)] <- 28
idparslist[[3]][3,c(6)] <- 29
idparslist[[3]][3,c(9)] <- 30
idparslist[[3]][3,c(4,5,7,8)] <- 0
idparslist[[3]][4,c(1)] <- 31
idparslist[[3]][4,c(5)] <- 32
idparslist[[3]][4,c(6)] <- 33
idparslist[[3]][4,c(7)] <- 34
idparslist[[3]][4,c(2,3,8,9)] <- 0
idparslist[[3]][5,c(2)] <- 35
idparslist[[3]][5,c(4)] <- 36
idparslist[[3]][5,c(6)] <- 37
idparslist[[3]][5,c(8)] <- 38
idparslist[[3]][5,c(1,3,7,9)] <- 0
idparslist[[3]][6,c(3)] <- 39
idparslist[[3]][6,c(4)] <- 40
idparslist[[3]][6,c(5)] <- 41
idparslist[[3]][6,c(9)] <- 42
idparslist[[3]][6,c(1,2,7,8)] <- 0
idparslist[[3]][7,c(1)] <- 43
idparslist[[3]][7,c(4)] <- 44
idparslist[[3]][7,c(8)] <- 45
idparslist[[3]][7,c(9)] <- 46
idparslist[[3]][7,c(2,3,5,6)] <- 0
idparslist[[3]][8,c(2)] <- 47
idparslist[[3]][8,c(5)] <- 48
idparslist[[3]][8,c(7)] <- 49
idparslist[[3]][8,c(9)] <- 50
idparslist[[3]][8,c(1,3,4,6)] <- 0
idparslist[[3]][9,c(3)] <- 51
idparslist[[3]][9,c(6)] <- 52
idparslist[[3]][9,c(7)] <- 53
idparslist[[3]][9,c(8)] <- 54
idparslist[[3]][9,c(1,2,4,5)] <- 0
diag(idparslist[[3]]) <- NA

## ------------------------------------------------------------------------
idparslist

## ------------------------------------------------------------------------
initparsopt <- c(rep(1.2,9), rep(0.1,9), rep(0.25,36))

## ------------------------------------------------------------------------
idparslist

## ------------------------------------------------------------------------
idparsopt <- c(1:9)

## ------------------------------------------------------------------------
#this would optimize speciation and extinction in the above setup
#idparsopt <- c(1:18)

## ------------------------------------------------------------------------
idparslist[[2]][] <- 10
idparslist[[3]][1,c(2,3,4,7)] <- 11
idparslist[[3]][1,c(5,6,8,9)] <- 0
idparslist[[3]][2,c(1,3,5,8)] <- 11
idparslist[[3]][2,c(4,6,7,9)] <- 0
idparslist[[3]][3,c(1,2,6,9)] <- 11
idparslist[[3]][3,c(4,5,7,8)] <- 0
idparslist[[3]][4,c(1,5,6,7)] <- 11
idparslist[[3]][4,c(2,3,8,9)] <- 0
idparslist[[3]][5,c(2,4,6,8)] <- 11
idparslist[[3]][5,c(1,3,7,9)] <- 0
idparslist[[3]][6,c(3,4,5,9)] <- 11
idparslist[[3]][6,c(1,2,7,8)] <- 0
idparslist[[3]][7,c(1,4,8,9)] <- 11
idparslist[[3]][7,c(2,3,5,6)] <- 0
idparslist[[3]][8,c(2,5,7,9)] <- 11
idparslist[[3]][8,c(1,3,4,6)] <- 0
idparslist[[3]][9,c(3,6,7,8)] <- 11
idparslist[[3]][9,c(1,2,4,5)] <- 0
diag(idparslist[[3]]) <- NA

## ------------------------------------------------------------------------
idparsopt <- c(1:9,11)

## ------------------------------------------------------------------------
idparsfix <- c(0,10)

## ------------------------------------------------------------------------
parsfix <- c(0,0.0001)

## ------------------------------------------------------------------------
library(DDD)
startingpoint <- bd_ML(brts = ape::branching.times(phylo_Vign))
intGuessLamba <- startingpoint$lambda0
intGuessMu <- startingpoint$mu0
#Make sure that the dimensions of initparsopt agree with those of idparslist and idparsopt, especially in the case of the initial guesses for rates supplied to the Q matrix. Rules of thumb are that if n=number of examined states, both intGuessLamba and intGuessMu should be replicated 2n times, and (intGuessLamba/5) should be replicated (2n)2/2 times:
initparsopt <- c(rep(intGuessLamba,9), rep((intGuessLamba/5),1))

## ------------------------------------------------------------------------
#secsse_ml(phylo_Vign,traits, num_concealed_states=3, idparslist, idparsopt, initparsopt, idparsfix, parsfix, cond="maddison_cond",root_state_weight = "maddison_weights", tol = c(1e-04, 1e-05, 1e-07), sampling_fraction=c(1,1,1), maxiter = 1000 * round((1.25)^length(idparsopt)), use_fortran=TRUE,methode="ode45", optimmethod = "simplex", num_cycles = 1, run_parallel=FALSE)

## ------------------------------------------------------------------------
#out<-secsse_ml(phylo_Vign,traits, num_concealed_states=3, idparslist, idparsopt, initparsopt, idparsfix, parsfix, cond="maddison_cond",root_state_weight = "maddison_weights", sampling_fraction=c(1,1,1), tol = c(1e-04, 1e-05, 1e-07), maxiter = 1000 * round((1.25)^length(idparsopt)), use_fortran=TRUE,methode="ode45", optimmethod = "simplex", num_cycles = 1, run_parallel=FALSE)
#saveRDS(out, file="output.RDS")

## ------------------------------------------------------------------------
#readRDS("output.RDS")

## ------------------------------------------------------------------------
#$MLpars[[1]]
#          1A           2A           3A           1B           2B           3B 
#4.842634e-16 1.080409e-01 7.843821e-02 4.029147e-09 3.018863e-02 3.018863e-02 

#$MLpars[[2]]
#         1A          2A          3A          1B          2B          3B 
#0.002000000 0.002000109 0.002734071 0.001988593 0.002169052 0.003969142 

#$MLpars[[3]]
#     1A   2A   3A   1B   2B   3B
#1A   NA 0.01 0.01 0.01 0.01 0.01
#2A 0.01   NA 0.01 0.01 0.01 0.01
#3A 0.01 0.01   NA 0.01 0.01 0.01
#1B 0.01 0.01 0.01   NA 0.01 0.01
#2B 0.01 0.01 0.01 0.01   NA 0.01
#3B 0.01 0.01 0.01 0.01 0.01   NA


#$ML
#[1] -848.0895

## ------------------------------------------------------------------------
masterBlock<-matrix(99,ncol=3,nrow=3,byrow=T) 

## ------------------------------------------------------------------------
diag(masterBlock) <- NA
masterBlock[1,2] <- 6
masterBlock[1,3] <- 7

masterBlock[2,1] <- 8
masterBlock[2,3] <- 9

masterBlock[3,1] <- 10
masterBlock[3,2] <- 11


## ------------------------------------------------------------------------
diff.conceal <- FALSE

## ------------------------------------------------------------------------
myQ<-q_doubletrans(traits,masterBlock,diff.conceal)
idparslist[[3]] <- myQ

## ------------------------------------------------------------------------
idparslist[[3]]

## ------------------------------------------------------------------------
diff.conceal <- TRUE
myQ <- q_doubletrans(traits,masterBlock,diff.conceal)
idparslist[[3]] <- myQ
idparslist[[3]]

## ------------------------------------------------------------------------
#shareFactors <- c(.1,.2)

## ------------------------------------------------------------------------
#initFactors <- c(1,1)

## ------------------------------------------------------------------------
# diag(masterBlock) <- NA
# masterBlock[1,2] <- 6
# masterBlock[1,3] <- 6.1  #factor 1: lobed to palmate
# 
# masterBlock[2,1] <- 7
# masterBlock[2,3] <- 8
# 
# masterBlock[3,1] <- 7.2  #factor 2: palmate to lobed
# masterBlock[3,2] <- 9

## ------------------------------------------------------------------------
#secsse_ml_struc(phylo_Vign..., shareFactors, initFactors)

## ------------------------------------------------------------------------
#       traits traits traits
# [1,]      2      2      2
# [2,]      1      1      1
# [3,]      2      2      2
# [4,]      3      1      1
# [5,]      1      2      3

