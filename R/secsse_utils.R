#' @title Parameter structure setting
#' Sets the parameters (speciation, extinction and transition) ids. Needed for 
#' ML calculation ([secsse_ml()]).
#' 
#' @inheritParams default_params_doc
#' 
#' @return A list that includes the ids of the parameters for ML analysis.
#' @examples
#' traits <- sample(c(0,1,2), 45,replace = TRUE) #get some traits
#' num_concealed_states <- 3
#' param_posit <- id_paramPos(traits,num_concealed_states)
#' @export
#' @rawNamespace useDynLib(secsse, .registration = TRUE)
#' @rawNamespace import(Rcpp)
#' @rawNamespace importFrom(RcppParallel, RcppParallelLibs)
id_paramPos <- function(traits, num_concealed_states) { #noLint
    idparslist <- list()
    if (is.matrix(traits)) {
        traits <- traits[, 1]
    }

    ly <- length(sort(unique(traits))) * 2 * num_concealed_states
    d <- ly / 2
    idparslist[[1]] <- 1:d
    idparslist[[2]] <- (d + 1):ly
    toMatrix <- 1
    matPos <- (ly + 1):(((d^2) - d) + d * 2)
    for (i in 1:d) {
        toMatrix <- c(toMatrix,
                      matPos[(i * d - (d - 1)):((i * d - (d - 1)) + d)])
    }
    toMatrix <- toMatrix[1:d^2]
    Q <- matrix(toMatrix, ncol = d, nrow = d, byrow = TRUE)
    diag(Q) <- NA
    idparslist[[3]] <- Q

    lab_states <- rep(as.character(sort(unique(traits))), num_concealed_states)

    lab_conceal <- NULL
    for (i in 1:num_concealed_states) {

        lab_conceal <- c(lab_conceal,
                         rep(LETTERS[i],
                             length(sort(unique(traits)))))
    }

    statesCombiNames <- character()
    for (i in seq_along(lab_states)) {
        statesCombiNames <- c(statesCombiNames,
                              paste0(lab_states[i],
                                     lab_conceal[i]))
    }
    colnames(idparslist[[3]]) <- statesCombiNames
    rownames(idparslist[[3]]) <- statesCombiNames
    names(idparslist) <- c("lambdas", "mus", "Q")
    names(idparslist[[1]]) <- statesCombiNames
    names(idparslist[[2]]) <- statesCombiNames
    return(idparslist)
}

#' @keywords internal
create_q_matrix_int <- function(masterBlock,
                                concealnewQMatr,
                                ntraits,
                                diff.conceal) {
    Q <- NULL
    for (i in 1:ntraits) {
        Qrow <- NULL
        for (ii in 1:ntraits) {
            entry <- masterBlock[i, ii]
            if (is.na(entry)) {
                Qrow <- cbind(Qrow, masterBlock)
            } else {
                if (diff.conceal == TRUE) {
                    entry <- concealnewQMatr[i, ii]
                }

                outDiagBlock <- matrix(0,
                                       ncol = ntraits,
                                       nrow = ntraits,
                                       byrow = TRUE)
                diag(outDiagBlock) <- entry
                Qrow <- cbind(Qrow, outDiagBlock)
            }
        }
        Q <- rbind(Q, Qrow)
    }
    return(Q)
}


#' @title Basic Qmatrix
#' Sets a Q matrix where double transitions are not allowed
#' 
#' @inheritParams default_params_doc
#' 
#' @return Q matrix that includes both examined and concealed states, it should
#' be declared as the third element of idparslist.
#' @description This function expands the Q_matrix, but it does so assuming
#' that the number of concealed traits is equal to the number of examined
#' traits, if you have a different number, you should consider looking at
#' the function [expand_q_matrix()].
#' @examples
#' traits <- sample(c(0,1,2), 45,replace = TRUE) #get some traits
#' # For a three-state trait
#' masterBlock <- matrix(99,ncol = 3,nrow = 3,byrow = TRUE)
#' diag(masterBlock) <- NA
#' masterBlock[1,2] <- 6
#' masterBlock[1,3] <- 7
#' masterBlock[2,1] <- 8
#' masterBlock[2,3] <- 9
#' masterBlock[3,1] <- 10
#' masterBlock[3,2] <- 11
#' myQ <- q_doubletrans(traits,masterBlock,diff.conceal = FALSE)
#' # now, it can replace the Q matrix from id_paramPos
#' num_concealed_states <- 3
#' param_posit <- id_paramPos(traits,num_concealed_states)
#' param_posit[[3]] <- myQ
#' @export
q_doubletrans <- function(traits, masterBlock, diff.conceal) {
    if (diff.conceal == TRUE &&
        all(floor(masterBlock) == masterBlock, na.rm = TRUE) == FALSE) {
        integersmasterBlock <- floor(masterBlock)
        factorBlock <- signif(masterBlock - integersmasterBlock, digits = 2)

        factorstoExpand <- unique(sort(c(factorBlock)))
        factorstoExpand <- factorstoExpand[factorstoExpand > 0]
        newshareFac <-
            (max(factorstoExpand * 10) + 1):(max(factorstoExpand * 10) +
                                                 length(factorstoExpand))
        newshareFac <- newshareFac / 10

        for (iii in seq_along(newshareFac)) {
            factorBlock[which(factorBlock == factorstoExpand[iii])] <-
                newshareFac[iii]
        }

        ntraits <- length(sort(unique(traits)))
        uniqParQ <- sort(unique(c(floor(masterBlock))))
        uniqParQ2 <- uniqParQ[which(uniqParQ > 0)]
        concealnewQ <- (max(uniqParQ2) + 1):(max(uniqParQ2) + length(uniqParQ2))

        for (iii in seq_along(concealnewQ)) {
            integersmasterBlock[which(integersmasterBlock == uniqParQ2[iii])] <-
                concealnewQ[iii]
        }
        concealnewQMatr <- integersmasterBlock + factorBlock

        Q <- create_q_matrix_int(masterBlock,
                                 concealnewQMatr,
                                 ntraits,
                                 diff.conceal)
    } else {
        ntraits <- length(sort(unique(traits)))
        uniqParQ <- sort(unique(c(masterBlock)))
        uniqParQ2 <- uniqParQ[which(uniqParQ > 0)]
        concealnewQ <- (max(uniqParQ2) + 1):(max(uniqParQ2) + length(uniqParQ2))
        concealnewQMatr <- masterBlock
        for (I in seq_along(uniqParQ2)) {
            uniqParQ2
            concealnewQMatr[concealnewQMatr == uniqParQ2[I]] <- concealnewQ[I]
        }

        Q <- create_q_matrix_int(masterBlock,
                                 concealnewQMatr,
                                 ntraits,
                                 diff.conceal)
    }
    uniq_traits <- unique(traits)
    uniq_traits <- uniq_traits[!is.na(uniq_traits)]
    all_names <- get_state_names(state_names = uniq_traits,
                                 num_concealed_states = length(uniq_traits))
    colnames(Q) <- all_names
    rownames(Q) <- all_names
    return(Q)
}


#' @title Data checking and trait sorting
#' In preparation for likelihood calculation, it orders trait data according
#' the tree tips
#' 
#' @inheritParams default_params_doc
#' 
#' @return Vector of traits
#' @examples
#' # Some data we have prepared
#' data(traits)
#' data('phylo_vignette')
#' traits <- sortingtraits(traits, phylo_vignette)
#' @export
sortingtraits <- function(trait_info, phy) {
    trait_info <- as.matrix(trait_info)
    if (length(phy$tip.label) != nrow(trait_info)) {
        stop("Number of species in the tree must be the same as
             in the trait file")
    }

    if (identical(as.character(sort(phy$tip.label)),
                  as.character(sort(trait_info[, 1]))) == FALSE) {
        mismatch <- match(as.character(sort(trait_info[, 1])),
                          as.character(sort(phy$tip.label)))
        mismatched <- (sort(trait_info[, 1]))[which(is.na(mismatch))]
        stop(
            paste(c("Mismatch on tip labels and taxa names, check the species:",
                    mismatched), collapse = " ")
        )
    }

    trait_info <- trait_info[match(phy$tip.label, trait_info[, 1]), ]
    trait_info[, 1] == phy$tip.label

    if (ncol(trait_info) == 2) {
        traits <- as.numeric(trait_info[, 2])
    }

    if (ncol(trait_info) > 2) {
        traits <- NULL
        for (i in 1:(ncol(trait_info) - 1)) {
            traits <- cbind(traits, as.numeric(trait_info[, 1 + i]))
        }
    }
    return(traits)
}

#' @title Parameter structure setting for cla_secsse
#' It sets the parameters (speciation, extinction and transition)
#' IDs. Needed for ML calculation with cladogenetic options (cla_secsse_ml)
#' 
#' @inheritParams default_params_doc
#' 
#' @return A list that includes the ids of the parameters for ML analysis.
#' @examples
#'traits <- sample(c(0,1,2), 45,replace = TRUE) #get some traits
#'num_concealed_states <- 3
#'param_posit <- cla_id_paramPos(traits, num_concealed_states)
#' @export
cla_id_paramPos <- function(traits, num_concealed_states) {
    idparslist <- list()
    if (is.matrix(traits)) {
        traits <- traits[, 1]
    }

    ly <- length(sort(unique(traits))) * 2 * num_concealed_states
    d <- ly / 2
    toMatrix <- 1
    matPos <- (ly + 1):(((d^2) - d) + d * 2)
    for (i in 1:d) {
        toMatrix <- c(toMatrix,
                      matPos[(i * d - (d - 1)):((i * d - (d - 1)) + d)])
    }
    toMatrix <- toMatrix[1:d^2]
    Q <- matrix(toMatrix, ncol = d, nrow = d, byrow = TRUE)
    diag(Q) <- NA
    lab_states <- rep(as.character(sort(unique(traits))), num_concealed_states)

    lab_conceal <- NULL
    for (i in 1:num_concealed_states) {
        lab_conceal <- c(lab_conceal,
                         rep(LETTERS[i],
                             length(sort(unique(traits)))))
    }

    statesCombiNames <- character()
    for (i in seq_along(lab_states)) {
        statesCombiNames <- c(statesCombiNames,
                              paste0(lab_states[i],
                                     lab_conceal[i]))
    }

    idparslist[[1]] <- matrix(0, ncol = d, nrow = 4)
    idparslist[[2]] <- (d + 1):ly
    idparslist[[3]] <- Q

    rownames(idparslist[[1]]) <- c("dual_inheritance",
                                   "single_inheritance",
                                   "dual_symmetric_transition",
                                   "dual_asymmetric_transition")

    colnames(idparslist[[1]]) <- statesCombiNames
    colnames(idparslist[[3]]) <- statesCombiNames
    rownames(idparslist[[3]]) <- statesCombiNames
    names(idparslist) <- c("lambdas", "mus", "Q")
    names(idparslist[[2]]) <- statesCombiNames
    return(idparslist)
}

#' @title Prepares the entire set of lambda matrices for cla_secsse.
#' It provides the set of matrices containing all the speciation rates
#' 
#' @inheritParams default_params_doc
#' 
#' @return A list of lambdas, its length would be the same than the number of
#' trait states * num_concealed_states..
#' @export
#' @examples
#' set.seed(13)
#' phylotree <- ape::rcoal(12, tip.label = 1:12)
#' traits <- sample(c(0, 1, 2),
#'                  ape::Ntip(phylotree), replace = TRUE)
#' num_concealed_states <- 3
#' # the idparlist for a ETD model (dual state inheritance model of evolution)
#' # would be set like this:
#' idparlist <- secsse::cla_id_paramPos(traits, num_concealed_states)
#' lambd_and_modeSpe <- idparlist$lambdas
#' lambd_and_modeSpe[1, ] <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
#' idparlist[[1]] <- lambd_and_modeSpe
#' idparlist[[2]][] <- 0
#' masterBlock <- matrix(4, ncol = 3, nrow = 3, byrow = TRUE)
#' diag(masterBlock) <- NA
#' idparlist[[3]] <- q_doubletrans(traits, masterBlock, diff.conceal = FALSE)
#' # Now, internally, clasecsse sorts the lambda matrices, so they look like
#' #  a list with 9 matrices, corresponding to the 9 states
#' # (0A,1A,2A,0B, etc)
#'
#' parameter <- idparlist
#' lambda_and_modeSpe <- parameter$lambdas
#' lambda_and_modeSpe[1, ] <- c(0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 0.01, 0.01, 0.01)
#' parameter[[1]] <- prepare_full_lambdas(traits, num_concealed_states,
#'                                        lambda_and_modeSpe)
prepare_full_lambdas <- function(traits,
                                 num_concealed_states,
                                 lambd_and_modeSpe) {
    if (is.list(lambd_and_modeSpe)) return(lambd_and_modeSpe)

    num_exami <- length(sort(unique(unlist(traits))))
    mat_size <- num_exami * num_concealed_states
    posib_trans <- matrix(1,
                          ncol = num_exami,
                          nrow = num_exami,
                          byrow = TRUE)
    diag(posib_trans) <- NA
    posib_trans <- q_doubletrans(unique(unlist(traits)),
                                 masterBlock = posib_trans,
                                 diff.conceal = FALSE)

    full_lambdas <- list()
    for (jj in 1:mat_size) {
        # dual_state_inhe
        m1 <- matrix(0, ncol = mat_size, nrow = mat_size)
        m1[jj, jj] <- as.numeric(lambd_and_modeSpe[, jj][1])

        # single_state_inhe
        m2 <- matrix(0, ncol = mat_size, nrow = mat_size)
        m2[, jj] <- posib_trans[jj, ]
        m2[jj, jj] <- 0
        m2[m2 == 1] <- as.numeric(lambd_and_modeSpe[, jj][2])
        # symet_state_emerge

        m3 <- matrix(0, ncol = mat_size, nrow = mat_size)

        diag(m3) <- posib_trans[jj, ]
        m3[jj, jj] <- 0
        m3[m3 == 1] <- as.numeric(lambd_and_modeSpe[, jj][3])
        # symet_state_emerge

        m4 <- matrix(0, ncol = mat_size, nrow = mat_size)
        for (i in seq_along(which(posib_trans[jj, ] == 1))) {
            m4[which(posib_trans[jj, ] == 1)[i], ] <- posib_trans[jj, ]
        }
        m4[, jj] <- 0
        m4[upper.tri(m4)] <- 0
        diag(m4) <- 0
        m4[is.na(m4)] <- 0
        m4[m4 == 1] <- as.numeric(lambd_and_modeSpe[, jj][4])
        full_lambdas[[jj]] <- m1 + m2 + m3 + m4
    }
    return(full_lambdas)
}

#' @keywords internal
penalty <- function(pars, loglik_penalty = 0) {
  if (loglik_penalty == 0) return(0)
  pars <- unlist(unlist(pars))
  return(loglik_penalty * sum(pars^2) / (2 * length(pars)))
}

#' @keywords internal
calc_mus <- function(is_complete_tree,
                     idparslist,
                     idparsfix,
                     parsfix,
                     idparsopt,
                     initparsopt) {
    mus <- NULL
    if (is_complete_tree) {
        mus <- rep(NA, length(idparslist[[2]]))
        for (i in seq_along(idparslist[[2]])) {
            mus[i] <- c(parsfix[which(idparsfix == idparslist[[2]][i])],
                        initparsopt[which(idparsopt == idparslist[[2]][i])])
        }
    }
    return(mus)
}

#' @keywords internal
check_tree <- function(phy, is_complete_tree) {
    if (ape::is.rooted(phy) == FALSE) {
        stop("The tree needs to be rooted.")
    }

    if (ape::is.binary(phy) == FALSE) {
        stop("The tree needs to be fully resolved.")
    }
    # using option = 2, which uses the variance, the default method until ape
    # 3.5. This seems to be less sensitive.
    if (ape::is.ultrametric(phy, option = 2) == FALSE && 
        is_complete_tree == FALSE) {
        stop("The tree needs to be ultrametric.")
    }
    if (any(phy$edge.length == 0)) {
      stop("The tree must have internode distancs that are all larger than 0.")
    }
}

#' @keywords internal
check_traits <- function(traits, sampling_fraction) {
    if (is.matrix(traits)) {
        if (length(sampling_fraction) != length(sort(unique(traits[, 1])))) {
            stop("Sampling_fraction must have as many elements 
           as the number of traits.")
        }

        if (all(sort(unique(as.vector(traits))) == sort(unique(traits[, 1]))) ==
            FALSE) {
            stop(
                "Check your trait argument; if you have more than one column,
        make sure all your states are included in the first column."
            )
        }
    } else {
        if (length(sampling_fraction) != length(sort(unique(traits)))) {
            stop("Sampling_fraction must have as many elements as 
           the number of traits.")
        }
    }

    if (length(sort(unique(as.vector(traits)))) < 2) {
        stop("The trait has only one state.")
    }
}

#' @keywords internal
check_root_state_weight <- function(root_state_weight, traits) {
    if (is.numeric(root_state_weight)) {
        #if (length(root_state_weight) != length(sort(unique(traits)))) {
        #    stop("There need to be as many elements in root_state_weight 
        #   as there are traits.")
        #}
        if (length(which(root_state_weight == 1)) != 1) {
            stop("The root_state_weight needs only one 1.")
        }
      if (sum(root_state_weight) > 1) {
        stop("Root state weights need to sum to 1")
      }
    } else {
        if (any(root_state_weight == "maddison_weights" |
                root_state_weight == "equal_weights" |
                root_state_weight == "proper_weights") == FALSE) {
            stop("The root_state_weight must be any of 
           maddison_weights, equal_weights, or proper_weights.")
        }
    }
}

#' @keywords internal
check_input <- function(traits,
                        phy,
                        sampling_fraction,
                        root_state_weight,
                        is_complete_tree) {
    check_root_state_weight(root_state_weight, sampling_fraction)

    check_tree(phy, is_complete_tree)

    # check_traits(traits, sampling_fraction)
}


#' @keywords internal
transf_funcdefpar <- function(idparsfuncdefpar,
                              functions_defining_params,
                              idfactorsopt,
                              trparsfix,
                              trparsopt,
                              idparsfix,
                              idparsopt) {
    trparfuncdefpar <- NULL
    ids_all <- c(idparsfix, idparsopt)

    values_all <- c(trparsfix / (1 - trparsfix),
                    trparsopt / (1 - trparsopt))
    a_new_envir <- new.env()
    x <- as.list(values_all)  ## To declare all the ids as variables

    if (is.null(idfactorsopt)) {
        names(x) <- paste0("par_", ids_all)
    } else {
        names(x) <- c(paste0("par_", ids_all), paste0("factor_", idfactorsopt))
    }
    list2env(x, envir = a_new_envir)

    for (jj in seq_along(functions_defining_params)) {
        myfunc <- functions_defining_params[[jj]]
        environment(myfunc) <- a_new_envir
        value_func_defining_parm <- local(myfunc(), envir = a_new_envir)

        ## Now, declare the variable that is just calculated, so it is available
        ## for the next calculation if needed
        y <- as.list(value_func_defining_parm)
        names(y) <- paste0("par_", idparsfuncdefpar[jj])
        list2env(y, envir = a_new_envir)

        if (is.numeric(value_func_defining_parm) == FALSE) {
            stop("Something went wrong with the calculation of 
                 parameters in 'functions_param_struct'")
        }
        trparfuncdefpar <- c(trparfuncdefpar, value_func_defining_parm)
    }
    trparfuncdefpar <- trparfuncdefpar / (1 + trparfuncdefpar)
    rm(a_new_envir)
    return(trparfuncdefpar)
}

#' @keywords internal
update_values_transform_cla <- function(trpars,
                                        idparslist,
                                        idpars,
                                        parvals) {
    for (i in seq_along(idpars)) {
        for (j in seq_len(nrow(trpars[[3]]))) {
            id <- which(idparslist[[1]][[j]] == idpars[i])
            trpars[[1]][[j]][id] <- parvals[i]
        }
        for (j in 2:3) {
            id <- which(idparslist[[j]] == idpars[i])
            trpars[[j]][id] <- parvals[i]
        }
    }
    return(trpars)
}

#' @keywords internal
transform_params_cla <- function(idparslist,
                                 idparsfix,
                                 trparsfix,
                                 idparsopt,
                                 trparsopt,
                                 structure_func,
                                 idparsfuncdefpar,
                                 trparfuncdefpar) {
    trpars1 <- idparslist
    for (j in seq_len(nrow(trpars1[[3]]))) {
        trpars1[[1]][[j]][, ] <- NA
    }

    for (j in 2:3) {
        trpars1[[j]][] <- NA
    }

    if (length(idparsfix) != 0) {
        trpars1 <- update_values_transform_cla(trpars1,
                                               idparslist,
                                               idparsfix,
                                               trparsfix)
    }

    trpars1 <- update_values_transform_cla(trpars1,
                                           idparslist,
                                           idparsopt,
                                           trparsopt)
    ## structure_func part
    if (!is.null(structure_func)) {
        trpars1 <- update_values_transform_cla(trpars1,
                                               idparslist,
                                               idparsfuncdefpar,
                                               trparfuncdefpar)
    }

    pre_pars1 <- list()
    pars1 <- list()

    for (j in seq_len(nrow(trpars1[[3]]))) {
        pre_pars1[[j]] <- trpars1[[1]][[j]][, ] / (1 - trpars1[[1]][[j]][, ])
    }

    pars1[[1]] <- pre_pars1
    for (j in 2:3) {
        pars1[[j]] <- trpars1[[j]] / (1 - trpars1[[j]])
    }

    return(pars1)
}

#' @keywords internal
update_values_transform <- function(trpars,
                                    idparslist,
                                    idpars,
                                    parvals) {
    for (i in seq_along(idpars)) {
        for (j in 1:3) {
            id <- which(idparslist[[j]] == idpars[i])
            trpars[[j]][id] <- parvals[i]
        }
    }
    return(trpars)
}

#' @keywords internal
transform_params_normal <- function(idparslist,
                                    idparsfix,
                                    trparsfix,
                                    idparsopt,
                                    trparsopt,
                                    structure_func,
                                    idparsfuncdefpar,
                                    trparfuncdefpar) {
    trpars1 <- idparslist
    for (j in 1:3) {
        trpars1[[j]][] <- NA
    }
    if (length(idparsfix) != 0) {
        trpars1 <- update_values_transform(trpars1,
                                           idparslist,
                                           idparsfix,
                                           trparsfix)
    }

    trpars1 <- update_values_transform(trpars1,
                                       idparslist,
                                       idparsopt,
                                       trparsopt)

    ## if structure_func part
    if (is.null(structure_func) == FALSE) {
        trpars1 <- update_values_transform(trpars1,
                                           idparslist,
                                           idparsfuncdefpar,
                                           trparfuncdefpar)
    }
    pars1 <- list()
    for (j in 1:3) {
        pars1[[j]] <- trpars1[[j]] / (1 - trpars1[[j]])
    }
    return(pars1)
}

#' @keywords internal
secsse_transform_parameters <- function(trparsopt,
                                        trparsfix,
                                        idparsopt,
                                        idparsfix,
                                        idparslist,
                                        structure_func) {
    if (!is.null(structure_func)) {
        idparsfuncdefpar <- structure_func[[1]]
        functions_defining_params <- structure_func[[2]]

        if (length(structure_func[[3]]) > 1) {
            idfactorsopt <- structure_func[[3]]
        } else {
            if (structure_func[[3]] == "noFactor") {
                idfactorsopt <- NULL
            } else {
                idfactorsopt <- structure_func[[3]]
            }
        }

        trparfuncdefpar <- transf_funcdefpar(idparsfuncdefpar =
                                                 idparsfuncdefpar,
                                             functions_defining_params =
                                                 functions_defining_params,
                                             idfactorsopt = idfactorsopt,
                                             trparsfix = trparsfix,
                                             trparsopt = trparsopt,
                                             idparsfix = idparsfix,
                                             idparsopt = idparsopt)
    }

    if (is.list(idparslist[[1]])) {
        # when the ml function is called from cla_secsse
        pars1 <- transform_params_cla(idparslist,
                                      idparsfix,
                                      trparsfix,
                                      idparsopt,
                                      trparsopt,
                                      structure_func,
                                      idparsfuncdefpar,
                                      trparfuncdefpar)
    } else {
        # when non-cla option is called
        pars1 <- transform_params_normal(idparslist,
                                         idparsfix,
                                         trparsfix,
                                         idparsopt,
                                         trparsopt,
                                         structure_func,
                                         idparsfuncdefpar,
                                         trparfuncdefpar)
    }
    return(pars1)
}

condition_root_edge <- function(mergeBranch2,
                                nodeM) {
  # note that we don't do maddison conditioning here.
  
}


#' @keywords internal
condition <- function(cond,
                      mergeBranch2,
                      weight_states,
                      lambdas,
                      nodeM,
                      is_root_edge = FALSE) {

    if (cond == "no_cond") {
      return(mergeBranch2)
    }
  
    lmb <- length(mergeBranch2)
    d <- length(lambdas)
    
    if (is_root_edge) {
      return(mergeBranch2 / (1 - nodeM[1:d]))
    }
    
    if (is.list(lambdas)) {
        if (cond == "maddison_cond") {
            pre_cond <- rep(NA, lmb) # nolint
            for (j in 1:lmb) {
                pre_cond[j] <- sum(weight_states[j] *
                                       lambdas[[j]] *
                                       (1 - nodeM[1:d][j]) ^ 2)
            }
            mergeBranch2 <- mergeBranch2 / sum(pre_cond) # nolint
        }

        if (cond == "proper_cond") {
            pre_cond <- rep(NA, lmb) # nolint
            prefactor <- ((1 - nodeM[1:d]) %o% (1 - nodeM[1:d]))
            for (j in 1:lmb) {
                pre_cond[j] <- sum(lambdas[[j]] * prefactor)
            }
            mergeBranch2 <- mergeBranch2 / pre_cond # nolint
        }

    } else {
        if (cond == "maddison_cond") {
            mergeBranch2 <-
                mergeBranch2 / sum(weight_states * lambdas *
                                       (1 - nodeM[1:d]) ^ 2)
        }

        if (cond == "proper_cond") {
            mergeBranch2 <- mergeBranch2 / (lambdas * (1 - nodeM[1:d]) ^ 2)
        }
    }
    return(mergeBranch2)
}

#' @keywords internal
update_complete_tree <- function(phy,
                                 lambdas,
                                 mus,
                                 q_matrix,
                                 method,
                                 atol,
                                 rtol,
                                 lmb,
                                 use_normalization) {
    time_inte <- max(abs(ape::branching.times(phy))) # nolint

    if (is.list(lambdas)) {
        y <- rep(0, lmb)
        nodeM <- ct_condition_cpp(rhs = "ode_cla",
                                  y, # nolint
                                  time_inte,
                                  lambdas,
                                  mus,
                                  q_matrix,
                                  method,
                                  atol,
                                  rtol,
                                  use_normalization)
        nodeM <- c(nodeM, y) # nolint
    } else {
        y <- rep(0, 2 * lmb)
        nodeM <- ct_condition_cpp(rhs = "ode_standard",
                                  y, # nolint
                                  time_inte,
                                  lambdas,
                                  mus,
                                  q_matrix,
                                  method,
                                  atol,
                                  rtol,
                                  use_normalization)
    }
    return(nodeM)
}


#' @keywords internal
create_states <- function(usetraits,
                          traits,
                          states,
                          sampling_fraction,
                          num_concealed_states,
                          d,
                          traitStates,
                          is_complete_tree,
                          phy,
                          ly,
                          mus,
                          nb_tip) {
    if (anyNA(usetraits)) {
        nas <- which(is.na(traits))
        for (iii in seq_along(nas)) {
            states[nas[iii], ] <- c(1 - rep(sampling_fraction,
                                            num_concealed_states),
                                rep(sampling_fraction, num_concealed_states))
        }
    }

    for (iii in seq_along(traitStates)) {  # Initial state probabilities
        StatesPresents <- d + iii
        toPlaceOnes <- StatesPresents + 
                       length(traitStates) * (0:(num_concealed_states - 1))
        tipSampling <- 1 * sampling_fraction
        states[which(usetraits == traitStates[iii]), 
               toPlaceOnes] <- tipSampling[iii]
    }

  if (is_complete_tree) {
    extinct_species <- geiger::is.extinct(phy)
    if (!is.null(extinct_species)) {
      for (i in seq_along(extinct_species)) {
        
        entry <- mus * states[which(phy$tip.label == extinct_species[i]), 
                              (d + 1):ly]
        
        states[which(phy$tip.label == extinct_species[i]), (d + 1):ly] <-
          mus * states[which(phy$tip.label == extinct_species[i]), (d + 1):ly]
      }
    }
    for (iii in 1:nb_tip) {
      states[iii, 1:d] <- 0
    }
  } else {
    for (iii in 1:nb_tip) {
      states[iii, 1:d] <- rep(1 - sampling_fraction, num_concealed_states)
    }
  }
  return(states)
}

#' @keywords internal
build_states <- function(phy,
                         traits,
                         num_concealed_states,
                         sampling_fraction,
                         is_complete_tree = FALSE,
                         mus = NULL,
                         num_unique_traits = NULL,
                         first_time = FALSE,
                         traitStates = NULL) {
    if (!is.matrix(traits)) {
        traits <- matrix(traits, nrow = length(traits), ncol = 1, byrow = FALSE)
    }

    if (length(phy$tip.label) != nrow(traits)) {
     stop("Number of species in the tree must be the same as in the trait file")
    }
  
    # if there are traits that are not in the observed tree,
    # the user passes these themselves.
    # yes, this is a weird use-case

    if (is.null(traitStates)) traitStates <- sort(unique(traits[, 1]))

    if (!is.null(num_unique_traits)) {
        if (num_unique_traits > length(traitStates)) {
            if (first_time)
                message("found un-observed traits, expanding state space")
            traitStates <- 1:num_unique_traits
        }
    }
    
    obs_traits <- unique(traits[, 1])
    obs_traits <- obs_traits[!is.na(obs_traits)]
    if (sum(obs_traits %in% traitStates) != length(obs_traits)) {
      stop("Tip traits are not in idparslist")
    }

    nb_tip <- ape::Ntip(phy)
    nb_node <- phy$Nnode
    ly <- length(traitStates) * 2 * num_concealed_states
    states <- matrix(ncol = ly, nrow = nb_tip + nb_node)
    d <- ly / 2
    ## In a example of 3 states, the names of the colums would be like:
    ##
    ## colnames(states) <- c("E0A","E1A","E2A","E0B","E1B","E2B",
    ##                   "D0A","D1A","D2A","D0B","D1B","D2B")
    states[1:nb_tip, ] <- 0
    ## I repeat the process of state assignment as many times as columns I have
    for (iv in seq_len(ncol(traits))) {
        states <- create_states(traits[, iv],
                                traits,
                                states,
                                sampling_fraction,
                                num_concealed_states,
                                d,
                                traitStates,
                                is_complete_tree,
                                phy,
                                ly,
                                mus,
                                nb_tip)
    }

    return(states)
}

#' @keywords internal
build_initStates_time <- function(phy,
                                  traits,
                                  num_concealed_states,
                                  sampling_fraction,
                                  is_complete_tree = FALSE,
                                  mus = NULL,
                                  num_unique_traits = NULL,
                                  first_time = FALSE,
                                  traitStates = NULL) {
    states <- build_states(phy,
                           traits,
                           num_concealed_states,
                           sampling_fraction,
                           is_complete_tree,
                           mus,
                           num_unique_traits,
                           first_time,
                           traitStates)
    phy$node.label <- NULL
    split_times <- sort(event_times(phy), decreasing = FALSE)
    ances <- as.numeric(names(split_times))

    forTime <- cbind(phy$edge, phy$edge.length)

    return(list(
        states = states,
        ances = ances,
        forTime = forTime
    ))
}

#' @keywords internal
get_weight_states <- function(root_state_weight,
                              num_concealed_states,
                              mergeBranch,
                              lambdas,
                              nodeM,
                              d,
                              is_cla = FALSE) {

    if (is.numeric(root_state_weight)) {
        weight_states <- rep(root_state_weight / num_concealed_states,
                             num_concealed_states)
    } else {
        if (root_state_weight == "maddison_weights") {
            weight_states <- (mergeBranch) / sum((mergeBranch))
        }

        if (root_state_weight == "proper_weights") {
            if (is_cla) {
                lmb <- length(mergeBranch)
                numerator <- rep(NA, lmb)
                for (j in 1:lmb) {
                    numerator[j] <-
                        mergeBranch[j] / sum(lambdas[[j]] *
                                    ((1 - nodeM[1:d]) %o% (1 - nodeM[1:d])))
                }
                weight_states <- numerator / sum(numerator) # nolint
            } else {
                weight_states <- (mergeBranch /
                                      (lambdas * (1 - nodeM[1:d]) ^ 2)) /
                    sum((mergeBranch / (lambdas * (1 - nodeM[1:d]) ^ 2)))
            }
        }

        if (root_state_weight == "equal_weights") {
            weight_states <- rep(1 / length(mergeBranch), length(mergeBranch))
        }
    }

    return(weight_states)
}

#' Times at which speciation or extinction occurs
#' @title Event times of a (possibly non-ultrametric) phylogenetic tree
#' @param phy phylogenetic tree of class phylo, without polytomies, rooted and
#' with branch lengths. Need not be ultrametric.
#' @return times at which speciation or extinction happens.
#' @note This script has been modified from BAMMtools' internal function
#' NU.branching.times
#' @export
event_times <- function(phy) {
    if (ape::is.ultrametric(phy)) {
        return(ape::branching.times(phy))
    } else {
        if (ape::is.binary(phy) == FALSE) {
            stop("error. Need fully bifurcating (resolved) tree\n")
        }
        phy$begin <- rep(0, nrow(phy$edge))
        phy$end <- rep(0, nrow(phy$edge))
        fx <- function(phy, node) {
            cur_time <- 0
            root <- length(phy$tip.label) + 1
            if (node > root) {
                cur_time <- phy$end[which(phy$edge[, 2] == node)]
            }
            dset <- phy$edge[, 2][phy$edge[, 1] == node]
            i1 <- which(phy$edge[, 2] == dset[1])
            i2 <- which(phy$edge[, 2] == dset[2])
            phy$end[i1] <- cur_time + phy$edge.length[i1]
            phy$end[i2] <- cur_time + phy$edge.length[i2]
            if (dset[1] > length(phy$tip.label)) {
                phy$begin[phy$edge[, 1] == dset[1]] <- phy$end[i1]
                phy <- fx(phy, node = dset[1])
            }
            if (dset[2] > length(phy$tip.label)) {
                phy$begin[phy$edge[, 1] == dset[2]] <- phy$end[i2]
                phy <- fx(phy, node = dset[2])
            }
            return(phy)
        }
        phy <- fx(phy, node = length(phy$tip.label) + 1)
        maxbt <- max(phy$end)
        nodes <- (length(phy$tip.label) + 1):(2 * length(phy$tip.label) - 1)
        bt <- numeric(length(nodes))
        names(bt) <- nodes
        for (i in seq_along(bt)) {
            tt <- phy$begin[phy$edge[, 1] == nodes[i]][1]
            bt[i] <- maxbt - tt
        }
        return(bt)
    }
}

#' Print likelihood for initial parameters
#'
#' @inheritParams default_params_doc
#'
#' @return Invisible `NULL`. Prints a `message()` to the console with the
#'   initial loglikelihood if `verbose >= 1`
#' @noRd
print_init_ll <- function(initloglik) {
    init_ll_msg1 <- "Calculating the likelihood for the initial parameters."
    init_ll_msg2 <-
        paste0("The loglikelihood for the initial parameter values is ",
               initloglik)
    init_ll_msg3 <- c("Optimizing the likelihood - this may take a while.")
    message(paste(init_ll_msg1, init_ll_msg2, init_ll_msg3, sep = "\n"))
    
    invisible(NULL)
}

#' @keywords internal
set_and_check_structure_func <- function(idparsfuncdefpar,
                                         functions_defining_params,
                                         idparslist,
                                         idparsopt,
                                         idfactorsopt,
                                         idparsfix,
                                         initfactors) {
  structure_func <- list()
  structure_func[[1]] <- idparsfuncdefpar
  structure_func[[2]] <- functions_defining_params
  
  # checks specific to when the user has specified factors:
  
  if (is.null(idfactorsopt) == FALSE) {
    if (length(initfactors) != length(idfactorsopt)) {
      stop("idfactorsopt should have the same length as initfactors.")
    }
  }
  
  if (is.list(functions_defining_params) == FALSE) {
    stop(
      "The argument functions_defining_params should be a list of 
            functions. See example and vignette"
    )
  }
  
  if (length(functions_defining_params) != length(idparsfuncdefpar)) {
    stop(
      "The argument functions_defining_params should have the same 
            length than idparsfuncdefpar"
    )
  }
  
  if (anyDuplicated(c(idparsopt, idparsfix, idparsfuncdefpar)) != 0) {
    stop("At least one element was asked to be fixed, 
             estimated or a function at the same time")
  }
  
  if (identical(as.numeric(sort(
    c(idparsopt, idparsfix, idparsfuncdefpar)
  )), as.numeric(sort(unique(
    unlist(idparslist)
  )))) == FALSE) {
    stop(
      "All elements in idparslist must be included in either 
            idparsopt or idparsfix or idparsfuncdefpar "
    )
  }
  if (is.null(idfactorsopt)) {
    structure_func[[3]] <- "noFactor"
  } else {
    structure_func[[3]] <- idfactorsopt
  }
  
  return(structure_func)
}

#' @keywords internal
check_ml_conditions <- function(traits,
                                idparslist,
                                initparsopt,
                                idparsopt,
                                idparsfix,
                                parsfix) {
  if (is.matrix(traits)) {
    warning("you are setting a model where some species have more
            than one trait state")
  }
  
  if (length(initparsopt) != length(idparsopt)) {
    stop("initparsopt must be the same length as idparsopt.
             Number of parameters to optimize does not match the number of
             initial values for the search")
  }
  
  if (length(idparsfix) != length(parsfix)) {
    stop("idparsfix and parsfix must be the same length.
             Number of fixed elements does not match the fixed figures")
  }
  
  if (anyDuplicated(c(idparsopt, idparsfix)) != 0) {
    stop("at least one element was asked to be both fixed and estimated ")
  }
  
  if (anyDuplicated(c(unique(sort(as.vector(idparslist[[3]]))),
                      idparsfix[which(parsfix == 0)])) != 0) {
    warning("Note: you set some transitions as impossible to happen.")
  }
  
  if (min(initparsopt) <= 0.0) {
    stop("All elements in init_parsopt need to be larger than 0")
  }
}

#' @keywords internal
get_trait_states <- function(idparslist,
                             num_concealed_states,
                             display_warning = TRUE) {
  trait_names <- names(idparslist[[1]])
  if (is.null(trait_names)) trait_names <- names(idparslist[[2]])
  if (is.null(trait_names)) trait_names <- colnames(idparslist[[3]])
  
  
  if (is.null(trait_names)) return(NULL)

  # by convention, hidden states are appended A,B,C letters
  num_traits <- length(idparslist[[2]]) / num_concealed_states
  focal_names <- trait_names[1:num_traits]
  
  for (i in seq_along(focal_names)) {
    focal_names[i] <- substr(focal_names[i], 1, nchar(focal_names[i]) - 1)
  }
  
  if (display_warning) {
    output <- "Deduced names and order of used states to be: "
    
    for (i in seq_along(focal_names)) {
      if (i != length(focal_names)) {
        output <- paste0(output, focal_names[i], ", ")
      } else {
        output <- paste0(output, focal_names[i])
      }
    }
    
    rlang::warn(message = paste0(output, "\n", "if this is incorrect, consider passing states as matching numeric 
  ordering, e.g. 1 for the first state, 2 for the second etc."))
  }

  return(focal_names)
}

#' @keywords internal
extract_data <- function(mat) {
  to_plot <- c()
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      entry <- mat[i, j]
      if (!is.na(entry)) {
        if (entry > 0) {
          x <- colnames(mat)[i]
          y <- rownames(mat)[j]
          to_plot <- rbind(to_plot, c(x, y, entry))
        }
      }
    }
  }
  return(to_plot)
}

#' @keywords internal
get_rates <- function(lambda_list, all_states, lambda_order) {
  to_plot <- c()
  for (i in seq_along(all_states)) {
    local_v <- extract_data(lambda_list[[i]])
    colnames(local_v) <- c("x", "y", "rate")
    local_v <- as.data.frame(local_v)
    
    local_v$focal_rate <- lambda_order[i]
    to_plot <- rbind(to_plot, local_v)
  }
  
  to_plot$x <- factor(to_plot$x, levels = all_states)
  to_plot$y <- factor(to_plot$y, levels = rev(all_states))
  to_plot$focal_rate <- factor(to_plot$focal_rate, levels = all_states)
  return(to_plot)
}


#' function to visualize the structure of the idparslist
#' @param idparslist idparslist setup list
#' @param state_names names of all states (including concealed states)
#' @return list with two ggplot objects: plot_qmat which visualizes the 
#' q_matrix, and plot_lambda which visualizes the lambda matrices. 
#' @export
#' @examples
#' idparslist <- list()
#' focal_matrix <-
#' secsse::create_default_lambda_transition_matrix(state_names = c("1", "2"),
#'                                                 model = "CR")
#' idparslist[[1]] <- 
#'   secsse::create_lambda_list(state_names = c("1", "2"),
#'                              num_concealed_states = 2,
#'                              transition_matrix = focal_matrix,
#'                              model = "CR")
#' idparslist[[2]] <- secsse::create_mu_vector(state_names = c("1", "2"),
#'                                             num_concealed_states = 2,
#'                                             model = "CR",
#'                                             lambda_list = idparslist[[1]])
#' shift_mat <- secsse::create_default_shift_matrix(state_names = c("1", "2"),
#'                                                  num_concealed_states = 2,
#'                                                  mu_vector = idparslist[[2]])
#' idparslist[[3]] <- secsse::create_q_matrix(state_names = c("1", "2"),
#'                                            num_concealed_states = 2,
#'                                            shift_matrix = shift_mat,
#'                                            diff.conceal = FALSE)
#' p <- plot_idparslist(idparslist, state_names = names(idparslist[[1]]))
#' p$plot_lambda
#' p$plot_qmat
plot_idparslist <- function(idparslist, 
                            state_names) {
  
  v <- extract_data(idparslist[[3]])
  colnames(v) <- c("x", "y", "rate")
  v <- as.data.frame(v)
  
  v$x <- factor(v$x, levels = state_names)
  v$y <- factor(v$y, levels = rev(state_names))
  
  plot_qmat <- 
    ggplot2::ggplot(v, 
                    ggplot2::aes(x = .data[["x"]], 
                                 y = .data[["y"]], 
                                 fill = .data[["rate"]])) +
    ggplot2::geom_tile()
  
  to_plot <- get_rates(idparslist[[1]], state_names, state_names)
  
  to_plot$x <- factor(to_plot$x, levels = state_names)
  to_plot$y <- factor(to_plot$y, levels = rev(state_names))
  to_plot$focal_rate <- factor(to_plot$focal_rate, levels = state_names)
  
  plot_lambda <- 
    ggplot2::ggplot(to_plot, 
                    ggplot2::aes(x = .data[["x"]], 
                                 y = .data[["y"]], 
                                 fill = .data[["rate"]])) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~.data[["focal_rate"]])
  
  return(list("plot_qmat" = plot_qmat,
              "plot_lambda" = plot_lambda))
}