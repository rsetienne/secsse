#' It sets the parameters (speciation, extinction and transition)
#' ids. Needed for ML calculation (secsse_ml)
#' @title Parameter structure setting
#' @param traits vector with trait states, order of states must be the same as tree tips, for help, see vignette.
#' @param num_concealed_states number of concealed states, generally equivalent to number of examined states.
#' @return A list that includes the ids of the parameters for ML analysis. 
#' @examples
#'traits  <-  sample(c(0,1,2), 45,replace=TRUE) #get some traits
#'num_concealed_states  <-  3 
#'param_posit <- id_paramPos(traits,num_concealed_states)
#' @export

id_paramPos  <-  function(traits,num_concealed_states){
  idparslist <- list()
  if(is.matrix(traits)){
    traits <- traits[,1]
  }
  
  ly  <-  length(sort(unique(traits))) * 2 * num_concealed_states
  d  <-  ly/2
  idparslist[[1]]  <-  1:d
  idparslist[[2]]  <-  (d+1):ly
  toMatrix <- 1
  matPos <- (ly+1): (((d^2)-d) + d * 2)
  for (i in 1:d){
    toMatrix <- c(toMatrix,matPos[(i*d-(d-1)):((i*d-(d-1))+d)])
    
  }
  toMatrix <- toMatrix[1:d^2]
  Q  <-  matrix(toMatrix,ncol=d,nrow=d,byrow=T)
  diag(Q)  <-  NA
  idparslist[[3]]  <-  Q
  
  lab_states <- rep(as.character(sort(unique(traits))),num_concealed_states)
  
  lab_conceal <- NULL
  for (i in 1:num_concealed_states){
    
    lab_conceal  <-  c(lab_conceal,rep(LETTERS[i],length(sort(unique(traits)))))
  }
  
  
  statesCombiNames  <-  character()
  for ( i in 1:length(lab_states)){
    statesCombiNames  <-  c(statesCombiNames,paste0(lab_states[i],lab_conceal[i]))
    
  }
  colnames(idparslist[[3]]) <- statesCombiNames
  rownames(idparslist[[3]]) <- statesCombiNames
  names(idparslist) <- c("lambdas","mus","Q")
  names(idparslist[[1]]) <- statesCombiNames
  names(idparslist[[2]]) <- statesCombiNames
  return(idparslist)
}



#' Sets a Q matrix where double transitions are not allowed
#' @title Basic Qmatrix 
#' @param traits vector with trait states, order of states must be the same as tree tips, for help, see vignette.
#' @param masterBlock matrix of transitions among only examined states, NA in the main diagonal, used to build the full transition rates matrix. 
#' @param diff.conceal should the concealed states be different? Normally it should be FALSE. 
#' @return Q matrix that includes both examined and concealed states, it should be declared as the third element of idparslist.
#' @examples
#' traits  <-  sample(c(0,1,2), 45,replace=TRUE) #get some traits
#' masterBlock <- matrix(99,ncol=3,nrow=3,byrow=TRUE) #For a three-state trait
#' diag(masterBlock) <- NA
#' masterBlock[1,2] <- 6
#' masterBlock[1,3] <- 7
#' masterBlock[2,1] <- 8
#' masterBlock[2,3] <- 9
#' masterBlock[3,1] <- 10
#' masterBlock[3,2] <- 11
#' myQ <- q_doubletrans(traits,masterBlock,diff.conceal=FALSE)
#' # now, it can replace the Q matrix from id_paramPos  
#' num_concealed_states  <-  3 
#' param_posit <- id_paramPos(traits,num_concealed_states)
#' param_posit[[3]] <-  myQ
#' @export

q_doubletrans <- function(traits,masterBlock,diff.conceal){
  
  if(diff.conceal==TRUE && all(floor(masterBlock)==masterBlock,na.rm = TRUE)==FALSE){
    integersmasterBlock <- floor(masterBlock)
    factorBlock <- signif(masterBlock-integersmasterBlock,digits = 2)
    
    factorstoExpand <- unique(sort(c(factorBlock)))
    factorstoExpand <- factorstoExpand[factorstoExpand>0]
    newshareFac <- (max(factorstoExpand*10)+1):(max(factorstoExpand*10)+length(factorstoExpand))
    newshareFac <- newshareFac/10
    
    for(iii in 1:length(newshareFac)){  
      factorBlock[which(factorBlock==factorstoExpand[iii])] <- newshareFac[iii]  
      
    }
    
    ntraits <- length(sort(unique(traits)))
    uniqParQ <- sort(unique(c(floor(masterBlock))))
    uniqParQ2 <- uniqParQ[which(uniqParQ>0)]
    concealnewQ <- (max(uniqParQ2)+1):(max(uniqParQ2)+length(uniqParQ2))
    
    for(iii in 1:length(concealnewQ)){  
      integersmasterBlock[which(integersmasterBlock==uniqParQ2[iii])] <- concealnewQ[iii]  
      
    }
    concealnewQMatr <- integersmasterBlock+factorBlock
    
    Q <- NULL
    for(i in 1:ntraits){
      Qrow <- NULL
      for(ii in 1:ntraits){
        entry <- masterBlock[i,ii]
        if(is.na(entry)){
          Qrow <- cbind(Qrow,masterBlock)
        } else{
          entry <- concealnewQMatr[i,ii]
          
          outDiagBlock <- matrix(0,ncol=ntraits,nrow=ntraits,byrow=T)
          diag(outDiagBlock) <- entry
          Qrow <- cbind(Qrow,outDiagBlock)
        }
        
      }
      Q <- rbind(Q,Qrow)
    }
  } else {
    
    ntraits <- length(sort(unique(traits)))
    uniqParQ <- sort(unique(c(masterBlock)))
    uniqParQ2 <- uniqParQ[which(uniqParQ>0)]
    concealnewQ <- (max(uniqParQ2)+1):(max(uniqParQ2)+length(uniqParQ2))
    concealnewQMatr <- masterBlock
    for( I in 1: length(uniqParQ2)){
      uniqParQ2
      concealnewQMatr[concealnewQMatr==uniqParQ2[I]] <- concealnewQ[I]
    }
    
    Q <- NULL
    for(i in 1:ntraits){
      Qrow <- NULL
      for(ii in 1:ntraits){
        entry <- masterBlock[i,ii]
        if(is.na(entry)){
          
          Qrow <- cbind(Qrow,masterBlock)
        } else{
          
          if(diff.conceal==TRUE){
            entry <- concealnewQMatr[i,ii]
          }
          outDiagBlock <- matrix(0,ncol=ntraits,nrow=ntraits,byrow=T)
          diag(outDiagBlock) <- entry
          Qrow <- cbind(Qrow,outDiagBlock)
        }
        
      }
      Q <- rbind(Q,Qrow)
    }
  }
  return(Q)
}


#' In preparation for likelihood calculation, it orders trait data according the tree tips
#' @title Data checking and trait sorting
#' @param traitinfo data frame where first column has species ids and the second one is the trait associated information.
#' @param phy phy phylogenetic tree of class phylo, ultrametric, fully-resolved, rooted and with branch lengths.
#' @return Vector of traits
#' @examples
#' # Some data we have prepared
#' data(traitinfo)
#' data("phylo_Vign")
#' traits <- sortingtraits(traitinfo,phylo_Vign)
#' @export


sortingtraits <- function(traitinfo,phy){
  traitinfo <- as.matrix(traitinfo)
  if(length(phy$tip.label)!=nrow(traitinfo)){
    stop("Number of species in the tree must be the same as in the trait file")
  }
  
  if(identical(as.character(sort(phy$tip.label)),
               as.character(sort(traitinfo[,1])))==FALSE){
    mismatch <- match( as.character(sort(traitinfo[,1])),as.character(sort(phy$tip.label)))
    mismatched <- (sort(traitinfo[,1]))[which(is.na(mismatch))]
    stop(cat("Mismatch on tip labels and taxa names, check the species:",
             mismatched))
  }
  
  traitinfo <- traitinfo[match(phy$tip.label,traitinfo[,1]),]
  traitinfo[,1]==phy$tip.label
  
  if(ncol(traitinfo)==2){
    traits <- as.numeric(traitinfo[,2])
  }
  
  if(ncol(traitinfo)>2){
    traits <- NULL
    for (i in 1:(ncol(traitinfo)-1)){
      
      traits <- cbind(traits,as.numeric(traitinfo[,1+i]))
      
    }
    
  }
  return(traits)
}

#' It sets the parameters (speciation, extinction and transition)
#' ids. Needed for ML calculation with cladogenetic options (cla_secsse_ml)
#' @title Parameter structure setting for cla_secsse
#' @param traits vector with trait states, order of states must be the same as tree tips, for help, see vignette.
#' @param num_concealed_states number of concealed states, generally equivalent to number of examined states.
#' @return A list that includes the ids of the parameters for ML analysis. 
#' @examples
#'traits  <-  sample(c(0,1,2), 45,replace=TRUE) #get some traits
#'num_concealed_states  <-  3 
#'param_posit <- cla_id_paramPos(traits,num_concealed_states)
#' @export

cla_id_paramPos  <-  function(traits,num_concealed_states){
  idparslist <- list()
  if(is.matrix(traits)){
    traits <- traits[,1]
  }
  
  ly  <-  length(sort(unique(traits))) * 2 * num_concealed_states
  d  <-  ly/2
  #idparslist[[1]]  <-  1:d
  toMatrix <- 1
  matPos <- (ly+1): (((d^2)-d) + d * 2)
  for (i in 1:d){
    toMatrix <- c(toMatrix,matPos[(i*d-(d-1)):((i*d-(d-1))+d)])
    
  }
  toMatrix <- toMatrix[1:d^2]
  Q  <-  matrix(toMatrix,ncol=d,nrow=d,byrow=T)
  diag(Q)  <-  NA
  lab_states <- rep(as.character(sort(unique(traits))),num_concealed_states)
  
  lab_conceal <- NULL
  for (i in 1:num_concealed_states){
    
    lab_conceal  <-  c(lab_conceal,rep(LETTERS[i],length(sort(unique(traits)))))
  }
  
  
  statesCombiNames  <-  character()
  for ( i in 1:length(lab_states)){
    statesCombiNames  <-  c(statesCombiNames,paste0(lab_states[i],lab_conceal[i]))
    
  }
  
  

  idparslist[[1]]<- matrix(0,ncol=d,nrow=4)
  idparslist[[2]]  <-  (d+1):ly
  idparslist[[3]]  <-  Q

  rownames(idparslist[[1]])<-c("dual_inheritance",
                      "single_inheritance",
                      "dual_symmetric_transition",
                      "dual_asymmetric_transition")
  
  colnames(idparslist[[1]])<-statesCombiNames
  colnames(idparslist[[3]]) <- statesCombiNames
  rownames(idparslist[[3]]) <- statesCombiNames
  names(idparslist) <- c("lambdas","mus","Q")
  names(idparslist[[2]]) <- statesCombiNames
  return(idparslist)
}


#' It provides the set of matrices containing all the speciation rates
#' @title Prepares the entire set of lambda matrices for cla_secsse.
#' @param traits vector with trait states, order of states must be the same as tree tips, for help, see vignette.
#' @param num_concealed_states number of concealed states, generally equivalent to number of examined states.
#' @param lambd_and_modeSpe a matrix with the 4 models of speciation possible.
#' @return A list of lambdas, its lenght would be the same than the number of trait states * num_concealed_states.. 
#' @export

prepare_full_lambdas<-function(traits,num_concealed_states,lambd_and_modeSpe){
num_exami<-length(sort(unique(traits)))
  mat_size<-num_exami*num_concealed_states
posib_trans<-matrix(1,ncol=num_exami,nrow=num_exami,byrow=TRUE)
diag(posib_trans)<-NA
posib_trans<-q_doubletrans(traits,masterBlock=posib_trans,diff.conceal=FALSE)

full_lambdas<-list()

for(jj in 1:mat_size){
  #dual_state_inhe
  m1<-matrix(0,ncol=mat_size,nrow=mat_size)
  m1[jj,jj] <- as.numeric(lambd_and_modeSpe[,jj][1])
  
  #single_state_inhe
  m2<-matrix(0,ncol=mat_size,nrow=mat_size)
  m2[,jj]<-posib_trans[jj,]
  m2[jj,jj] <- 0
  m2[m2==1] <-as.numeric(lambd_and_modeSpe[,jj][2])
  #symet_state_emerge
  
  m3<-matrix(0,ncol=mat_size,nrow=mat_size)
  
  diag(m3)<-posib_trans[jj,]
  m3[jj,jj] <- 0
  m3[m3==1] <-as.numeric(lambd_and_modeSpe[,jj][3])
  #symet_state_emerge
  
  m4<-matrix(0,ncol=mat_size,nrow=mat_size)
  for(i in 1:length(which(posib_trans[jj,]==1))){
    m4[which(posib_trans[jj,]==1)[i],]<-posib_trans[jj,] 
  }
  m4[,jj]<-0
  m4[upper.tri(m4)] <- 0
  diag(m4)<-0
  m4[is.na(m4)] <-0
  m4[m4==1] <-as.numeric(lambd_and_modeSpe[,jj][4])
  
  full_lambdas[[jj]]<-m1+m2+m3+m4
}
return(full_lambdas)
}

