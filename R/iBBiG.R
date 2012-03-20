# Author     : Daniel Gusenleitner
# eMail      : gusef@jimmy.harvard.edu
# Date       : 01-20-11
# Description: Genetic algorithm version of the bi-clustering algorithm
#*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
#############################################################
      

calculaterowScore<-function(colVector,binaryMatrix,alpha){
  ##calculates the rowScores for each term
  colSize<-sum(colVector)
  currentCluster<-binaryMatrix[,colVector==1]
  pScore<-rowSums(currentCluster)
  p1<-pScore/colSize
  p0<-1-p1
  entropy<--p1*log2(p1)-p0*log2(p0)
  entropy[is.na(entropy)]<-0
  eScore<-(1-entropy)
  pScore[p1<0.5]<-0
  score<-pScore*(eScore^alpha)
  return(score)
}

removeRowInformation<-function(rowScore,colVector,alpha){
  ##weights the pValues in the covariate matrix, that produced the score for the
  ##module. the weight is 1 minus the entropy part of the score
  colSize<-sum(colVector)
  currentCluster<-rowScore[colVector==1]
  pScore<-sum(currentCluster)
  p1<-pScore/colSize
  if (p1>=0.5){
     p0<-1-p1
     entropy<--p1*log2(p1)-p0*log2(p0)
     entropy[is.na(entropy)]<-0
     weight<-1-(1-entropy)^alpha
     rowScore[colVector==1]<-rowScore[colVector==1]*weight
   }
   return(rowScore)
}

removeInformation<-function(colVector,binaryMatrix,alpha){
   ##weights the used rows and thereby removes the information
   ##that was used in the top cluster

   binaryMatrix<-apply(binaryMatrix,1,removeRowInformation,colVector,alpha)
   binaryMatrix<-t(binaryMatrix)
   return(binaryMatrix)
}
    

linearRanking<-function(pos,nInd,SP){
   ##determine the selection probability (linear ranking) 2-SP+2*(SP-1)*(Pos-1)/(Nind-1)
   probab<-2-SP+2*(SP-1)*(pos-1)/(nInd-1)
   return(probab)
}       


getCumProbabilities<-function(probabilities){
   ##calculates the cummulated propabilites   
   probabilities<-probabilities/sum(probabilities)
   for (i in 2:length(probabilities)){
      probabilities[i]<-probabilities[i]+probabilities[i-1]
   }
   return(probabilities)
}       

#discretizerowScores<-function(rowScores){
     ## In order to get a binary vector of membership of rows in clusters need to discretize the score
    ## Will do a simple approach a fit a bimodal distribution using mclust.
 
#     require(mclust)
     #RN<-try(Mclust(rowScores, G=2)$classification-1)    
     ##work around to convert lists returned with try catch error
#       RN <- tryCatch(Mclust(rowScores, G = 2)$classification - 1, error = function(e) return((x)))
     #if (inherits(RN, "list")) RN<-do.call(rbind, lapply(RN, unlist))
#     return(RN)
#     }



iBBiG<-function(binaryMatrix, 
                nModules,     
                alpha=0.3,
                pop_size=100,
                mutation=0.08,
                stagnation=50,
                selection_pressure=1.2,
                max_sp=15,
                success_ratio=0.6){

   #save the raw data for later
   rawMatrix<-binaryMatrix
   
   selectionP<-sapply(1:pop_size,linearRanking,pop_size,selection_pressure)
   sp<-getCumProbabilities(selectionP)

   ##annotation for the results
   ##binary matrix where 1 represents an association between a row and column
   NumberxCol<-matrix(nrow=0, ncol=ncol(binaryMatrix))
   colnames(NumberxCol)<-colnames(binaryMatrix)
   

   rowScorexNumber<-matrix(nrow=nrow(binaryMatrix),ncol=0)
   rownames(rowScorexNumber)<-rownames(binaryMatrix)

   RowxNumber <-matrix(nrow= nrow(binaryMatrix), ncol=0)
   rownames(RowxNumber)<-rownames(binaryMatrix)


   clusterScores<-numeric(0)

   ## Number of rows (Sigs) and columns (Covariates)
   colSize<-ncol(binaryMatrix)
   rowSize<-nrow(binaryMatrix)

   ## nModules is the number of clusters to be discovered. Normally this is set to an excess number as the true number of clusters can be estimated
   ## from score and size of resulting clusters

   for (i in 1:nModules){
      cat('Module: ',i)
      colVector<-rep(0,colSize)
         
      
      out <- .C("clusterCovsC",
                as.double (binaryMatrix),
                as.integer(colVector),
                as.integer(colSize),
                as.integer(rowSize),
                as.double (alpha),
                as.integer(pop_size),
                as.integer(stagnation),
                as.double (mutation),
                as.double (success_ratio),
                as.integer(max_sp),
                as.double (sp),
                PACKAGE="iBBiG")
     
      colVector<-out[[2]] 

	  if (sum(colVector)==0){
           rowScore<-rep(FALSE,rowSize)
           rowVector<-rowScore
           clusterScore<-0
       }
	  else{
      #calculate row scores
      rowScore<-calculaterowScore(colVector,binaryMatrix,alpha)
      rowVector<-rowScore>0
      clusterScore<-sum(rowScore)
	    #remove information from binaryMatrix that was used in this cluster
      binaryMatrix<-removeInformation(colVector,binaryMatrix,alpha)    
	   }
	   
      
      #add metainformation (score,size) 
      clusterSize<-sum(colVector)
      clusterScores<-c(clusterScores,clusterScore)

      #add cluster to result matrices
      rowScorexNumber<-cbind(rowScorexNumber,rowScore)
      RowxNumber<-cbind(RowxNumber,rowVector)
      NumberxCol<-rbind(NumberxCol,colVector>0)
      cat(' ... done\n')
   }
    colnames(RowxNumber) <-paste("M", 1:ncol(RowxNumber))
    rownames(NumberxCol) <-paste("M", 1:nrow(NumberxCol))
    colnames(rowScorexNumber) <-paste("M", 1:ncol(rowScorexNumber))
    names(clusterScores)<-paste("M", 1:length(clusterScores))
 
  res<- new("iBBiG", Seeddata=rawMatrix, Clusterscores=clusterScores, RowScorexNumber=rowScorexNumber,  NumberxCol = NumberxCol, RowxNumber=RowxNumber, Number = nModules)
   return(res)
}

# library(iBBiG)
# simData<-makeArtificial(verbose=FALSE)
# binaryMatrix<-simData@Seeddata
# clusters<-iBBiG(binaryMatrix,8
