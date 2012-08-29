##*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
## Filename   : iBBiG_makeSimulatedData.R
## Author     : Aedin Culhane
## eMail      : aedin@jimmy.harvard.edu
## Date       : Nov 22nd 2011
## Description: Functions to create Simulated data and run iBBiG, Gusenleitner et al., submitted to Bioinformatics, 2011
##*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*



addSignal<-function(arti,startC,startR,endC,endR,densityH,densityL){
    # simple function to add noise to binary matrix
	numRows<-endR-startR+1
	numCols<-endC-startC+1
	diffDens<-densityH-densityL
	stepD<-diffDens/numRows
	for (i in startR:endR){
	  random<-as.integer(runif(numCols)<(densityH-(i-startR)*stepD))
	  arti[i,startC:endC]<-random
	  }  
	return(arti)
	}


makeSimDesignMat<-function(verbose=TRUE) {
   ## Simple function to create a design matrix for the simulated data used by GusenLeitner et al.,
   
   designMat<-rbind(c(251,51,275,300,0.9,0.4),
                    c(51,251,225,325,0.8,0.4),
					c(1,1,50,50,0.8,0.5),
                    c(46,46,85,85,0.9,0.4), 
                    c(81,81,110,110,0.8,0.4),
                    c(106,106,125,125,0.9,0.6),
                    c(151,151,190,190,0.6,0.5)
                    )
   colnames(designMat)<-c("startC","startR","endC","endR","densityH","densityL") 
   rownames(designMat)<-paste("M", 1:7,sep="")

   if (verbose){
     print("***** Summary of Design Matrix ******")
     print(cbind(Rows=(designMat[,4]-designMat[,2])+1,Cols=(designMat[,3]-designMat[,1])+1,DensityLow=designMat[,6], DenistyHigh=designMat[,5]))
     }
  # print(designMat)
   return(designMat)
   }


makeArtificial<-function(nRow=400,nCol=400,noise=0.1, verbose=TRUE, dM=makeSimDesignMat(verbose=verbose), seed=123){
   ##creates an artificial dataset 
   ## Given a design Matrix(dM) which specifies the cluster size and density
   ## where the nrow = the number of clusters and the ncol = 6 with the columns
   ## headings "startC","startR","endC","endR","densityH","densityL". This function
   ## creates an matrix of simulated  (artifical,arti) data.
   require(biclust)

   if (!is.null(seed)) set.seed(seed)

   arti<-matrix(0,nRow,nCol) 

   arti[1:nRow,1:nCol]<-as.integer(runif(nRow*nCol)<noise)
   rownames(arti)<-paste('sig',1:nRow,sep='_')
   colnames(arti)<-paste('cov',1:nCol,sep='_')
   for (i in 1:nrow(dM)) {
     arti<-addSignal(arti,dM[i,1],dM[i,2],dM[i,3],dM[i,4],dM[i,5],dM[i,6])
     }

  ## Format the data in the iBBiG class that extends biclust

  nClust=nrow(dM)
  RClust<-apply(dM[,c(2,4)], 1,function(x)x[1]:x[2])
  CClust<-apply(dM[,c(1,3)], 1,function(x)x[1]:x[2])
  
  RN<-matrix(FALSE, nrow(arti), nClust)
  for (i in 1:nClust) RN[RClust[[i]],i] <-TRUE
  NC <-matrix(FALSE, nClust, ncol(arti))
  for (i in 1:nClust) NC[i,CClust[[i]]] <-TRUE
  if (verbose) {
  	cat("\nCluster sizes in new iBBiG (Biclust) data object \n")
  	cat(paste("Number of Modules:", nClust, "\n"))
  	cat(c("Rows",apply(RN,2,sum), "\n"), sep="\t")
  	cat(c("Columns", apply(NC,1, sum),"\n"), sep="\t")
        }
  
  new("iBBiG", Parameters=list(designMatrix=dM, nRow=nRow, nCol=nCol), Seeddata = arti, RowxNumber = RN, NumberxCol = NC, Number = nClust)

}




