##*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
## Filename   : iBBiG_Class.R
## Author     : Aedin Culhane
## eMail      : aedin@jimmy.harvard.edu
## Date       : Nov 22nd 2011
## Description: Define the class and methods for iBBiG
##*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*



setClass("iBBiG", representation(
               Seeddata="matrix", 
               RowScorexNumber="matrix", 
               Clusterscores ="numeric"), 
               contains="Biclust")

#RowScorexNumber accessor
setGeneric("RowScorexNumber", function(object) standardGeneric("RowScorexNumber"))
setGeneric("RowScorexNumber<-", function(object, value) standardGeneric("RowScorexNumber<-"))
setMethod(RowScorexNumber, "iBBiG", function(object) slot(object, "RowScorexNumber"))
setReplaceMethod("RowScorexNumber", "iBBiG", function(object, value)
{
  slot(object, "RowScorexNumber") <- value
  validObject(object)
  object
})

#Clusterscores accessor
setGeneric("Clusterscores", function(object) standardGeneric("Clusterscores"))
setGeneric("Clusterscores<-", function(object, value) standardGeneric("Clusterscores<-"))
setMethod(Clusterscores, "iBBiG", function(object) slot(object, "Clusterscores"))
setReplaceMethod("Clusterscores", "iBBiG", function(object, value)
{
  slot(object, "Clusterscores") <- value
  validObject(object)
  object
})

#Seeddata accessor
setGeneric("Seeddata", function(object) standardGeneric("Seeddata"))
setGeneric("Seeddata<-", function(object, value) standardGeneric("Seeddata<-"))
setMethod(Seeddata, "iBBiG", function(object) slot(object, "Seeddata"))
setReplaceMethod("Seeddata", "iBBiG", function(object, value)
{
  slot(object, "Seeddata") <- value
  validObject(object)
  object
})

#Parameters accessor
setGeneric("Parameters", function(object) standardGeneric("Parameters"))
setGeneric("Parameters<-", function(object, value) standardGeneric("Parameters<-"))
setMethod(Parameters, "iBBiG", function(object) slot(object, "Parameters"))
setReplaceMethod("Parameters", "iBBiG", function(object, value)
{
  slot(object, "Parameters") <- value
  validObject(object)
  object
})

#RowxNumber accessor
setGeneric("RowxNumber", function(object) standardGeneric("RowxNumber"))
setGeneric("RowxNumber<-", function(object, value) standardGeneric("RowxNumber<-"))
setMethod(RowxNumber, "iBBiG", function(object) slot(object, "RowxNumber"))
setReplaceMethod("RowxNumber", "iBBiG", function(object, value)
{
  slot(object, "RowxNumber") <- value
  validObject(object)
  object
})

#NumberxCol accessor
setGeneric("NumberxCol", function(object) standardGeneric("NumberxCol"))
setGeneric("NumberxCol<-", function(object, value) standardGeneric("NumberxCol<-"))
setMethod(NumberxCol, "iBBiG", function(object) slot(object, "NumberxCol"))
setReplaceMethod("NumberxCol", "iBBiG", function(object, value)
{
  slot(object, "NumberxCol") <- value
  validObject(object)
  object
})

#Number accessor
setGeneric("Number", function(object) standardGeneric("Number"))
setGeneric("Number<-", function(object, value) standardGeneric("Number<-"))
setMethod(Number, "iBBiG", function(object) slot(object, "Number"))
setReplaceMethod("Number", "iBBiG", function(object, value)
{
  slot(object, "Number") <- value
  validObject(object)
  object
})

#info accessor
setGeneric("info", function(object) standardGeneric("info"))
setGeneric("info<-", function(object, value) standardGeneric("info<-"))
setMethod(info, "iBBiG", function(object) slot(object, "info"))
setReplaceMethod("info", "iBBiG", function(object, value)
{
  slot(object, "info") <- value
  validObject(object)
  object
})


### PLOT *******************************

setMethod("plot", "iBBiG", 
   function(x, type="color", n=1:x@Number, col="default", reorder=F, main=NULL, ...)  {

     if(length(x@Clusterscores)==0){
       layout(matrix(c(1,1,1,1,2,2), 3, 2, byrow=TRUE), respect=TRUE)
     }else{
       layout(matrix(c(1,1,1,1,1,1,2,3,4), 3, 3, byrow=TRUE))
     }

    if (nrow(x@Seeddata)>0){
      iBBiGmat<-x@Seeddata
    }else{
      iBBiGmat<-matrix(0, ncol=ncol(x@NumberxCol), nrow=nrow(x@RowxNumber))
    } 

    ## white = #FFFFFF (no association), gray = #666666 (noise)
    colBW= c("#FFFFFF","#666666")

    if (type=="bw") { col= colBW; legendcol=colBW[2]}

    if (type=="color") {
        if ("default"%in%col) col<-colors()[c(367,565,506,645,639,657,373,258,619,573,631,543)]
        #col<-c('#FF0000','#D0006E','#25D500','#B48A00','#FF9900','#1049A9','#33CDC7')
        #col=c('#FF0000','#FF69B4','#B48A00','#7fff00','#1049A9','#FF9900','#33CDC7','#FFFF00', "#6600FF" ,'#D0006E','#eecbad')
          
          
        while (length(col)<length(n)) col=c(col,col)
        legendcol=col[n]
        col= c(colBW, legendcol)

        for (i in 1:length(n)){
	        rInd<-x@RowxNumber[,n[i]]==1
            cInd<-x@NumberxCol[n[i],]==1
            iBBiGmat[rInd,cInd][iBBiGmat[rInd,cInd]>0] <-(i+1)
	      }
     }
     
     if (reorder){
        #recursive function to reorder the seeding matrix
        orderRec<-function(col,module){
           if (col<ncol(module)){
             #all covariates that belong to this module
             currentCovs<-module[,col]==1
             sub1<-sum(currentCovs)>0
             if (sub1){
               subMod<-module
               subMod[!currentCovs,]<-2
               covOrder<-orderRec(col+1,subMod)                             
             }
         
             rest<-module[,col]==0 
             sub2<-sum(rest)>0
             if (sub2){
                subMod<-module
                subMod[!rest,]<-2
                restOrder<-orderRec(col+1,subMod)
             }
         
              if (sub1 && sub2){
                 combOrder<-c(covOrder,restOrder)
              }else if (sub1){
                 combOrder<-covOrder
              }else{
                 combOrder<-restOrder
              }
           }else{
              combOrder<-(1:nrow(module))[module[,1]!=2]
           }
           return(combOrder)
        }
     
        covOrder<-orderRec(1,t(x@NumberxCol))
        sigOrder<-orderRec(1,x@RowxNumber)
        iBBiGmat<-iBBiGmat[sigOrder,covOrder]
     }

    iBBiGmat<-t(iBBiGmat)[,ncol(iBBiGmat):1]
    image(iBBiGmat,col=col,xlab='Phenotypes',ylab='Gene Signatures',axes=F,...)
    if (length(legendcol)>1) legend("topright", legend=as.character(paste("M", n, sep="")), fill=legendcol, col=legendcol)
    if (!is.null(title)) title(main = main, font.main = 4, ...)

    barplot(apply(x@NumberxCol, 1, sum), col= legendcol, main="Module Size", ylab="Number (Phenotypes)", las=2)
	
    #barplot(apply(x@RowxNumber, 3, sum), col= legendcol, main="Module Size", ylab="Number of Gene Sets", las=2)
     
    if(length(x@Clusterscores)>0) { 
	   barplot(x@Clusterscores, col= legendcol, main="Module Score", ylab="Pairwise Test Score", las=2)
	   clusterWeight <-log(apply(x@NumberxCol,1,sum)/ncol(x@NumberxCol)*(x@Clusterscores))
	   clusterWeight[is.na(clusterWeight)]<-0
       barplot(clusterWeight, col= legendcol, main="Weighted Score", ylab="Weighted Score", las=2)
      } 
     })


## as Convert BiClust to iBBiG

setAs("Biclust" , "iBBiG",
  function (from , to ){
         new(to, Seeddata = from@Parameters$seeddata, RowScorexNumber=matrix(),Clusterscores =numeric(), RowxNumber = from@RowxNumber, NumberxCol = from@NumberxCol, Number = from@Number)
       })
### show and summary ******************************


setMethod("show", "iBBiG",
  function(object){
  cat("\nAn object of class",class(object),"\n")

    n<-object@Number
    n<-min(c(n,5))
    if(n>1)
    {
    cat("\nNumber of Clusters found: ",object@Number, "\n")
    cat("\nFirst ",n," Cluster scores and sizes:\n")

    rowcolsizes<-rbind(object@Clusterscores[1:n], colSums(object@RowxNumber[,1:n]),rowSums(object@NumberxCol[1:n,]))
    rownames(rowcolsizes)<-c("Cluster Score", "Number of Rows:","Number of Columns:")
    colnames(rowcolsizes)<-colnames(object@RowxNumber)[1:n]
    #print.default(format(rowcolsizes, print.gap = 2, quote = FALSE))
    print(rowcolsizes)
    }
    else
    {
    if(n==1) cat("\nThere was one cluster found with Score", object@Clusterscores[1],"and \n ",sum(object@RowxNumber[,1]), "Rows and ", sum(object@NumberxCol), "columns")
    if(n==0) cat("\nThere was no cluster found")
    }

    cat("\n\n")
})


setMethod("summary", "iBBiG",
  function(object)
  {
  cat("\nAn object of class",class(object),"\n")

    n<-object@Number
    n<-min(c(n,5))
    if(n>1)
    {
    cat("\nNumber of Clusters found: ",object@Number, "\n")
    cat("\nFirst ",n," Cluster scores and sizes:\n")

    rowcolsizes<-rbind(object@Clusterscores[1:n], colSums(object@RowxNumber[,1:n]),rowSums(object@NumberxCol[1:n,]))
    rownames(rowcolsizes)<-c("Cluster Score", "Number of Rows:","Number of Columns:")
    colnames(rowcolsizes)<-paste("M", 1:n)
    #print.default(format(rowcolsizes, print.gap = 2, quote = FALSE))
    print(rowcolsizes)
    }
    else
    {
    if(n==1) cat("\nThere was one cluster found with Score", object@Clusterscores[1],"and \n ",sum(object@RowxNumber[,1]), "Rows and ", sum(object@NumberxCol), "columns")
    if(n==0) cat("\nThere was no cluster found")
    }
})


#---------------------------
# Reorder/Subset
setMethod("[", 'iBBiG',
          function(x, i, drop=FALSE) {
            ind<-1:x@Number         
            if (!all(abs(i)%in%ind)) stop(paste("Invalid select"))
            x@info<-list("originalOrder"=ind, "select"=i)         
            x@RowxNumber<-x@RowxNumber[,i, drop=drop]
            x@NumberxCol<-x@NumberxCol[i, , drop=drop]
            x@RowScorexNumber<-x@RowScorexNumber[,i, drop=drop]
            x@Clusterscores<-x@Clusterscores[i]
            x@Number<-length(x@Clusterscores)                   
            #   print("cluster order updated, information about order is available in slot info")
            return(x)
          })

#setClassUnion("iBBiGorBiclust", c("iBBiG","Biclust")) 

setGeneric("JIdist", function(clustObj, GS,...) standardGeneric("JIdist"))
#setMethod("JIdist", c('iBBiGorBiclust','iBBiGorBiclust'), 
setMethod("JIdist", c('Biclust','Biclust'),
function(clustObj, GS, margin="col", best=TRUE) {
  ## Calculates Jaccard Index between, each Module and modules discovered
  ## if best is TRUE, it will report the highest JI for each module
  ## Both clustObj and GS have class biclust. 
  ## clustObj is a result of a cluster analysis
  ## GS is an simulated dataset in biclust format 
  ## GS is generated tt<-makeArtifical() 
  ## dist.binary report distance sqrt(1-S), so convert to JI
  ## margin is whether to calculate the JI distance on columns or rows. 
  
  #require(ade4)
  nC= GS@Number
  if (clustObj@Number==0) return(NULL)
  GSnames<-colnames(GS@RowxNumber)
  if (is.null(GSnames)) GSnames<- paste("GS", 1:nC, sep="_")
  
  if (margin=="row") JImat <-1-as.matrix(dist.binary(t(cbind(GS@RowxNumber*1, clustObj@RowxNumber*1)), method=1))^2
  if (margin=="col") JImat <-1-as.matrix(dist.binary((rbind(GS@NumberxCol*1,  clustObj@NumberxCol*1)), method=1))^2 
  if (margin=="both") JImat <-(1-as.matrix(dist.binary((rbind(GS@NumberxCol*1,  clustObj@NumberxCol*1)), method=1))^2)+(1-as.matrix(dist.binary(t(cbind(GS@RowxNumber, clustObj@RowxNumber*1)), method=1))^2)
  JImat <-JImat[-c(1:nC),1:nC]
  colnames(JImat) = GSnames
  #print(JImat)
  
  if (best) {
    maxJIind=vector()
    maxJIind<-apply(JImat, 2, which.max)
    #print(maxJIind)
    JI=vector()
    for (i in 1:nC) JI[i]<- JImat[maxJIind[i],i]  
    #print(JI)
    JI<-round(JI,3)
    res=cbind(n=maxJIind, JI=JI)
    rownames(res) = GSnames
    return(as.data.frame(res))
  }
  
  if(!best) return(as.data.frame(JImat))
})


setGeneric("analyzeClust", function(clustObj, GS,...) standardGeneric("analyzeClust"))
setMethod("analyzeClust", c('iBBiG','iBBiG'),
  function(clustObj,GS, ...) {
  clustObj <-list(clustObj)
  print(class(clustObj))
	analyzeClust(clustObj, GS, ...)
	})

setMethod("analyzeClust", c('Biclust','iBBiG'),
  function(clustObj,GS, ...) {
  clustObj <-list(as(clustObj, "iBBiG"))
  print(class(clustObj))
  analyzeClust(clustObj, GS, ...)
	})
setMethod("analyzeClust", c('list','iBBiG'),
function(clustObj,GS,  margin="col",stats=TRUE, type="bestRun",n=NULL,parameters=NULL) {
  ## This will report the best cluster, JE, size and sensivity, specificity
  ## PPV and NPV of a given cluster (clustObj) and the truth GS
  ## Both clustObj and GS are of class biclust
  ## n is optional and not required, but can be used to give the function res[[20]]
  ## Paramters is any optional text regarding the run, for example test paramters eg parameters=names(res)[n]

  if (!is.null(n))  {
    clustObj<-clustObj[[n]]
    n = rep(n, nC)
  }
  # if (clustObj@Number==0) return(NULL)
  
  
  predStats<-function(rTrue, rPred){
    ## Expect 2 binary vectors with 1=TRUE, 0=FALSE
    ## A true positive (TP)
    ## B false positve (FP)
    ## C false negative (FN)
    ## D true negative (TN)
    ## Precision is the same as PPV
    ## See http://en.wikipedia.org/wiki/Specificity_(tests)
    
    
    rTab<- table(rTrue, rPred)
    #  print(rTab)
    A<-rTab["1","1"]
    B<-rTab["0","1"]
    C<-rTab["1","0"]
    D<-rTab["0","0"]
    #  print(c(A,B,C,D))
    accuracy = sum(A,D)/sum(A,B,C,D)
    precision = A/sum(A,B)    ## Same at the PPV
    sensitivity = A/sum(A,C)
    specificity  = D/sum(B,D)
    PPV = A/sum(A,B)
    NPV =  D/sum(C,D)
    statsRes<-c(accuracy=accuracy,sensitivity=sensitivity,specificity =specificity ,PPV=PPV,NPV=NPV)
    statsRes<-round(statsRes,3)
    
    return(statsRes)
  }
  
  
  clustSize<-function(clustObj, clustInd) {
    nRow<-apply(clustObj@RowxNumber[,clustInd],2,sum)
    nCol<-apply(clustObj@NumberxCol[clustInd,],1, sum)
    return(cbind(nRow=nRow, nCol=nCol))
  }
  
  #print(Sres)
  
  reportPredStats<-function(clustObj, clustInd, GS) {
    nC= GS@Number
    ## ClustInd is the $n, index of cluster from JIdist          
    colStats<-t(sapply(1:nC, function(i) predStats(GS@NumberxCol[i,]*1,(clustObj@NumberxCol[clustInd[i],]>0)*1)))
    colnames(colStats) <- paste("col", colnames(colStats), sep="-")
    
    
    rowStats<-t(sapply(1:nC, function(i) predStats(GS@RowxNumber[,i]*1,clustObj@RowxNumber[,clustInd[i]]*1)))
    colnames(rowStats) <- paste("row", colnames(rowStats), sep="-")
    
    return(cbind(colStats, rowStats))
  }
  
  # Parse List
  if (length(clustObj)==1) unlist(clustObj)
  
  if(is.list(clustObj)) {
    
    resultsJI<-sapply(clustObj, function(x) JIdist(x, GS, margin=margin)$JI)
    resultsJIn<-sapply(clustObj, function(x) JIdist(x, GS, margin=margin)$n)
    
    
    #print(resultsJI)
    
    if(type=="best") {
      ## Merge best clusters from each run
      bC<-apply(resultsJI, 1, which.max)
      n<-bC  
      bCm<- mapply(function(a,b) resultsJIn[a,b], a=1:length(bC), b=bC)
      dupInd<- duplicated(cbind(bC, bCm))
      bC<-bC[!dupInd]
      bCm<-bCm[!dupInd]
      Number = length(bC)
      
      
      RowxNumber<-mapply(function(a,b) clustObj[[a]]@RowxNumber[,b], a=bC, b=bCm, SIMPLIFY=TRUE)
      
      NumberxCol<-t(mapply(function(a,b) clustObj[[a]]@NumberxCol[b,], a=bC, b=bCm, SIMPLIFY=TRUE))
      
      Clusterscores<-mapply(function(a,b) clustObj[[a]]@Clusterscores[b], a=bC, b=bCm, SIMPLIFY=TRUE)
      
      RowScorexNumber<-mapply(function(a,b) clustObj[[a]]@RowScorexNumber[,b], a=bC, b=bCm, SIMPLIFY=TRUE) 
      # print(RowScorexNumber)
      
      
      
      newClustObj<-   res<- new("iBBiG", Seeddata=clustObj[[1]]@Seeddata, Clusterscores=Clusterscores, RowScorexNumber=RowScorexNumber,  NumberxCol = NumberxCol, RowxNumber=RowxNumber, Number = Number, info= list(bC, resultsJI))
      #return(newClustObj)
    }
    
    #if(type=="mean") {
    ## Not implmented would be done, by running lapply analyzeClust on  all object in the list clustObj and then taking the mean of it
    #}
    
    if(type=="bestRun") {
      n<-which.max(apply(resultsJI, 2, sum))
      #  print(n)
      newClustObj<-clustObj[[n]]   
      n<-rep(n, nrow(resultsJI))
      #print(newClustObj)
    }
    
    clustObj<-newClustObj
    
  }
  
  nC= GS@Number
  # print(clustObj)
  
  # Parse individual Biclust obj
  Sres<- JIdist(clustObj, GS, margin=margin, best=TRUE)
  
  clustInd<-Sres$n  
  Sres<-cbind(Sres, clustSize(clustObj, clustInd))
  
  if(stats) Sres<-cbind(Sres, reportPredStats(clustObj, clustInd, GS))
  
  
  if (!is.null(n))  Sres<-cbind(Run=n,Sres)
  if (!is.null(parameters))  Sres<-cbind(Parameters=rep(parameters,nC), Sres)
  #rownames(Sres)= paste("M", 1:nC, sep="")
  return(Sres)
})



