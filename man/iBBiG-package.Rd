\name{iBBiG-package}
\alias{iBBiG-package}
\docType{package}
\title{
iBBiG performs bi-clustering of binary matrices
}
\description{
iBBiG is a bi-clustering algorithm, optimized for module discovery in sparse noisy binary genomics data.  
We designed iBBiG to have high specificity and thereby minimize the false positive rate when discovering new classes; 
the iterative approach employed in iBBiG is able to discover weak signals, even if they are potentially masked by stronger ones. 

}
\details{
\tabular{ll}{
Package: \tab iBBiG\cr
Type: \tab Package\cr
Version: \tab 0.99.1\cr
Date: \tab 2012-03-15\cr
License: \tab Free Artistic\cr
LazyLoad: \tab yes\cr
Depends: \tab methods\cr
}
The main functions is iBBiG.  This is the biclustering algorithm.  
}
\author{
Aedin Culhane, Daniel Gusenleitner

Maintainer: Aedin <aedin@jimmy.harvard.edu>

}
\references{
Daniel Gusenleitner, Eleanor A Howe, Stefan Bentink, John Quackenbush and Aedin C Culhane
iBBiG: Iterative Binary Bi-clustering of Gene Sets
Bioinformatics. In review.
}

\keyword{ package }
\keyword{ clustering }
\keyword{ GSEA }
\keyword{ metaanalysis }
\keyword{ biclustering }% __ONLY ONE__ keyword per line

\seealso{
Also see \code{\link[biclust]{biclust}} ~~
}
\examples{


#create simulated datasets
binMat<-makeArtificial()
binMat
plot(binMat)
res<- try(iBBiG(binMat@Seeddata, nModules=10))
plot(res)
res

## Subset a cluster

res[4]
res[1:2]


## As iBBiG extends the class Biclust can use Biclust functions on it
## View the rows and columns of an iBBiG object

## Create a list of matrices, one for each cluster
Modules<-bicluster(res@Seeddata, res)
length(Modules)
lapply(Modules, dim)


# Or extract a list of a specific cluster
M1<-bicluster(res@Seeddata, res, 1)
dim(M1[[1]])
str(M1)
M1[[1]][1:5,1:3]

}
