\name{makeArtificial}
\alias{makeSimDesignMat}
\alias{addSignal}
\alias{makeArtificial}


\title{Create a 400x400 simulated binary matrix for testing iBBiG and other binary biclustering methods}
\description{
Create a binary matrix of 400 rows x 400 columns, where 1 is a positive association. This matrix is seeded with 7 modules of various size and with various levels of noise as described by Gusenleitner et al., 

}
\usage{
makeArtificial(nRow = 400, nCol = 400, noise = 0.1, verbose = TRUE, dM = makeSimDesignMat(verbose = verbose))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nRow}{Numeric \code{nRow} number of rows
}
  \item{nCol}{
Numeric \code{nRow} number of columns
}
  \item{noise}{
Numeric. Value between 0-1. Default is 10 percent random noise (1) introduced into the spare binary matrix
}
  \item{verbose}{
Verbose output. Default is TRUE
}
  \item{dM}{
A design matrix specifying where the columns are. The function makeSimDesignMat create the matrix which specifies the design matrix
}
}
\details{
See Guesnleitner et al, for more information
}
\value{
Output is a class of Biclust.
}
\references{
Daniel Gusenleitner, Eleanor A Howe, Stefan Bentink, John Quackenbush and Aedin C Culhane
iBBiG: Iterative Binary Bi-clustering of Gene Sets
Bioinformatics. In review.
}
\author{
Aedin Culhane, Daniel Gusenleitner
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Further functions for viewing and clustering binaray data are available in the package biclust.  We have written iBBiG and its classes so that it is compatible with biclust, and the class iBBiG inherits Biclust-class.
}


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
arti<-makeArtificial()
plot(arti)
 }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{clustering}
\keyword{GSEA}
\keyword{metaanalysis }
\keyword{biclustering}