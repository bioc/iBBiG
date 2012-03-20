\name{iBBiG-class}
\Rdversion{1.1}
\docType{class}
\alias{iBBiG-class}
\alias{RowScorexNumber}
\alias{RowScorexNumber<-}
\alias{Clusterscores}
\alias{Clusterscores<-}
\alias{Seeddata}
\alias{Seeddata<-}
\alias{Parameters}
\alias{Parameters<-}
\alias{RowxNumber}
\alias{RowxNumber<-}
\alias{NumberxCol}
\alias{NumberxCol<-}
\alias{Number}
\alias{Number<-}
\alias{info}
\alias{info<-}
\alias{analyzeClust}
\alias{JIdist}
\alias{RowScorexNumber,iBBiG-method}
\alias{Clusterscores,iBBiG-method}
\alias{Seeddata,iBBiG-method}
\alias{Parameters,iBBiG-method}
\alias{RowxNumber,iBBiG-method}
\alias{NumberxCol,iBBiG-method}
\alias{Number,iBBiG-method}
\alias{info,iBBiG-method}
\alias{plot,iBBiG,ANY-method}
\alias{show,iBBiG-method}
\alias{summary,iBBiG-method}
\alias{[,iBBiG-method}
\alias{analyzeClust,iBBiG,iBBiG-method}
\alias{JIdist,iBBiG,iBBiG-method}
  
\title{Class \code{"iBBiG"}}
\description{
Class to contain and describe result of iBBiG Anlaysis
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("iBBiG", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Seeddata}:}{Input binaryMatrix}
    \item{\code{RowScorexNumber}:}{Matrix. Score for each signature (row) in each cluster. Matrix with dimensions, Number of Rows in Seeddata x Number of clusters}
    \item{\code{Clusterscores}:}{Vector. Score for each cluster. It has length equal to the number of clusters} 
    \item{\code{Parameters}:}{List of Input Parameters (if provided)}
    \item{\code{RowxNumber}:}{Binary or Logical Matrix with dimensions, Number of Rows in Seeddata x Number of clusters, where 1 represents cluster membership}
    \item{\code{NumberxCol}:}{Binary or Logical Matrix with dimensions, Number of clusters x Number of Columns in Seeddata ,where 1 represents cluster membership}
    \item{\code{Number}:}{Numeric. Number of modules(clusters)}
    \item{\code{info}:}{list. which is a general contained for other information. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{Biclust}"}, directly.
}
\section{Methods}{
  \describe{
    \item{RowScorexNumber}{\code{signature(x = "iBBiG")}: Returns the row scores fore each cluster. }
    \item{Clusterscores}{\code{signature(x = "iBBiG")}: Returns the overall score for each cluster. }
    \item{Seeddata}{\code{signature(x = "iBBiG")}: Returns the original binary matrix, the clustering is based on. }
    \item{Parameters}{\code{signature(x = "iBBiG")}: Returns parameter sets, inhereted from biclust. }
    \item{RowxNumber}{\code{signature(x = "iBBiG")}: Returns a logical matrix indicating, which rows are included in each bicluster. }
    \item{NumberxCol}{\code{signature(x = "iBBiG")}: Returns a logical matrix indicating, which columns are included in each bicluster. }
    \item{Number}{\code{signature(x = "iBBiG")}: Returns the number of biclusters contained in the iBBiG object.}
    \item{info}{\code{signature(x = "iBBiG")}: Returns additional information on the particular iBBiG object, inhereted from biclust.}
    \item{plot}{\code{signature(x = "iBBiG")}: Plot the iBBiG clustering.}
    \item{show}{\code{signature(object = "iBBiG")}: Shows the Biclusters. }
    \item{summary}{\code{signature(object = "iBBiG")}: Summary of found bi-clusters. }
    \item{[}{\code{signature(object = "iBBiG")}: ... }
    \item{JIdist}{\code{signature(object = "iBBiG")}: ... }
    \item{analyzeClust}{\code{signature(object = "iBBiG")}: ... }
    }

}

\references{
Daniel Gusenleitner, Eleanor A Howe, Stefan Bentink, John Quackenbush and Aedin C Culhane
iBBiG: Iterative Binary Bi-clustering of Gene Sets
Bioinformatics. In review.
}
\author{
Aedin Culhane, Daniel Gusenleitner
}

\seealso{
Further functions for viewing and clustering binary data are available in the package biclust.  We have written iBBiG and its classes so that it is compatible with biclust, and the class iBBiG inherits Biclust-class.
}

\examples{
showClass("iBBiG")
}
\keyword{clustering}
\keyword{GSEA}
\keyword{metaanalysis}
\keyword{biclustering}