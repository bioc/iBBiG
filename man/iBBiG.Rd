\name{iBBiG}
\alias{iBBiG}

\title{Iterative Binary Bi-Clustering for GeneSets}
\description{
iBBiG is a bi-clustering algorithm which is optimized for clustering binary data resulting from discretized p-values of genomic analyses
}
\usage{
iBBiG(binaryMatrix, nModules, alpha = 0.3, pop_size = 100, mutation = 0.08, stagnation = 50, selection_pressure = 1.2, max_sp = 15, success_ratio = 0.6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{binaryMatrix}{Matrix. A binary or logical matrix.
}
  \item{nModules}{Numeric. The number of expected  modules. As iBBiG is optimized to find a miminal number, nModules can be a larger than expected value
}
  \item{alpha}{Numeric, weighting factor, that will balances the tradeoff between specificity and sensitivity. Default 0.3.  Simulated studies indicate range 0.3-0.5 is appropriate
}
  \item{pop_size}{Numeric. Default 100. Population size establishes the genetic diversity of solutions in Genetic Algorithm. Simulated studies show that it has marginal effect on performance.
}
  \item{mutation}{Numeric. Default 0.08. Mutation rate of GA. Simulated studies show that it has little effect on performance.
}
  \item{stagnation}{Numeric. Default is stop criterion of 50 iterations of stagnation. Simulated studies show that it has little effect on performance.
}
  \item{selection_pressure}{Numeric. Default is 1.2. Selection pressure for parent selection. Simulated studies show that it has little effect on performance}
  \item{max_sp}{ Numeric. Default is 15. Simulated studies show that it has little effect on performance
}


  \item{success_ratio}{
Numeric. Deafult 0.6. Success ratio determines how many children have to outperform at least one of their parents. Simulated studies show that it has little effect on performance
}
}
\details{
iBBiG is a bi-clustering algorithm, optimized for module discovery in sparse noisy binary genomics data.  We designed iBBiG to have high specificity and thereby minimize the false positive rate when discovering new classes; the iterative approach employed in iBBiG is able to discover weak signals, even if they are potentially masked by stronger ones. For a compairions with global clustering approaches (K-means, hierarchical cluster analysis) and bi-clustering approaches (Bimax, FABIA, COALESCE) see our manuscript Gusenleitner et al., 2012.
An advantage of iBBiG relative to other methods is that it does not require a priori knowledge of the true number of clusters. Following the application of iBBiG, the number of true clusters can be estimated from the weighted cluster scores and RowScorexNumber of the extracted modules. In some cases, we observed that a module may represent the residue or remaining signal of a stronger, previously extracted module. This residue remains because iBBiG only removes information from the data matrix that is actually used for the entropy based score in a module. However, we do not consider these residual modules to be a shortcoming of the method as their existence facilitates discovery of the true overlap between modules and, further, these modules can be easily  detected by looking at the overlap of clinical covariates and gene sets.\\
Although iBBiG includes several parameters, we have shown that most impact only computation time, and do not effect cluster discovery. The only parameter that had an impact on cluster discovery was alpha, which is a weighting factor that balances the cost of increasing cluster size (number of rows) against cluster homogeneity. In generating small homogeneous clusters, one might miss information. Conversely, large hetergeneous clusters may contain more false positives. Although alpha does not regulate the number of clusters, decreasing stringency, by increasing alpha values may produce greater numbers of clusters. As a results the alpha parameter is useful in adjusting the sensitivity-specificity ratio. Alpha has a range 0.1-1 where 0.1 will generate fewer, smaller homogeneous clusters whereas 0.9 is less stringent and results in more hetergeneous clusters (with greater potential for false positives). Increasing alpha will generate more clusters of greater size, with potentially greater specificity at the expense of decreased sensitivity. Following tests on simulated data we recommended alpha values between 0.3-0.5 (Gusenleitner et al., 2012).  The default alpha is 0.3
}
\value{
Returns an object with class iBBiG, which extents the class Biclust.

 \item{Seeddata }{Input binaryMatrix}
 \item{RowScorexNumber}{Matrix. Score for each signature (row) in each cluster. Matrix with dimensions, Number of Rows in Seeddata x Number of clusters}
 \item{Clusterscores}{Vector. Score for each cluster. It has length equal to the number of clusters.  } 
 \item{ Parameters}{List of Input Parameters (if provided)}
 \item{RowxNumber}{Binary or Logical Matrix with dimensions, Number of Rows in Seeddata x Number of clusters, where 1 represents cluster membership}
 \item{ NumberxCol}{Binary or Logical Matrix with dimensions, Number of clusters x Number of Columns in Seeddata ,where 1 represents cluster membership}
  \item{Number}{Numeric. Number of modules(clusters)}
  \item{info}{list. which is a general contained for other information. }
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
binMat<-makeArtificial()
plot(binMat)
res<- iBBiG(binMat@Seeddata, nModules=10)
plot(res)
res
analyzeClust(res,binMat)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clustering }
\keyword{ GSEA }
\keyword{ metaanalysis }
\keyword{ biclustering }% __ONLY ONE__ keyword per line
