\name{dimRestrict}
\alias{dimRestrict}

\title{
Significant dimensions identification
}

\description{
Evaluate the number of significant dimensions in the data.
}

\usage{
dimRestrict(res, file = "", rand = NULL)
}


\arguments{
  \item{res}{
an object of class PCA, CA or MCA.
}
  \item{file}{
the file path where to write the function execution in Rmarkdown language. If not specified, the description is written in the console.
}
  \item{rand}{
an optional vector of eigenvalues to compare the observation with. If NULL, use the result of the \code{\link{eigenRef}} function for comparison.
}
}

\value{
  \item{ncp}{the number of significant dimensions.}
}

\author{
Simon Thuleau and Francois Husson
}

\seealso{
  \code{\link{eigenRef}}, \code{\link{inertiaDistrib}}
}

\examples{
\dontrun{
require(FactoMineR)
data(decathlon)
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
dimRestrict(res.pca, file = "PCA.Rmd")
}
}
\keyword{inertia}
