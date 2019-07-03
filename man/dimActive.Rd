\name{dimActive}
\alias{dimActive}

\title{
Number of active dimensions
}

\description{
Give the number of active elements used to build the factorial analysis : individuals (or rows) and variables (or columns)
}

\usage{
dimActive(res)
}

\arguments{
  \item{res}{
an object of class PCA, CA or MCA.
}
}

\author{
Simon Thuleau and Francois Husson
}

\examples{
\dontrun{
require(FactoMineR)
data(decathlon)
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
dimActive(res.pca)
}
}
\keyword{dimension}
