\name{whichFacto}
\alias{whichFacto}

\title{
Analysis class
}
\description{
Return the class of the factorial object (ie. the kind of analysis performed)
}

\usage{
whichFacto(res)
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
require(FactoMineR)
data(decathlon)
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
whichFacto(res.pca)
}

\keyword{class}
