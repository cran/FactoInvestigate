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
Simon THULEAU, with the participation of François HUSSON
}

\examples{
\dontrun{
data(decathlon)
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
whichFacto(res.pca)

data(children)
res.ca = CA(children, row.sup = 15:18, col.sup = 6:8, graph = FALSE)
whichFacto(res.ca)

data(tea)
res.mca = MCA(tea, quanti.sup = 19,quali.sup = 20:36, graph = FALSE)
whichFacto(res.mca)
}
}

\keyword{class}
