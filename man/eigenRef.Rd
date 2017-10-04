\name{eigenRef}
\alias{eigenRef}

\title{
Reference eigen value
}

\description{
Compute the eigen values of random datasets, with the hypothesis of independence.
}

\usage{
eigenRef(res, dim = NULL, q = 0.95, time = "10000L", parallel = TRUE)
}

\arguments{
  \item{res}{
an object of class PCA, CA or MCA
}
  \item{dim}{
a numerical vector giving the factorial dimensions for with to compute the eigenvalues calculation.
}
  \item{q}{
the quantile of computed values to use as reference value (ie. the confidence about the signification of dimensions)
}
  \item{time}{
a character indicating the loop condition. This string is made of a number and a letter coupled. The number X with letter \code{L} means to compute X datasets exactly. The number X with letter \code{s} means to compute as many datasets as possible during approximativley X seconds.
}
  \item{parallel}{
a boolean : if \code{TRUE}, the computation uses map reduce on the processor cores to increase the performance. Useful for huge datasets
}
}

\value{
  \item{datasets}{the number of random datasets simulated.}
  \item{quantile}{the quantile used for the reference definition.}
  \item{inertia}{the reference inertia for the dimensions declared.}
}

\author{
Simon Thuleau and Francois Husson
}

\seealso{
 \code{\link{dimRestrict}}, \code{\link{inertiaDistrib}}
}

\examples{
\dontrun{
data(decathlon)
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
eigenRef(res.pca, q = 0.95, time = "10s")

data(children)
res.ca = CA(children, row.sup = 15:18, col.sup = 6:8, graph = FALSE)
eigenRef(res.ca, q = 0.99, time = "10000L")

data(tea)
res.mca = MCA(tea, quanti.sup = 19,quali.sup = 20:36, graph = FALSE)
eigenRef(res.mca, dim = 1:8, q = 0.90, time = "10s")
}
}

\keyword{inertia}
