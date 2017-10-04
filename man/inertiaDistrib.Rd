\name{inertiaDistrib}
\alias{inertiaDistrib}

\title{
Inertia distribution analysis
}

\description{
Analysis of the inertia distribution among each axis, the amount and the significativity
}

\usage{
inertiaDistrib(res, file = "", ncp = NULL, q = 0.95, time = "10000L", 
        parallel = TRUE, figure.title = "Figure", graph = TRUE, options = NULL)
}

\arguments{
  \item{res}{
an object of class PCA, CA or MCA.
}
  \item{file}{
the file path where to write the description in Rmarkdown language. If not specified, the description is written in the console.
}
  \item{ncp}{
an integer to force the number of dimension to analyse.
}
  \item{q}{
the quantile of computed values to use as reference value (ie. the confidence about the signification of dimensions).
}
  \item{time}{
a character indicating the loop condition. This string is made of a number and a letter coupled. The number X with letter \code{L} means to compute X datasets exactly. The number X with letter \code{s} means to compute as many datasets as possible during approximativley X seconds.
}
  \item{parallel}{
a boolean : if \code{TRUE}, the computation uses map reduce on the processor cores to increase the performance. Useful for huge datasets.
}
  \item{figure.title}{
the text label to add before graph title.
}
  \item{graph}{
a boolean : if \code{TRUE}, graphs are plotted.
}
\item{options}{a character string that gives the output options fir the figures. 
If NULL, \code{options="r, echo = FALSE, fig.align = 'center', fig.height = 3.5, fig.width = 5.5"} for linux and Mac
and \code{options="r, echo = FALSE, fig.height = 3.5, fig.width = 5.5"} for Windows
}
}


\value{
  \item{ncp}{the number of significant dimensions (or the dimensions kept).}
}

\author{
Simon Thuleau and Francois Husson
}

\seealso{
 \code{\link{dimRestrict}}, \code{\link{eigenRef}}
}

\examples{
\dontrun{
data(decathlon)
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
inertiaDistrib(res.pca, q = 0.95, time = "10s")

data(children)
res.ca = CA(children, row.sup = 15:18, col.sup = 6:8, graph = FALSE)
inertiaDistrib(res.ca, q = 0.99, time = "10000L")

data(tea)
res.mca = MCA(tea, quanti.sup = 19,quali.sup = 20:36, graph = FALSE)
inertiaDistrib(res.mca, dim = 1:8, q = 0.90, time = "10s")
}
}

\keyword{inertia}