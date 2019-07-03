\name{description}
\alias{description}

\title{
Factorial dimensions description
}

\description{
Describe a couple of dimensions, giving the individuals specific to each dimension, and the variables that characterize each group of individuals.
}

\usage{
description(res, file = "", dim = 1:2, desc = dim, Iselec = "contrib", 
            Vselec = "cos2", Rselec = "cos2", Cselec = "cos2", Icoef = 1, 
            Vcoef = 1, Rcoef = 1, Ccoef = 1, mmax = 10, nmax = 10)
}

\arguments{
  \item{res}{
an object of class PCA, CA or MCA.
}
  \item{file}{
the file path where to write the description in Rmarkdown language. If not specified, the description is written in the console.
}
  \item{dim}{
a 2 dimensional numerical vector giving the factorial dimensions to use for the representation (by default the first plane).
}
  \item{desc}{
a 2 dimensional numerical vector giving the factorial dimensions to describe (by default the \code{dim} value).
}
  \item{Iselec}{
the individuals to select; see the details section.
}
  \item{Vselec}{
the variables to select; see the details section.
}
  \item{Rselec}{
the rows to select (for a \code{\link{CA} res} object); see the details section.
}
  \item{Cselec}{
the columns to select (for a \code{\link{CA} res} object); see the details section.
}
  \item{Icoef}{
a numerical coefficient to adjust the individuals selection rule; see the details section.
}
  \item{Vcoef}{
a numerical coefficient to adjust the variables selection rule; see the details section.
}
  \item{Rcoef}{
a numerical coefficient to adjust the rows selection rule (for a \code{\link{CA} res} object); see the details section.
}
  \item{Ccoef}{
a numerical coefficient to adjust the columns selection rule (for a \code{\link{CA} res} object); see the details section.
}
  \item{mmax}{
an integer giving the maximum number of individuals to illustrate each cluster (by default 10).
}
  \item{nmax}{
an integer giving the maximum number of variables to illustrate each cluster (by default 10).
}
}

\details{
The \code{Iselec} argument (respectively \code{Vselec}, \code{Rselec} or \code{Cselec}) is used in order to select a part of the elements that are drawn and described. For example, you can use either :
\cr - \code{Iselec = 1:5} then the individuals (respectively the variables, the rows or the columns) numbered 1 to 5 are drawn.
\cr - \code{Iselec = c("name1","name5")} then the individuals (respectively the variables, the rows or the columns) named \code{name1} and \code{name5} are drawn.
\cr - \code{Iselec = "contrib 10"} then the 10 active or illustrative individuals (respectively the variables, the rows or the columns) that have the highest contribution on the 2 dimensions of the plane are drawn.
\cr - \code{Iselec = "contrib"} then the optimal number of active or illustrative individuals (respectively the variables, the rows or the columns) that have the highest contribution on the 2 dimensions of the plane are drawn.
\cr - \code{Iselec = "cos2 5"} then the 5 active or illustrative individuals (respectively the variables, the rows or the columns) that have the highest cos2 on the 2 dimensions of the plane are drawn.
\cr - \code{Iselec = "cos2 0.8"} then the active or illustrative individuals (respectively the variables, the rows or the columns) that have a \code{cos2} higher to \code{0.8} on the plane are drawn.
\cr - \code{Iselec = "cos2"} then the optimal number of active or illustrative individuals (respectively the variables, the rows or the columns) that have the highest cos2 on the 2 dimensions of the plane are drawn.
\cr 
\cr The \code{Icoef} argument (respectively \code{Vcoef}, \code{Rcoef} or \code{Ccoef}) is used in order to adjust the selection of the elements when based on \code{Iselec = "contrib"} or \code{Iselec = "cos2"}. For example :
\cr - if \code{Icoef = 2}, the threshold is 2 times higher, and thus 2 times more restrictive.
\cr - if \code{Icoef = 0.5}, the threshold is 2 times lower, and thus 2 times less restrictive.
}

\author{
Simon Thuleau and Francois Husson
}

\seealso{
  \code{\link{classif}}
}

\examples{
\dontrun{
require(FactoMineR)
data(decathlon)
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
description(res.pca, file = "PCA.Rmd", dim = 1:2)
}
}
\keyword{description}
