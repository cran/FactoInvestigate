\name{selection}
\alias{selection}

\title{
Graphical elements selection
}

\description{
Select the best elements to plot in a graph
}

\usage{
selection(res, dim = 1:2, margin = 1, selec = "cos2", coef = 1)
}

\arguments{
  \item{res}{
an object of class PCA, CA or MCA.
}
  \item{dim}{
a 2 dimensional numerical vector giving the factorial dimensions to use for the representation (by default the first plane).
}
  \item{margin}{
an integer (by default 1). If equals 1, the function computes on the individuals (or rows). If equals 2, the function computes on the active variables (or columns). If equals 3, the function computes on the supplementary variables.
}
  \item{selec}{
the elements to select ; see the details section.
}
  \item{coef}{
a numerical coefficient to adjust the elements selection rule ; see the details section.
}
}

\details{
The \code{selec} argument is used in order to select a part of the elements that are drawn and described. For example, you can use either :
\cr - \code{selec = 1:5} then the elements numbered 1 to 5 are drawn.
\cr - \code{selec = c("name1","name5")} then the elements named \code{name1} and \code{name5} are drawn.
\cr - \code{selec = "contrib 10"} then the 10 active or illustrative elements that have the highest contribution on the 2 dimensions of the plane are drawn.
\cr - \code{selec = "contrib"} then the optimal number of active or illustrative elements that have the highest contribution on the 2 dimensions of the plane are drawn.
\cr - \code{selec = "cos2 5"} then the 5 active or illustrative elements that have the highest cos2 on the 2 dimensions of the plane are drawn.
\cr - \code{selec = "cos2 0.8"} then the active or illustrative elements that have a \code{cos2} higher to \code{0.8} on the plane are drawn.
\cr - \code{selec = "cos2"} then the optimal number of active or illustrative elements that have the highest cos2 on the 2 dimensions of the plane are drawn.
\cr 
\cr The \code{coef} argument is used in order to adjust the selection of the elements when based on \code{selec = "contrib"} or \code{selec = "cos2"}. For example :
\cr - if \code{coef = 2}, the threshold is 2 times higher, and thus 2 times more restrictive.
\cr - if \code{coef = 0.5}, the threshold is 2 times lower, and thus 2 times less restrictive.
}

\value{
  \item{drawn}{the elements selected.}
  \item{what.drawn}{the criterion of selection (as a sentence).}
}

\author{
Simon Thuleau and Francois Husson
}

\seealso{
  \code{\link{description}}
}

\examples{
\dontrun{
require(FactoMineR)
data(decathlon)
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
selection(res.pca, margin = 1, selec = "contrib 10")
}
}
\keyword{selection}
