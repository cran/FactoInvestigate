\name{classif}
\alias{classif}

\title{
Classification description
}

\description{
Realise the hierachical ascending classification (\code{\link{HCPC}} function) of the individuals (or rows) and describe the specifications of each cluster.
}

\usage{
classif(res, file = "", dim = 1:2, nclust = -1, selec = "contrib", coef = 1, 
        mmax = 10, nmax = 10, figure.title = "Figure", graph = TRUE, options = NULL)
}

\arguments{
  \item{res}{
an object of class PCA, CA or MCA.
}
  \item{file}{
the file path where to write the description in Rmarkdown langage. If not specified, the description is written in the console.
}
  \item{dim}{
a 2 dimensional numerical vector giving the factorial dimensions to describe (by default the first plane).
}
  \item{nclust}{
an integer forcing the number of clusters desired. The value -1 return the clustering evaluated as the most appropriate (default).
}
  \item{selec}{
the selection criterion of individuals to plot on the graph.
}
  \item{coef}{
a numerical coefficient to adjust the selection rule (exemple : if equals 2, the threshold is 2 times higher, and thus more restrictive)
}
  \item{mmax}{
an integer giving the maximum number of individuals to illustrate each cluster (by defaut 10).
}
  \item{nmax}{
an integer giving the maximum number of variables to illustrate each cluster (by defaut 10).
}
  \item{figure.title}{
the text label to add before graph title.
}
  \item{graph}{
a boolean : if \code{TRUE}, the graph is ploted into the console.
}
\item{options}{a character string that gives the output options for the figures. 
If NULL, \code{options="r, echo = FALSE, fig.align = 'center', fig.height = 3.5, fig.width = 5.5"} for linuw and Mac
and \code{options="r, echo = FALSE, fig.height = 3.5, fig.width = 5.5"} for Windows
}
}

\details{
The \code{selec} argument is used in order to select a part of the individuals that are drawn and described. For example, you can use either :
\cr - \code{selec = 1:5} then the individuals numbered 1 to 5 are drawn.
\cr - \code{selec = c("name1","name5")} then the individuals named \code{name1} and \code{name5} are drawn.
\cr - \code{selec = "contrib 10"} then the 10 active or illustrative individuals that have the highest contribution on the 2 dimensions of the plane are drawn.
\cr - \code{selec = "contrib"} then the optimal number of active or illustrative individuals that have the highest contribution on the 2 dimensions of the plane are drawn.
\cr - \code{selec = "cos2 5"} then the 5 active or illustrative individuals that have the highest cos2 on the 2 dimensions of the plane are drawn.
\cr - \code{selec = "cos2 0.8"} then the active or illustrative individuals that have a \code{cos2} higher to \code{0.8} on the plane are drawn.
\cr - \code{selec = "cos2"} then the optimal number of active or illustrative individuals that have the highest cos2 on the 2 dimensions of the plane are drawn.
\cr 
\cr The \code{coef} argument is used in order to adjust the selection of the individuals when based on \code{selec = "contrib"} or \code{selec = "cos2"}. For example :
\cr - if \code{coef = 2}, the threshold is 2 times higher, and thus 2 times more restrictive.
\cr - if \code{coef = 0.5}, the threshold is 2 times lower, and thus 2 times less restrictive.
}

\value{
\item{res.hcpc}{the result of the \code{\link{HCPC}} function on the dimensions specified.}
}

\author{
Simon Thuleau and Francois Husson
}

\seealso{
  \code{\link{description}}
}

\examples{
\dontrun{
data(decathlon)
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
classif(res.pca, file = "PCA.Rmd")

data(children)
res.ca = CA(children, row.sup = 15:18, col.sup = 6:8, graph = FALSE)
classif(res.ca, file = "CA.Rmd")

data(tea)
res.mca = MCA(tea, quanti.sup = 19,quali.sup = 20:36, graph = FALSE)
classif(res.mca, file = "MCA.Rmd")
}
}

\keyword{classification}