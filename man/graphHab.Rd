\name{graphHab}
\alias{graphHab}

\title{
Colored factor graph
}
\description{
Realised the graph of individuals colored after a variable categories
}

\usage{
graphHab(res, file = "", dim = 1:2, hab = NULL, ellipse = TRUE, Iselec = "contrib", 
         Rselec = "cos2", Cselec = "contrib", Icoef = 1, Rcoef = 1, Ccoef = 1, 
         figure.title = "Figure", graph = TRUE, cex = 0.7, options = NULL)
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
  \item{hab}{
a variable name or index to use to color the individuals (or rows) among the variable categories.
}
  \item{ellipse}{
a boolean : if \code{TRUE}, ellipses are plotted with the coloration of individuals (or rows).
}
  \item{Iselec}{
the individuals to select ; see the details section.
}
  \item{Rselec}{
the rows to select (for a \code{\link[FactoMineR]{CA} res} object) ; see the details section.
}
  \item{Cselec}{
the columns to select (for a \code{\link[FactoMineR]{CA} res} object) ; see the details section.
}
  \item{Icoef}{
a numerical coefficient to adjust the individuals selection rule ; see the details section.
}
  \item{Rcoef}{
a numerical coefficient to adjust the rows selection rule (for a \code{\link[FactoMineR]{CA} res} object) ; see the details section.
}
  \item{Ccoef}{
a numerical coefficient to adjust the columns selection rule (for a \code{\link[FactoMineR]{CA} res} object) ; see the details section.
}
  \item{figure.title}{
the text label to add before graph title.
}
  \item{graph}{
a boolean : if \code{TRUE}, graphs are ploted.
}
  \item{cex}{
an optional argument for the generic plot functions, used to adjust the size of the elements plotted.
}
\item{options}{a character string that gives the output options fir the figures. 
If NULL, \code{options="r, echo = FALSE, fig.align = 'center', fig.height = 3.5, fig.width = 5.5"} for linux and Mac
and \code{options="r, echo = FALSE, fig.height = 3.5, fig.width = 5.5"} for Windows
}
}

\details{
The \code{Iselec} argument (respectively \code{Rselec} or \code{Cselec}) is used in order to select a part of the elements that are drawn and described. For example, you can use either :
\cr - \code{Iselec = 1:5} then the individuals (respectively the rows or the columns) numbered 1 to 5 are drawn.
\cr - \code{Iselec = c("name1","name5")} then the individuals (respectively the rows or the columns) named \code{name1} and \code{name5} are drawn.
\cr - \code{Iselec = "contrib 10"} then the 10 active or illustrative individuals (respectively the rows or the columns) that have the highest contribution on the 2 dimensions of the plane are drawn.
\cr - \code{Iselec = "contrib"} then the optimal number of active or illustrative individuals (respectively the rows or the columns) that have the highest contribution on the 2 dimensions of the plane are drawn.
\cr - \code{Iselec = "cos2 5"} then the 5 active or illustrative individuals (respectively the rows or the columns) that have the highest cos2 on the 2 dimensions of the plane are drawn.
\cr - \code{Iselec = "cos2 0.8"} then the active or illustrative individuals (respectively the rows or the columns) that have a \code{cos2} higher to \code{0.8} on the plane are drawn.
\cr - \code{Iselec = "cos2"} then the optimal number of active or illustrative individuals (respectively the rows or the columns) that have the highest cos2 on the 2 dimensions of the plane are drawn.
\cr 
\cr The \code{Icoef} argument (respectively \code{Rcoef} or \code{Ccoef}) is used in order to adjust the selection of the elements when based on \code{Iselec = "contrib"} or \code{Iselec = "cos2"}. For example :
\cr - if \code{Icoef = 2}, the threshold is 2 times higher, and thus 2 times more restrictive.
\cr - if \code{Icoef = 0.5}, the threshold is 2 times lower, and thus 2 times less restrictive.
}


\author{
Simon Thuleau and Francois Husson
}

\seealso{
  \code{\link{factoGraph}}, \code{\link{graphInd}}, \code{\link{graphCA}}, \code{\link{graphVar}}, \code{\link{graphSup}}
}

\examples{
\dontrun{
require(FactoMineR)
data(decathlon)
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
graphHab(res.pca)
}
}
\keyword{graphs}
