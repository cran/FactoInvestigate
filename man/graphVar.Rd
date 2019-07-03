\name{graphVar}
\alias{graphVar}

\title{
Variables factor map
}
\description{
Realise the optimised variables graph
}

\usage{
graphVar(res, file = "", dim = 1:2, Vselec = "cos2", Vcoef = 1, 
         figure.title = "Figure", graph = TRUE, cex = 0.7, options=NULL)
}



\arguments{
  \item{res}{
an object of class PCA or MCA.
}
  \item{file}{
the file path where to write the description in Rmarkdown language. If not specified, the description is written in the console.
}
  \item{dim}{
a 2 dimensional numerical vector giving the factorial dimensions to use for the representation (by default the first plane).
}
  \item{Vselec}{
the variables to select ; see the details section.
}
  \item{Vcoef}{
a numerical coefficient to adjust the variables selection rule ; see the details section.
}
  \item{figure.title}{
the text label to add before graph title.
}
  \item{graph}{
a boolean : if \code{TRUE}, graphs are plotted.
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
The \code{Vselec} argument is used in order to select a part of the variables that are drawn and described. For example, you can use either :
\cr - \code{Vselec = 1:5} then the variables numbered 1 to 5 are drawn.
\cr - \code{Vselec = c("name1","name5")} then the variables named \code{name1} and \code{name5} are drawn.
\cr - \code{Vselec = "contrib 10"} then the 10 active or illustrative variables that have the highest contribution on the 2 dimensions of the plane are drawn.
\cr - \code{Vselec = "contrib"} then the optimal number of active or illustrative variables that have the highest contribution on the 2 dimensions of the plane are drawn.
\cr - \code{Vselec = "cos2 5"} then the 5 active or illustrative variables that have the highest cos2 on the 2 dimensions of the plane are drawn.
\cr - \code{Vselec = "cos2 0.8"} then the active or illustrative variables that have a \code{cos2} higher to \code{0.8} on the plane are drawn.
\cr - \code{Vselec = "cos2"} then the optimal number of active or illustrative variables that have the highest cos2 on the 2 dimensions of the plane are drawn.
\cr 
\cr The \code{Vcoef} argument is used in order to adjust the selection of the variables when based on \code{Vselec = "contrib"} or \code{Vselec = "cos2"}. For example :
\cr - if \code{Vcoef = 2}, the threshold is 2 times higher, and thus 2 times more restrictive.
\cr - if \code{Vcoef = 0.5}, the threshold is 2 times lower, and thus 2 times less restrictive.
}


\author{
Simon Thuleau and Francois Husson
}

\seealso{
  \code{\link{factoGraph}}, \code{\link{graphInd}}, \code{\link{graphHab}}, \code{\link{graphCA}}, \code{\link{graphSup}}
}

\examples{
\dontrun{
require(FactoMineR)
data(decathlon)
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
graphVar(res.pca)
}
}
\keyword{graphs}
