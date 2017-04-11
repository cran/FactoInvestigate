\name{graphSup}
\alias{graphSup}

\title{
Supplementary variables factor map
}
\description{
Realise the optimised graph of supplementary variables
}

\usage{
graphSup(res, file = "", dim = 1:2, Mselec = "cos2", Mcoef = 1, 
         figure.title = "Figure", graph = TRUE, cex = 0.7)
}

\arguments{
  \item{res}{
an object of class PCA, CA or MCA.
}
  \item{file}{
the file path where to write the description in Rmarkdown langage. If not specified, the description is written in the console.
}
  \item{dim}{
a 2 dimensional numerical vector giving the factorial dimensions to use for the representation (by default the first plane).
}
  \item{Mselec}{
the supplementary variables to select ; see the details section.
}
  \item{Mcoef}{
a numerical coefficient to adjust the supplementary variables selection rule ; see the details section.
}
  \item{figure.title}{
the text label to add before graph title.
}
  \item{graph}{
a boolean : if \code{TRUE}, graphs are ploted.
}
  \item{cex}{
an optional argument for the generic plot functions, used to adjust the size of the elements ploted.
}
}

\details{
The \code{Mselec} argument is used in order to select a part of the illustrative variables that are drawn and described. For example, you can use either :
\cr - \code{Mselec = 1:5} then the illustrative variables numbered 1 to 5 are drawn.
\cr - \code{Mselec = c("name1","name5")} then the illustrative variables named \code{name1} and \code{name5} are drawn.
\cr - \code{Mselec = "cos2 5"} then the 5 illustrative variables that have the highest cos2 on the 2 dimensions of the plane are drawn.
\cr - \code{Mselec = "cos2 0.8"} then the illustrative variables that have a \code{cos2} higher to \code{0.8} on the plane are drawn.
\cr - \code{Mselec = "cos2"} then the optimal number of illustrative variables that have the highest cos2 on the 2 dimensions of the plane are drawn.
\cr 
\cr The \code{Mcoef} argument is used in order to adjust the selection of the illustrative variables when based on \code{Mselec = "cos2"}. For example :
\cr - if \code{Mcoef = 2}, the threshold is 2 times higher, and thus 2 times more restrictive.
\cr - if \code{Mcoef = 0.5}, the threshold is 2 times lower, and thus 2 times less restrictive.
}


\author{
Simon THULEAU, with the participation of François HUSSON
}

\seealso{
  \code{\link{factoGraph}}, \code{\link{graphInd}}, \code{\link{graphHab}}, \code{\link{graphCA}}, \code{\link{graphVar}}
}

\examples{
\dontrun{
data(decathlon)
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
graphSup(res.pca)

data(children)
res.ca = CA(children, row.sup = 15:18, col.sup = 6:8, graph = FALSE)
graphSup(res.ca)

data(tea)
res.mca = MCA(tea, quanti.sup = 19,quali.sup = 20:36, graph = FALSE)
graphSup(res.mca)
}
}

\keyword{graphs}
