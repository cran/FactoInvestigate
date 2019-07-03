\name{outliers}
\alias{outliers}

\title{
Outliers detection
}

\description{
Detection of singular individuals that concentrates too much inertia.
}

\usage{
outliers(res, file = "", Vselec = "cos2", Vcoef = 1, nmax = 10, 
         figure.title = "Figure", graph = TRUE, cex = 0.7, options = NULL)
}

\arguments{
  \item{res}{
an object of class PCA or MCA.
}
  \item{file}{
a numerical vector giving the factorial dimensions for with to compute the eigen values calculation.
}
  \item{Vselec}{
the variables to select ; see the details section.
}
  \item{Vcoef}{
a numerical coefficient to adjust the variables selection rule ; see the details section.
}
  \item{nmax}{
an integer giving the maximum number of variables to illustrate each outlier (by default 10).
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
\item{options}{a character string that gives the output options for the figures. 
If NULL, \code{options="r, echo = FALSE, fig.align = 'center', fig.height = 3.5, fig.width = 5.5"} for linux and Mac
and \code{options="r, echo = FALSE, fig.height = 3.5, fig.width = 5.5"} for Windows
}
}

\details{
The algorithm detects an individual as an outlier if its contribution to the plane if higher to 3 standard deviation.
\cr 
\cr The \code{Vselec} argument is used in order to select a part of the variables that are drawn and described. For example, you can use either :
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

\value{
  \item{new.res}{
the res object without the outliers (they are completely eliminated).
}
  \item{res.out}{
the res object with the outliers as supplementary individuals.
}
  \item{memory}{
the original res object.
}
  \item{N}{
the number of outliers.
}
  \item{ID}{
the label of outliers.
}
}

\author{
Simon Thuleau and Francois Husson
}

\examples{
\dontrun{
require(FactoMineR)
data(decathlon)
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
outliers(res.pca, file = "PCA.Rmd")
}
}
\keyword{outlier}
