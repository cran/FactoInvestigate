\name{getParam}
\alias{getParam}

\title{
Factorial parameters
}
\description{
Get all the factorial object parameters
}

\usage{
getParam(res)
}


\arguments{
  \item{res}{
an object of class PCA, CA or MCA.
}
}


\value{
  \item{data}{
the dataset.
}
  \item{ind}{
the number of individuals.
}
  \item{var}{
the number of variables.
}
  \item{row}{
the number of rows (CA).
}
  \item{col}{
the number of columns (CA).
}
  \item{ind.sup}{
the number of supplementary individuals.
}
  \item{quanti.sup}{
the number of quantitative supplementary variables.
}
  \item{quali.sup}{
the number of qualitative supplementary variables.
}
  \item{row.sup}{
the number of supplementary rows (CA).
}
  \item{col.sup}{
the number of supplementary columns (CA).
}
  \item{row.w}{
the weights of each row.
}
  \item{col.w}{
the weights of each columns.
}
  \item{scale}{
a boolean indicatting if the variables are scaled or not.
}
  \item{ncp.mod}{
the number of component kept in the analysis object.
}
  \item{modalites}{
the list of factors for each qualitative variables.
}
}


\author{
Simon THULEAU, with the participation of François HUSSON
}

\seealso{
  \code{\link{whichFacto}}
}

\examples{
\dontrun{
data(decathlon)
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
getParam(res.pca)

data(children)
res.ca = CA(children, row.sup = 15:18, col.sup = 6:8, graph = FALSE)
getParam(res.ca)

data(tea)
res.mca = MCA(tea, quanti.sup = 19,quali.sup = 20:36, graph = FALSE)
getParam(res.mca)
}
}

\keyword{parameter}
