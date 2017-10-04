\name{scriptRmd}
\alias{scriptRmd}

\title{
File script
}

\description{
Read the script of a file and return each line as a character chain
}

\usage{
scriptRmd(file, output = "code.R")
}

\arguments{
  \item{file}{
the file path to read.
}
  \item{output}{
the file path to write the R code.
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
create.rmd(res.pca, file = "PCA.Rmd", document = "pdf_document")
scriptRmd(file = "PCA.Rmd")

data(children)
res.ca = CA(children, row.sup = 15:18, col.sup = 6:8, graph = FALSE)
create.rmd(res.ca, file = "CA.Rmd", document = "html_document")
scriptRmd(file = "CA.Rmd")

data(tea)
res.mca = MCA(tea, quanti.sup = 19,quali.sup = 20:36, graph = FALSE)
create.rmd(res.mca, file = "MCA.Rmd", document = c("word_document", "pdf_document"))
scriptRmd(file = "MCA.Rmd")
}
}

\keyword{script}
