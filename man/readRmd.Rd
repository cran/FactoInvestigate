\name{readRmd}
\alias{readRmd}

\title{
Read Rmarkdown file
}

\description{
Compile and open a Rmarkdown file.
}

\usage{
readRmd(file, document = "html_document")
}

\arguments{
  \item{file}{
the file path where to write the description in Rmarkdown langage. If not specified, the description is written in the console.
}
  \item{document}{
a character vector giving the document format desired between "word_document", "pdf_document" and "html_document". This have to be any of those indicated in the file config (by \code{\link{createRmd}}).
}
}

\author{
Simon Thuleau and Francois Husson
}

\seealso{
  \code{\link{createRmd}}, \code{\link{writeRmd}}
}
\examples{
\dontrun{
require(FactoMineR)
data(decathlon)
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
create.rmd(res.pca, file = "PCA.Rmd", document = "pdf_document")
readRmd(file = "PCA.Rmd", document = "pdf_document")

data(children)
res.ca = CA(children, row.sup = 15:18, col.sup = 6:8, graph = FALSE)
create.rmd(res.ca, file = "CA.Rmd", document = "html_document")
readRmd(file = "CA.Rmd", document = "html_document")

data(tea)
res.mca = MCA(tea, quanti.sup = 19,quali.sup = 20:36, graph = FALSE)
create.rmd(res.mca, file = "MCA.Rmd", document = c("word_document", "pdf_document"))
readRmd(file = "MCA.Rmd", document = "word_document")
}
}

\keyword{Rmarkdown}
