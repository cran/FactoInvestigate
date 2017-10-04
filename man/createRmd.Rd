\name{createRmd}
\alias{createRmd}

\title{
Create Rmarkdown file
}

\description{
Intialise a Rmarkdown file in which to write the results of the package functions.
}

\usage{
createRmd(res, file = "", 
          document = c("word_document", "pdf_document", "html_document"))
}

\arguments{
  \item{res}{
an object of class PCA, CA or MCA.
}
  \item{file}{
the file path where to write the description in Rmarkdown langage. If the file already exists, its content is overwritten. If not specified, the description is written in the console.
}
  \item{document}{
a character vector giving the document format desired between "word_document", "pdf_document" and "html_document".
}
}


\author{
Simon Thuleau and Francois Husson
}

\seealso{
  \code{\link{writeRmd}}, \code{\link{readRmd}}
}

\examples{
\dontrun{
data(decathlon)
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
createRmd(res.pca, file = "PCA.Rmd", document = "pdf_document")

data(children)
res.ca = CA(children, row.sup = 15:18, col.sup = 6:8, graph = FALSE)
createRmd(res.ca, file = "CA.Rmd", document = "html_document")

data(tea)
res.mca = MCA(tea, quanti.sup = 19,quali.sup = 20:36, graph = FALSE)
createRmd(res.mca, file = "MCA.Rmd", document = c("word_document", "pdf_document"))
}
}

\keyword{Rmarkdown}