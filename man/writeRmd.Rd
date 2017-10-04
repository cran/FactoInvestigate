\name{writeRmd}
\alias{writeRmd}

\title{
Write Rmarkdown file
}

\description{
Writes text or dumps a variable in a Rmarkdown file, and declares the utilisation and the configuration of a chunk.
}

\usage{
writeRmd(..., file = "", append = TRUE, sep = " ", end = "\n", dump = FALSE, 
         start = FALSE, stop = FALSE, options = NULL)
}

\arguments{
  \item{\dots}{
some R objects or other arguments to pass to the \code{\link{cat}} function.
}
  \item{file}{
the file path where to write the description in Rmarkdown langage. If not specified, the description is written in the console.
}
  \item{append}{
a boolean, if TRUE the text is written at the end of the file. Else it is overwritten.
}
  \item{sep}{
a character chain to insert between each element written in the file  (by default a blank space).
}
  \item{end}{
a character chain to add at the end of the text written in the file (by default a line break).
}
  \item{dump}{
a boolean : if \code{TRUE}, the text send to the function is interpreted as a variable name. A \code{dump} as to be written in a chunck declaration.
}
  \item{start}{
a boolean : if \code{TRUE}, the text written is preceded by a beginning chunk declaration.
}
  \item{stop}{
a boolean : if \code{TRUE}, the text written is preceded by a ending chunk declaration.
}
  \item{options}{
a character chain listing the options to declare for a chunk declaration.
}
}

\details{
To learn about all the possible chunck options, see \url{http://yihui.name/knitr/options}.
\cr Anyway, to declare a R langage chunk, write at least \code{"r"} as option.
}


\author{
Simon Thuleau and Francois Husson
}

\seealso{
  \code{\link{createRmd}}, \code{\link{readRmd}}
}

\examples{
\dontrun{
require(FactoMineR)
data(decathlon)
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
create.rmd(res.pca, file = "PCA.Rmd", document = "pdf_document")

drawn = selection(res.pca)$drawn

writeRmd(start = TRUE, options = "r, echo = FALSE, fig.align = 'center', fig.height = 3.5, 
         fig.width = 5.5", file = "PCA.Rmd", end = "")
writeRmd("drawn", file = file, dump = TRUE)
writeRmd("plot.PCA(res, select = drawn, choix = 'ind', invisible = 'quali', title = '')", 
         stop = TRUE, file = "PCA.Rmd")
           
writeRmd("**", figure.title, " - ", "Individuals factor map (PCA)", "**", file = "PCA.Rmd", 
         sep = "")
}
}
\keyword{Rmarkdown}