\name{Investigate}
\alias{Investigate}

\title{Resume factorial Analysis}

\description{Compute all the package functions : detection of outliers, evaluation of inertia distribution, dimensions description, classification and realisation of graphical views. All the results are written as Word, html or PDF documents.}

\usage{
Investigate(res, file = "Investigate.Rmd", document = c("html_document"), 
            Iselec = "contrib", Vselec = "cos2", Rselec = "contrib", 
            Cselec = "cos2", Mselec = "cos2", Icoef = 1, Vcoef = 1, Rcoef = 1, 
            Ccoef = 1, Mcoef = 1, ncp = NULL, time = "10s", nclust = -1, 
            mmax = 10, nmax = 10, hab = NULL, ellipse = TRUE, display.HCPC = TRUE, 
            out.selec = TRUE, remove.temp = TRUE, parallel = TRUE, cex = 0.7,
			openFile = TRUE, keepRmd = FALSE, codeGraphInd = NULL, 
			codeGraphVar=NULL, codeGraphCA = NULL, options = NULL, 
			language = "auto")
}

\arguments{
  \item{res}{
a PCA, CA or MCA object.
}
  \item{file}{
the file path where to write the description in Rmarkdown language. If the file already exists, its content is overwritten. If not specified, the description is written in the console.
}
  \item{document}{
a character vector giving the document format desired between "word_document", "pdf_document" and "html_document".
}
  \item{Iselec}{
the individuals to select ; see the details section.
}
  \item{Vselec}{
the variables to select ; see the details section.
}
  \item{Rselec}{
the rows to select (for a \code{\link{CA} res} object) ; see the details section.
}
  \item{Cselec}{
the columns to select (for a \code{\link{CA} res} object) ; see the details section.
}
  \item{Mselec}{
the supplementary variables to select ; see the details section.
}
  \item{Icoef}{
a numerical coefficient to adjust the individuals selection rule ; see the details section.
}
  \item{Vcoef}{
a numerical coefficient to adjust the variables selection rule ; see the details section.
}
  \item{Rcoef}{
a numerical coefficient to adjust the rows selection rule (for a \code{\link{CA} res} object) ; see the details section.
}
  \item{Ccoef}{
a numerical coefficient to adjust the columns selection rule (for a \code{\link{CA} res} object) ; see the details section.
}
  \item{Mcoef}{
a numerical coefficient to adjust the supplementary variables selection rule ; see the details section.
}
  \item{ncp}{
an integer to force the number of dimension to analyse.
}
  \item{time}{
a character indicating the loop condition. This string is made of a number and a letter coupled. The number X with letter \code{L} means to compute X datasets exactly. The number X with letter \code{s} means to compute as many datasets as possible during approximativley X seconds.
}
  \item{nclust}{
an integer to force the number of cluster for the classification.
}
  \item{mmax}{
an integer giving the maximum number of individuals (or rows) to illustrate each group (by defaut 10).
}
  \item{nmax}{
an integer giving the maximum number of variables (or columns) to illustrate each group of individuals (by defaut 10).
}
  \item{hab}{
a variable name or index to use to color the individuals (or rows) among the variable categories.
}
  \item{ellipse}{
a boolean : if \code{TRUE}, ellipses are plotted with the coloration of individuals (or rows).
}
  \item{display.HCPC}{
a boolean : if \code{TRUE}, the function performs the classification.
}
  \item{out.selec}{
a boolean : if \code{TRUE}, the function performs the detection of outliers.
}
  \item{remove.temp}{
a boolean : if \code{TRUE}, the temporary files created are deleted after the function execution.
}
  \item{parallel}{
a boolean : if \code{TRUE}, the computation uses map reduce on the processor cores to increase the performance. Useful for huge datasets.
}
  \item{cex}{
an optional argument for the generic plot functions, used to adjust the size of the elements plotted.
}
\item{openFile}{Open the file with the appropriate application; TRUE by default}
\item{keepRmd}{Keep the Rmd file; FALSE by default}
  \item{codeGraphInd}{
a character string corresponding to the code to use for the individuals graph.
}
  \item{codeGraphVar}{
a character string corresponding to the code to use for the variables graph.
}
  \item{codeGraphCA}{
a character string corresponding to the code to use for the CA graph.
}
\item{options}{a character string that gives the output options fir the figures. 
If NULL, \code{options="r, echo = FALSE, fig.align = 'center', fig.height = 3.5, fig.width = 5.5"} for linux and Mac
and \code{options="r, echo = FALSE, fig.height = 3.5, fig.width = 5.5"} for Windows}
\item{language}{possible values "auto", "en", or "fr": by default, "auto" detects the language (English or French), 
"en" for English and "fr" for "French"}
}

\details{
The \code{Iselec} argument (respectively \code{Vselec}, \code{Rselec} or \code{Cselec}) is used in order to select a part of the elements that are drawn and described. For example, you can use either :
\cr - \code{Iselec = 1:5} then the individuals (respectively the variables, the rows or the columns) numbered 1 to 5 are drawn.
\cr - \code{Iselec = c("name1","name5")} then the individuals (respectively the variables, the rows or the columns) named \code{name1} and \code{name5} are drawn.
\cr - \code{Iselec = "contrib 10"} then the 10 active or illustrative individuals (respectively the variables, the rows or the columns) that have the highest contribution on the 2 dimensions of the plane are drawn.
\cr - \code{Iselec = "contrib"} then the optimal number of active or illustrative individuals (respectively the variables, the rows or the columns) that have the highest contribution on the 2 dimensions of the plane are drawn.
\cr - \code{Iselec = "cos2 5"} then the 5 active or illustrative individuals (respectively the variables, the rows or the columns) that have the highest cos2 on the 2 dimensions of the plane are drawn.
\cr - \code{Iselec = "cos2 0.8"} then the active or illustrative individuals (respectively the variables, the rows or the columns) that have a \code{cos2} higher to \code{0.8} on the plane are drawn.
\cr - \code{Iselec = "cos2"} then the optimal number of active or illustrative individuals (respectively the variables, the rows or the columns) that have the highest cos2 on the 2 dimensions of the plane are drawn.
\cr 
\cr The \code{Icoef} argument (respectively \code{Vcoef}, \code{Rcoef} or \code{Ccoef}) is used in order to adjust the selection of the elements when based on \code{Iselec = "contrib"} or \code{Iselec = "cos2"}. For example :
\cr - if \code{Icoef = 2}, the threshold is 2 times higher, and thus 2 times more restrictive.
\cr - if \code{Icoef = 0.5}, the threshold is 2 times lower, and thus 2 times less restrictive.
}

\value{the function creates and opens a Word, html or PDF document that contains all the descriptions of analysis.}

\author{Simon Thuleau and Francois Husson}

\examples{
require(FactoMineR)
data(decathlon)
\dontrun{
res.pca = PCA(decathlon, quanti.sup = c(11:12), quali.sup = c(13), graph = FALSE)
Investigate(res.pca, file = "PCA.Rmd", document = "html_document", time = "1000L", 
            parallel = FALSE)

data(children)
res.ca = CA(children, row.sup = 15:18, col.sup = 6:8, graph = FALSE)
Investigate(res.ca, file = "CA.Rmd", document = "pdf_document")

data(tea)
res.mca = MCA(tea, quanti.sup = 19,quali.sup = 20:36, graph = FALSE)
Investigate(res.mca, file = "MCA.Rmd", document = c("word_document", "pdf_document"))
}
}

\keyword{multivariate analysis}
\keyword{report}
