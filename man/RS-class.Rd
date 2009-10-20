\name{RS-class}
\docType{class}
\alias{RS-class}
\alias{rs<-}
\alias{rs}
\alias{segBy}
\alias{rs,RS-method}
\alias{rs<-,RS-method}
\alias{segBy,RS-method}
\alias{show,RS-method}

\title{Class "RS" contains the Reduced Segment data matrix derived
    from the output of segmentation}
\description{The class contains the Reduced Segment data matrix 
    derived from the output of the segment method of DNACopy and
    provides functions to manipulate the data or perform other 
    computational operations}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RS", ...)}.
A constructor RS(rs, by, imput, XY) can also be used to intantiate 
an object of this class
}
\section{Slots}{
	 \describe{
    \item{\code{rs}:}{Object of class \code{"ANY"} eith a matrix if
        the RS is by region or gene or a list of matrix if the
        RS is by sample pairs}
    \item{\code{by}:}{Object of class \code{"character"}  a 
        charactere string indicating how the RS matrix is obtained. 
        Valid values include region, gene, or pair}
  }
}
\section{Methods}{
  \describe{
    \item{cor}{\code{signature(x = "RS")}: Calculates the sample wise
        correlation coefficients using the Reduced Segment matix}
    \item{dist}{\code{signature(x = "RS")}: Calculates the distance
        between samples contained in the Reduced Segment matrix}
    \item{genefilter}{\code{signature(expr = "RS")}: filters features 
        by calling the genefilter function}
    \item{madFilter}{\code{signature(object = "RS")}: filters features
        by mean absolute deviation}
    \item{rs}{\code{signature(object = "RS")}: get method for the by 
        slot of an RS object}
    \item{rs<-}{\code{signature(object = "RS")}: assignment method for 
        the rs slot}
    \item{segBy}{\code{signature(object = "RS")}: get method for the by
        slot of an RS object}
    \item{show}{\code{signature(object = "RS")}: show method for an RS 
        object}
	 }
}
\author{Jianhua Zhang}

\seealso{
	\code{\linkS4class{CNSeg}}
}
\examples{
    data("sampleData")
    # take a subset of the data for speed
    seg <- CNSeg(sampleData[which(is.element(sampleData[, "ID"], sample(unique(sampleData[, "ID"]), 20))), ])
    rsByregion <- getRS(seg, by = "region", imput = TRUE, XY = FALSE)
    rsByregion
}
\keyword{classes}
