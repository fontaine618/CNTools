\name{CNSeg-class}
\docType{class}
\alias{CNSeg-class}
\alias{segList}
\alias{sampleData}
\alias{CNSeg}
\alias{segList<-}
\alias{segList-method}
\alias{segList,CNSeg-method}
\alias{segList<-,CNSeg-method}
\alias{show,CNSeg-method}

\title{Class "CNSeg" contains the output of DNACopy segmentation data
   that can be operated on by the associated methods}
\description{Segmentation results of copy number data can not be 
    operated by other analysis tools such as classification. The 
    class "CNSeg" provides methods to convert segmentation data
    into a (what is called Reduced segments (RS))matrix format 
    so that other computations can be performed}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CNSeg",..)}.
A constructor \code{CNSeg(segList)} can also be used to instantiate 
CNSeg object
}
\section{Slots}{
	 \describe{
    \item{\code{segList}:}{Object of class \code{"data.frame"} is
    the "output" element of the object return by the segment 
    function of DNACopy}
  }
}
\section{Methods}{
  \describe{
    \item{getRS}{\code{signature(object = "CNSeg")}: get the 
        Reduced Segment data matrix}
    \item{segList}{\code{signature(object = "CNSeg")}: get the 
        value for slot "segList" of a CNSeg object}
    \item{segList<-}{\code{signature(object = "CNSeg")}: set the 
        value for slot "segList" of a CNSeg object}
    \item{show}{\code{signature(object = "CNSeg")}: print a CNSeg
        object}
	 }
}
\references{NA}
\author{Jianhua Zhang}
\seealso{
	\code{\linkS4class{RS}}
}
\examples{
   data("sampleData")
   # take a subset of the data for speed
   seg <- CNSeg(sampleData[which(is.element(sampleData[, "ID"], sample(unique(sampleData[, "ID"]), 10))), ])
    rsBypair <- getRS(seg, by = "pair", input = FALSE, XY = FALSE)
    rsBypair
}
\keyword{classes}
