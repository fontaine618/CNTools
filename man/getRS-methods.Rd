\name{getRS-methods}
\docType{methods}
\alias{getRS}
\alias{getRS-method}
\alias{getRS,CNSeg-method}
\alias{geneInfo}
\title{method that convert segment data into reduced segment matrix}
\description{
   \code{getRS} takes a CNSeg object containing the output of the segment function of DNAcopy and format the data into a matrix based on overlapping chromosome region (by  = region", gene (by = gene) or pair overlapping chromosome region (by = pair)
}
\section{Methods}{
\describe{

\item{object = "CNSeg"}{a reduced segment can be generated in three ways; by chromosomal regions that overlap across sample (by = region), by genes (by = gene), or by pair of samples with chromosome regions aligned (by = pair). User may choose to imput cells (by region or gene only) where a value can not be assigned by setting imput = TRUE. The X and Y chromosomes can dropped by stting XY = FALSE. }
}}
\keyword{methods}

