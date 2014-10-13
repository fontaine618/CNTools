\name{dist-methods}
\docType{methods}
\alias{getDist}
\alias{getDist-method}
\alias{getDist,ANY-method}
\alias{getDist,RS-method}
\title{A method that extends the generic function dist to handle reduced segemnts}
\description{
  \code{getDist} gets the data contained by an RS object ready to be passed to the generic functon for distance calculations
}
\section{Methods}{
\describe{

\item{x = "ANY"}{see dist of stats}

\item{x = "RS"}{takes an RS object and then call the dist function of stats for distance calculations}
}}
\keyword{methods}
