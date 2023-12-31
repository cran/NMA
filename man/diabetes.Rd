% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.r
\docType{data}
\name{diabetes}
\alias{diabetes}
\title{Elliott and Mayer (2007)'s network meta-analysis data}
\format{An arm-based dataset with 22 studies}
\usage{
data(diabetes)
}
\description{
A network meta-analysis data from Elliott and Mayer (2007) that compared 5 antihypertensive drug classes and placebo for occurrence of diabetes.
\itemize{
  \item \code{study}: Study ID
  \item \code{trt}: Treatment (Diuretic, ACEI (ACE inhibitor), ARB, Beta blocker, CCB (Calcium-channel blocker), Placebo)
  \item \code{n}: Sample size
  \item \code{d}: Number of events (occurrence of diabetes)
}
}
\references{
Elliott, W. J., and Meyer, P. M. (2007). 
Incident diabetes in clinical trials of antihypertensive drugs: a network meta-analysis.
\emph{Lancet.}
\strong{369}: 201-207.
}
\keyword{datasets}
