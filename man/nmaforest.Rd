\name{nmaforest}
\alias{nmaforest}
\title{Generating a ranked forest plot for the synthesis results of network meta-analysis}
\usage{
nmaforest(x,method="NH",col.plot="black",digits=3,ascending=TRUE)
}
\arguments{
\item{x}{Output object of \code{setup}}
\item{method}{Estimation and prediction method. \code{NH}: Noma-Hamura's improved REML-based methods (default). \code{REML}: The ordinary REML method. \code{fixed}: Fixed-effect model.}
\item{col.plot}{Color of the confidence interval plot (default: black)}
\item{digits}{Number of decimal places}
\item{ascending}{Type of order. Default is ascending order, but it can be changed to descending order changing to \code{FALSE}.}
}
\value{
A ranked forest plot for the synthesis results of network meta-analysis is generated.
}
\description{
A ranked forest plot for the synthesis results of network meta-analysis is generated based on the \code{forestplot} package by simple command. Details of the forestplot is customized by using the output objects of \code{obj.forest} function); see also \code{help(obj.forest)}. 
}
\examples{
data(heartfailure)

hf2 <- setup(study=study,trt=trt,d=d,n=n,measure="OR",ref="Placebo",data=heartfailure)

nmaforest(hf2)                           # Default setting
nmaforest(hf2, col.plot="blue")          # Change the color
nmaforest(hf2, ascending=FALSE)          # Change to the descending order
}
